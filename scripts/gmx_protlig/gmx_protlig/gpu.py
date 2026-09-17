"""
GPU acceleration hook for gromacs_py.

Dynamically intercepts GROMACS `mdrun` commands issued by gromacs_py via
os_command_py.os_command.Command, injecting hardware acceleration flags
(-pin, -pme, -bonded, -update) while distinguishing between energy minimization
(integrator=steep/cg) and dynamical MD (integrator=md/sd).
"""

import os
import logging
from typing import List, Optional

logger = logging.getLogger("gmx_protlig.gpu")

_ORIGINAL_COMMAND_INIT = None
_HOOK_STATE = {
    "enabled": False,
    "pin": "on",
    "pme": "gpu",
    "bonded": "gpu",
    "update": "gpu",
    "extra_flags": [],
    "verbose": True,
}


def _is_minimization_cmd(list_cmd: List[str]) -> bool:
    """
    Check whether an mdrun command invocation is running energy minimization.
    Minimization integrators (steep, cg, l-bfgs) do not support dynamical GPU features
    like -bonded gpu or -update gpu.
    """
    for i, arg in enumerate(list_cmd):
        if arg in ("-s", "-deffnm") and i + 1 < len(list_cmd):
            target = list_cmd[i + 1].lower()
            if any(k in target for k in ["em", "mini", "init_em"]):
                return True
            # Check matching .mdp file in current directory if exists
            mdp_file = f"{list_cmd[i + 1]}.mdp"
            if os.path.isfile(mdp_file):
                try:
                    with open(mdp_file, "r", encoding="utf-8", errors="replace") as fh:
                        for line in fh:
                            line = line.strip().lower()
                            if line.startswith("integrator"):
                                parts = line.split("=", 1)
                                if len(parts) == 2:
                                    alg = parts[1].split(";")[0].strip()
                                    if alg in ("steep", "cg", "l-bfgs", "nm"):
                                        return True
                except Exception:
                    pass
    return False


def _custom_command_init(self, list_cmd, *args, **kwargs):
    global _ORIGINAL_COMMAND_INIT, _HOOK_STATE

    if _HOOK_STATE.get("enabled") and len(list_cmd) >= 2 and list_cmd[1] == "mdrun":
        list_cmd = list(list_cmd)
        is_em = _is_minimization_cmd(list_cmd)

        flags_to_add = []
        # -pin is always safe and beneficial
        pin_val = _HOOK_STATE.get("pin")
        if pin_val and "-pin" not in list_cmd:
            flags_to_add.extend(["-pin", str(pin_val)])

        if is_em:
            if _HOOK_STATE.get("verbose"):
                logger.info(
                    f"[gromacs_py GPU Hook] Energy Minimization detected; "
                    f"using CPU bonded/update with flags: {' '.join(flags_to_add)}"
                )
        else:
            # Dynamical MD (equi/prod)
            pme_val = _HOOK_STATE.get("pme")
            bonded_val = _HOOK_STATE.get("bonded")
            update_val = _HOOK_STATE.get("update")

            if pme_val and "-pme" not in list_cmd:
                flags_to_add.extend(["-pme", str(pme_val)])
            if bonded_val and "-bonded" not in list_cmd:
                flags_to_add.extend(["-bonded", str(bonded_val)])
            if update_val and "-update" not in list_cmd:
                flags_to_add.extend(["-update", str(update_val)])

            extra = _HOOK_STATE.get("extra_flags") or []
            flags_to_add.extend(extra)

            if _HOOK_STATE.get("verbose"):
                logger.info(
                    f"[gromacs_py GPU Hook] Dynamical MD detected (equi/prod); "
                    f"injected GPU acceleration flags: {' '.join(flags_to_add)}"
                )

        list_cmd.extend(flags_to_add)

    _ORIGINAL_COMMAND_INIT(self, list_cmd, *args, **kwargs)


def enable_gromacs_py_gpu_hook(
    pin: str = "on",
    pme: str = "gpu",
    bonded: str = "gpu",
    update: str = "gpu",
    extra_flags: Optional[List[str]] = None,
    verbose: bool = True,
) -> bool:
    """
    Enable GPU acceleration flags injection into gromacs_py's mdrun calls.
    Returns True if successfully enabled or already enabled with updated settings.
    """
    global _ORIGINAL_COMMAND_INIT, _HOOK_STATE

    try:
        import os_command_py.os_command as osc
    except ImportError:
        logger.warning(
            "os_command_py is not available in current environment; "
            "cannot monkey-patch gromacs_py Command."
        )
        return False

    _HOOK_STATE["enabled"] = True
    _HOOK_STATE["pin"] = pin
    _HOOK_STATE["pme"] = pme
    _HOOK_STATE["bonded"] = bonded
    _HOOK_STATE["update"] = update
    _HOOK_STATE["extra_flags"] = list(extra_flags) if extra_flags else []
    _HOOK_STATE["verbose"] = verbose

    if _ORIGINAL_COMMAND_INIT is None:
        _ORIGINAL_COMMAND_INIT = osc.Command.__init__
        osc.Command.__init__ = _custom_command_init
        logger.info(
            f"gromacs_py GPU acceleration hook activated "
            f"(-pin {pin}, -pme {pme}, -bonded {bonded}, -update {update})"
        )

    return True


def disable_gromacs_py_gpu_hook() -> None:
    """Disable the GPU acceleration hook and restore original Command.__init__."""
    global _ORIGINAL_COMMAND_INIT, _HOOK_STATE
    _HOOK_STATE["enabled"] = False
    if _ORIGINAL_COMMAND_INIT is not None:
        try:
            import os_command_py.os_command as osc
            osc.Command.__init__ = _ORIGINAL_COMMAND_INIT
        except ImportError:
            pass
        _ORIGINAL_COMMAND_INIT = None
        logger.info("gromacs_py GPU acceleration hook disabled.")


def is_gpu_hook_enabled() -> bool:
    """Check if the GPU acceleration hook is currently active."""
    return bool(_HOOK_STATE.get("enabled"))
