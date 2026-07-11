#!/usr/bin/env python3
"""
gmx_logparser.py — Parse GROMACS MD .log files for summary (completed or running) and ETA.

Features:
- Summary:
  - Performance (ns/day) [last and mean]
  - Simulated time from nsteps × dt and from last checkpoint
  - Walltime from Started/Finished or checkpoint timestamps
  - Averages (temperature, pressure, density, energies) if present in the log
  - Notes / warnings / fatal errors
- ETA (optional):
  - Uses recent checkpoint entries to estimate time remaining
  - Reports ns/day from recent progress
  - Works for running logs; can refresh in a loop with --watch

Usage:
  python gmx_logparser.py md.log
  python gmx_logparser.py md.log --eta
  python gmx_logparser.py md.log --watch 60 --eta
  python gmx_logparser.py md.log --json
"""
import re
import json
import argparse
import unicodedata
from datetime import datetime, timedelta
from tabulate import tabulate
from statistics import mean, stdev

MONTH_MAP = {
    'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6,
    'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12
}

def normalize_ts(ts: str) -> str:
    ts_norm = unicodedata.normalize('NFKC', ts)
    ts_norm = re.sub(r'\s+', ' ', ts_norm).strip()
    ts_norm = re.sub(r'\bS?Sep\b', 'Sep', ts_norm)  # tolerate odd "S Sep"
    return ts_norm

def ts_to_epoch(ts: str):
    try:
        parts = normalize_ts(ts).split()
        if len(parts) != 5:
            return None
        mon = MONTH_MAP.get(parts[1])
        if mon is None:
            return None
        day = int(parts[2])
        h, m, s = map(int, parts[3].split(':'))
        year = int(parts[4])
        return int(datetime(year, mon, day, h, m, s).timestamp())
    except Exception:
        return None

def seconds_to_hms(sec: float) -> str:
    return str(timedelta(seconds=int(sec)))

def _parse_energy_sections(lines: list[str]) -> list[dict]:
    """
    Parse 'Energies (kJ/mol)' blocks printed in GROMACS logs.
    For each header row, align headers to the following numeric row by using
    the numeric tokens' start columns, so multi-word labels like
    'Pres. DC (bar)' and 'Pressure (bar)' stay separate.
    Returns a list of dicts (one per block).
    """
    num_re = re.compile(r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)')
    blocks: list[dict] = []
    i = 0
    while i < len(lines):
        if re.search(r'\bEnergies\s*\(kJ/mol\)', lines[i], re.IGNORECASE):
            i += 1
            block: dict = {}
            while i < len(lines):
                hdr = lines[i].rstrip('\n')
                if not hdr.strip():
                    i += 1
                    break  # blank line usually ends the energies section
                if re.search(r'\bEnergies\s*\(kJ/mol\)', hdr, re.IGNORECASE):
                    break  # next section
                has_letters = re.search(r'[A-Za-z]', hdr) is not None
                if not has_letters:
                    i += 1
                    continue

                # Find the values line after the header
                j = i + 1
                while j < len(lines) and not lines[j].strip():
                    j += 1
                if j >= len(lines):
                    i = j
                    break
                vals_line = lines[j].rstrip('\n')

                # Extract numeric tokens with positions
                nums = [(m.start(1), m.end(1), float(m.group(1))) for m in num_re.finditer(vals_line)]
                if not nums:
                    i = j + 1
                    continue

                # Build header segments aligned to numeric columns
                header_map: dict = {}
                for idx, (start, end, val) in enumerate(nums):
                    seg_start = start
                    seg_end = nums[idx + 1][0] if idx + 1 < len(nums) else len(vals_line)
                    # Clip to header line length
                    if seg_start >= len(hdr):
                        label = ''
                    else:
                        label = hdr[seg_start:min(seg_end, len(hdr))].strip()
                    if not label:
                        # Fallback: try splitting header on 2+ spaces and index
                        alt_labels = [h.strip() for h in re.split(r'\s{2,}', hdr.strip()) if h.strip()]
                        if idx < len(alt_labels):
                            label = alt_labels[idx]
                        else:
                            label = f'col{idx+1}'
                    header_map[label] = val

                # Merge into current block and advance
                block.update(header_map)
                i = j + 1
                continue

            if block:
                blocks.append(block)
            continue
        i += 1
    return blocks

def _find_last(lines: list[str], patterns: list[str], conv=float):
    """Find the last occurrence matching any pattern; return converted group(1)."""
    for ln in reversed(lines):
        for pat in patterns:
            m = re.search(pat, ln, re.IGNORECASE)
            if m:
                try:
                    return conv(m.group(1))
                except Exception:
                    continue
    return None

def _parse_last_tpd(lines: list[str]) -> dict:
    """
    Parse the last 'Temperature Pressure [Density]' mini-table in the log.
    Returns a dict with temperature_K, pressure_bar, density_kg_per_m3 when found.
    """
    last = {}
    # Two header variants: full and abbreviated
    pat_full = re.compile(r'^\s*Temperature\s+Pressure(?:\s+Density)?\s*$', re.IGNORECASE)
    pat_abbr = re.compile(r'^\s*Temp\w*\s+Pres\w*(?:\s+Dens\w*)?\s*$', re.IGNORECASE)
    i = 0
    while i < len(lines):
        ln = lines[i].rstrip()
        if pat_full.match(ln) or pat_abbr.match(ln):
            # Find next non-empty line with numbers
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines):
                nums = re.findall(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', lines[j])
                vals = [float(x) for x in nums] if nums else []
                if len(vals) >= 2:
                    last = {
                        "temperature_K": vals[0],
                        "pressure_bar": vals[1],
                    }
                    if len(vals) >= 3:
                        last["density_kg_per_m3"] = vals[2]
                # continue scanning; we want the last occurrence
                i = j
                continue
        i += 1
    return last

def _parse_all_tpd(lines: list[str]) -> list[tuple[float | None, float | None, float | None]]:
    """
    Return a list of (temperature_K, pressure_bar, density_kg_per_m3) tuples
    found in 'Temperature Pressure [Density]' mini-tables throughout the log.
    """
    out: list[tuple[float | None, float | None, float | None]] = []
    pat_full = re.compile(r'^\s*Temperature\s+Pressure(?:\s+Density)?\s*$', re.IGNORECASE)
    pat_abbr = re.compile(r'^\s*Temp\w*\s+Pres\w*(?:\s+Dens\w*)?\s*$', re.IGNORECASE)
    i = 0
    while i < len(lines):
        ln = lines[i].rstrip()
        if pat_full.match(ln) or pat_abbr.match(ln):
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines):
                nums = re.findall(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', lines[j])
                vals = [float(x) for x in nums] if nums else []
                t = vals[0] if len(vals) >= 1 else None
                p = vals[1] if len(vals) >= 2 else None
                d = vals[2] if len(vals) >= 3 else None
                out.append((t, p, d))
                i = j
                continue
        i += 1
    return out

def _collect_tp_series(lines: list[str]) -> tuple[list[float], list[float]]:
    """
    Collect temperature (K) and pressure (bar) samples from:
    - Energies (kJ/mol) blocks (keys: 'Temperature', 'Pressure (bar)', avoid 'Pres. DC')
    - Temperature/Pressure mini-tables
    """
    temps: list[float] = []
    press: list[float] = []
    # From energies blocks
    for blk in _parse_energy_sections(lines):
        # Temperature
        t_key = next((k for k in blk.keys() if k.lower().startswith('temperature')), None)
        if t_key is not None:
            try:
                temps.append(float(blk[t_key]))
            except Exception:
                pass
        # Pressure (prefer 'Pressure (bar)', exclude 'Pres. DC')
        p_val = None
        if 'Pressure (bar)' in blk:
            p_val = blk.get('Pressure (bar)')
        else:
            for k in blk.keys():
                kl = k.lower()
                if 'pressure' in kl and not kl.startswith('pres. dc'):
                    p_val = blk[k]
                    break
        if p_val is not None:
            try:
                press.append(float(p_val))
            except Exception:
                pass
    # From mini-tables
    for t, p, _ in _parse_all_tpd(lines):
        if t is not None:
            temps.append(t)
        if p is not None:
            press.append(p)
    return temps, press

def _fmt_val_sd(val: float | None, sd: float | None, digits: int = 2) -> str:
    if val is None:
        return "NA"
    if sd is None or sd == 0:
        return f"{val:.{digits}f}"
    return f"{val:.{digits}f} ± {sd:.{digits}f}"

def parse_summary(logfile: str) -> dict:
    with open(logfile, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    # Extract dt (ps) and nsteps
    dt_ps = None
    nsteps = None
    for ln in lines:
        if 'dt' in ln:
            m = re.search(r'\bdt\s*=\s*([0-9.]+)', ln)
            if m:
                try:
                    dt_ps = float(m.group(1))
                except ValueError:
                    pass
        if 'nsteps' in ln.lower():
            m = re.search(r'nsteps\s*=\s*(\d+)', ln, re.IGNORECASE)
            if m:
                try:
                    nsteps = int(m.group(1))
                except ValueError:
                    pass
        if dt_ps is not None and nsteps is not None:
            break
    if nsteps is None:
        for ln in lines:
            m = re.search(r'(\d+)\s+steps', ln, re.IGNORECASE)
            if m:
                nsteps = int(m.group(1))
                break

    # Performance lines
    perf_vals = []
    for ln in lines:
        m = re.search(r'Performance:\s*([0-9.]+)\s*(?:ns/day|ns\/day)', ln, re.IGNORECASE)
        if m:
            try:
                perf_vals.append(float(m.group(1)))
            except ValueError:
                pass

    # Timestamps
    start_ts = None
    finish_ts = None
    checkpoint_ts = []
    for ln in lines:
        s = ln.strip()
        m = re.search(r'(?:Started|Starting) mdrun .* at (.+)$', s, re.IGNORECASE)
        if m:
            start_ts = ts_to_epoch(m.group(1))
        m = re.search(r'(?:Finished|Stopping) mdrun .* at (.+)$', s, re.IGNORECASE)
        if m:
            finish_ts = ts_to_epoch(m.group(1))
        m = re.search(r'Writing checkpoint, step\s+[0-9,]+\s+at\s+(.+)$', s, re.IGNORECASE)
        if m:
            ts = ts_to_epoch(m.group(1))
            if ts:
                checkpoint_ts.append(ts)

    # Last step (from checkpoint lines)
    last_step = None
    for ln in reversed(lines):
        m = re.search(r'Writing checkpoint, step\s+([0-9,]+)', ln, re.IGNORECASE)
        if m:
            try:
                last_step = int(m.group(1).replace(',', ''))
                break
            except ValueError:
                pass

    # Averages (best-effort; dependent on log content)
    def find_avg(patterns, conv=float):
        for ln in lines:
            for pat in patterns:
                m = re.search(pat, ln, re.IGNORECASE)
                if m:
                    try:
                        return conv(m.group(1))
                    except ValueError:
                        continue
        return None

    avg_temp_K = find_avg([
        r'Average temperature[:=]\s*([0-9.]+)\s*K',
        r'Temperature[:=]\s*([0-9.]+)\s*K',
    ])
    avg_press_bar = find_avg([
        r'Average pressure[:=]\s*([0-9.\-]+)\s*bar',
        r'Pressure[:=]\s*([0-9.\-]+)\s*bar',
    ])
    avg_density = find_avg([
        r'Average density[:=]\s*([0-9.]+)\s*kg/m\^?3',
        r'Density[:=]\s*([0-9.]+)\s*kg/m\^?3',
    ])
    avg_potential_kjmol = find_avg([
        r'Average potential energy[:=]\s*([0-9.\-eE]+)\s*kJ/mol',
        r'Potential Energy[:=]\s*([0-9.\-eE]+)\s*kJ/mol',
        r'Potential\s*[:=]\s*([0-9.\-eE]+)\s*kJ/mol',
    ])
    avg_total_energy_kjmol = find_avg([
        r'Average total energy[:=]\s*([0-9.\-eE]+)\s*kJ/mol',
        r'Total Energy[:=]\s*([0-9.\-eE]+)\s*kJ/mol',
        r'Total\s*[:=]\s*([0-9.\-eE]+)\s*kJ/mol',
    ])

    # Notes / warnings / errors
    notes, warnings, errors = [], [], []
    for ln in lines:
        if re.search(r'^\s*NOTE\b', ln):
            notes.append(ln.strip())
        if re.search(r'^\s*WARNING\b', ln):
            warnings.append(ln.strip())
        if re.search(r'Fatal error', ln, re.IGNORECASE):
            errors.append(ln.strip())

    summary = {}
    if perf_vals:
        summary['performance_ns_per_day_last'] = perf_vals[-1]
        summary['performance_ns_per_day_mean'] = sum(perf_vals) / len(perf_vals)
    else:
        # No final "Performance:" line yet (likely still running). Estimate from checkpoints.
        try:
            total_steps, checkpoints, dt_ps_det = parse_checkpoints(logfile)
            if checkpoints and len(checkpoints) >= 2:
                (s1, t1), (s2, t2) = checkpoints[-2], checkpoints[-1]
                ts1, ts2 = ts_to_epoch(t1), ts_to_epoch(t2)
                if ts1 and ts2 and s1 is not None and s2 is not None and ts2 > ts1:
                    dt_ps_used = dt_ps_det if dt_ps_det is not None else 0.002
                    ns_per_step = dt_ps_used / 1000.0
                    ns_per_sec = ((s2 - s1) * ns_per_step) / (ts2 - ts1)
                    summary['performance_ns_per_day_est'] = ns_per_sec * 86400.0
        except Exception:
            pass

    if dt_ps is not None and nsteps is not None:
        summary['simulated_time_ns'] = (dt_ps * nsteps) / 1000.0
    if last_step is not None and dt_ps is not None:
        summary['simulated_time_ns_from_last_checkpoint'] = (dt_ps * last_step) / 1000.0

    summary['nsteps'] = nsteps
    summary['dt_ps'] = dt_ps

    wall_start = start_ts if start_ts else (min(checkpoint_ts) if checkpoint_ts else None)
    wall_end = finish_ts if finish_ts else (max(checkpoint_ts) if checkpoint_ts else None)
    if wall_start and wall_end and wall_end >= wall_start:
        wall_sec = wall_end - wall_start
        summary['walltime_seconds'] = wall_sec
        summary['walltime_hms'] = seconds_to_hms(wall_sec)

    averages = {}
    if avg_temp_K is not None:
        averages['temperature_K'] = avg_temp_K
    if avg_press_bar is not None:
        averages['pressure_bar'] = avg_press_bar
    if avg_density is not None:
        averages['density_kg_per_m3'] = avg_density
    if avg_potential_kjmol is not None:
        averages['potential_energy_kJ_per_mol'] = avg_potential_kjmol
    if avg_total_energy_kjmol is not None:
        averages['total_energy_kJ_per_mol'] = avg_total_energy_kjmol
    if averages:
        summary['averages'] = averages

    if notes:
        summary['notes'] = notes
    if warnings:
        summary['warnings'] = warnings
    if errors:
        summary['errors'] = errors

    # Parse energies (last block)
    energy_blocks = _parse_energy_sections(lines)
    if energy_blocks:
        summary['energies'] = energy_blocks[-1]
        summary['energies_terms'] = list(energy_blocks[-1].keys())

    # Estimated temperature/pressure:
    # - Prefer explicit averages parsed above.
    # - Else try last instantaneous reports in the log.
    # - Else try from the last energy block if keys exist.
    if avg_temp_K is None:
        temp_last = _find_last(lines, [
            r'\bTemperature[:=]\s*([0-9.]+)\s*K',
            r'\bTemp\.?\s*[:=]\s*([0-9.]+)\s*K',
        ], float)
        if temp_last is None and 'energies' in summary:
            for key in summary['energies'].keys():
                if key.lower().startswith('temperature'):
                    temp_last = summary['energies'][key]
                    break
        if temp_last is None:
            tpd = _parse_last_tpd(lines)
            if "temperature_K" in tpd:
                temp_last = tpd["temperature_K"]
        if temp_last is not None:
            summary['temperature_K_est'] = temp_last
    if avg_press_bar is None:
        press_last = _find_last(lines, [
            r'\bPressure[:=]\s*([0-9.\-]+)\s*bar',
            r'\bPres\.?\s*[:=]\s*([0-9.\-]+)\s*bar',
        ], float)
        if press_last is None and 'energies' in summary:
            for key in summary['energies'].keys():
                if key.lower().startswith('pressure'):
                    press_last = summary['energies'][key]
                    break
        if press_last is None:
            tpd = _parse_last_tpd(lines)
            if "pressure_bar" in tpd:
                press_last = tpd["pressure_bar"]
        if press_last is not None:
            summary['pressure_bar_est'] = press_last
    # Also capture density estimate if present
    if summary.get('averages', {}).get('density_kg_per_m3') is None:
        tpd = _parse_last_tpd(lines)
        if "density_kg_per_m3" in tpd:
            summary['density_kg_per_m3_est'] = tpd["density_kg_per_m3"]

    # Compute series mean and standard deviation for T/P from all occurrences
    t_series, p_series = _collect_tp_series(lines)
    if t_series:
        summary['temperature_series_n'] = len(t_series)
        summary['temperature_K_est_mean'] = mean(t_series)
        summary['temperature_K_est_sd'] = stdev(t_series) if len(t_series) > 1 else None
        # Backfill single-value estimate if missing
        summary.setdefault('temperature_K_est', t_series[-1])
    if p_series:
        summary['pressure_series_n'] = len(p_series)
        summary['pressure_bar_est_mean'] = mean(p_series)
        summary['pressure_bar_est_sd'] = stdev(p_series) if len(p_series) > 1 else None
        summary.setdefault('pressure_bar_est', p_series[-1])

    return summary

def parse_checkpoints(logfile: str):
    """Return (total_steps, checkpoints[(step, ts_str)], dt_ps_detected)."""
    with open(logfile, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()

    total_steps = None
    for line in lines:
        if "nsteps" in line.lower():
            m = re.search(r'nsteps\s*=\s*(\d+)', line, re.IGNORECASE)
            if m:
                total_steps = int(m.group(1))
                break
    if total_steps is None:
        for line in lines:
            m = re.search(r'(\d+)\s+steps', line, re.IGNORECASE)
            if m:
                total_steps = int(m.group(1))
                break

    dt_ps = None
    for line in lines:
        if re.search(r'\bdt\b', line):
            m = re.search(r'dt\s*=\s*([0-9.]+)', line)
            if m:
                try:
                    dt_ps = float(m.group(1))
                    break
                except ValueError:
                    pass

    checkpoints = []
    for line in lines:
        m = re.search(r'Writing checkpoint, step\s+([0-9,]+)\s+at\s+(.+)', line)
        if m:
            step_str = m.group(1).replace(",", "").strip()
            try:
                step = int(step_str)
            except ValueError:
                step = None
            ts_clean = normalize_ts(m.group(2).strip())
            checkpoints.append((step, ts_clean))

    return total_steps, checkpoints, dt_ps

def estimate_eta(logfile: str, smooth_n: int = 2, dt_ps_arg: float = 0.002) -> dict:
    total_steps, checkpoints, dt_ps_detected = parse_checkpoints(logfile)
    if total_steps is None:
        raise ValueError("Could not determine total steps from log file.")
    if len(checkpoints) < 2:
        raise ValueError("Not enough checkpoint entries found.")

    dt_ps_used = dt_ps_detected if dt_ps_detected is not None else dt_ps_arg
    n = min(smooth_n, len(checkpoints))
    recent = checkpoints[-n:]

    total_time = 0
    total_steps_done = 0
    for (s1, t1), (s2, t2) in zip(recent[:-1], recent[1:]):
        ts1 = ts_to_epoch(t1)
        ts2 = ts_to_epoch(t2)
        if ts1 is None or ts2 is None or s1 is None or s2 is None:
            continue
        total_time += (ts2 - ts1)
        total_steps_done += (s2 - s1)

    if total_steps_done == 0:
        raise ValueError("No progress detected in the selected smoothing window.")

    time_per_step = total_time / total_steps_done
    current_step = recent[-1][0]
    steps_remaining = total_steps - current_step
    time_remaining = steps_remaining * time_per_step

    ns_per_step = dt_ps_used / 1000.0
    ns_per_sec = (total_steps_done * ns_per_step) / total_time
    ns_per_day = ns_per_sec * 86400

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return {
        "timestamp": now,
        "total_steps": total_steps,
        "current_step": current_step,
        "eta_hms": seconds_to_hms(time_remaining),
        "eta_seconds": time_remaining,
        "progress": current_step / total_steps if total_steps else None,
        "ns_per_day": ns_per_day,
        "dt_ps_used": dt_ps_used,
        "dt_ps_detected": dt_ps_detected
    }

def print_summary(summary: dict, tablefmt: str, logfile: str, eta: dict | None = None):
    rows: list[tuple[str, str | float | int]] = [
        ("Log file", logfile),
        ("Performance (last) ns/day", f"{summary.get('performance_ns_per_day_last', float('nan')):.2f}" if summary.get('performance_ns_per_day_last') is not None else "NA"),
        ("Performance (mean) ns/day", f"{summary.get('performance_ns_per_day_mean', float('nan')):.2f}" if summary.get('performance_ns_per_day_mean') is not None else "NA"),
        ("Performance (est.) ns/day", f"{summary.get('performance_ns_per_day_est', float('nan')):.2f}" if summary.get('performance_ns_per_day_est') is not None else "NA"),
        ("Simulated time (nsteps×dt) ns", f"{summary.get('simulated_time_ns', float('nan')):.3f}" if summary.get('simulated_time_ns') is not None else "NA"),
        ("Sim time (from last checkpoint) ns", f"{summary.get('simulated_time_ns_from_last_checkpoint', float('nan')):.3f}" if summary.get('simulated_time_ns_from_last_checkpoint') is not None else "NA"),
        ("nsteps", summary.get('nsteps', "NA")),
        ("dt (ps)", summary.get('dt_ps', "NA")),
        ("Walltime (h:m:s)", summary.get('walltime_hms', "NA")),
        ("Walltime (s)", summary.get('walltime_seconds', "NA")),
    ]
    # Averages appended
    if 'averages' in summary and summary['averages']:
        for k, v in summary['averages'].items():
            # Attach sd as error if available for temp/pressure
            if k == 'temperature_K':
                rows.append((f"Average {k}", _fmt_val_sd(v, summary.get('temperature_K_est_sd'))))
            elif k == 'pressure_bar':
                rows.append((f"Average {k}", _fmt_val_sd(v, summary.get('pressure_bar_est_sd'))))
            else:
                rows.append((f"Average {k}", v))

    # Estimated T/P/D with ± sd when available
    if summary.get('temperature_K_est_mean') is not None or summary.get('temperature_K_est') is not None:
        rows.append(("Estimated temperature (K)", _fmt_val_sd(summary.get('temperature_K_est_mean') or summary.get('temperature_K_est'),
                                                               summary.get('temperature_K_est_sd'))))
    if summary.get('pressure_bar_est_mean') is not None or summary.get('pressure_bar_est') is not None:
        rows.append(("Estimated pressure (bar)", _fmt_val_sd(summary.get('pressure_bar_est_mean') or summary.get('pressure_bar_est'),
                                                             summary.get('pressure_bar_est_sd'))))
    if summary.get('density_kg_per_m3_est') is not None:
        rows.append(("Estimated density (kg/m^3)", summary['density_kg_per_m3_est']))

    # Energies (last block)
    if 'energies' in summary and summary['energies']:
        for k, v in summary['energies'].items():
            rows.append((f"Energy: {k} (kJ/mol)" if 'temperature' not in k.lower() and 'pressure' not in k.lower() else f"{k}", v))

    # ETA rows (if provided)
    if eta:
        rows.extend([
            ("ETA Timestamp", eta["timestamp"]),
            ("Step (current/total)", f"{eta['current_step']} / {eta['total_steps']}"),
            ("Progress", f"{eta['progress']*100:.2f}%" if eta['progress'] is not None else "NA"),
            ("ETA (h:m:s)", eta["eta_hms"]),
            ("ETA (s)", int(eta["eta_seconds"])),
            ("Perf (ns/day, est.)", f"{eta['ns_per_day']:.2f}"),
            ("dt used (ps)", eta["dt_ps_used"]),
            ("dt detected (ps)", eta["dt_ps_detected"] if eta["dt_ps_detected"] is not None else "NA"),
        ])

    print(tabulate(rows, headers=["Metric", "Value"], tablefmt=tablefmt))

    def print_msgs(name: str, msgs):
        if not msgs:
            return
        data = [(i + 1, m) for i, m in enumerate(msgs)]
        print("\n" + tabulate(data, headers=[name, "Message"], tablefmt=tablefmt))

    print_msgs("Note #", summary.get('notes', []))
    print_msgs("Warning #", summary.get('warnings', []))
    print_msgs("Error #", summary.get('errors', []))

def main():
    ap = argparse.ArgumentParser(description="Parse GROMACS MD .log for summary and ETA.")
    ap.add_argument("logfile", help="Path to GROMACS .log file")
    ap.add_argument("--json", action="store_true", help="Print JSON summary/ETA")
    ap.add_argument("--tablefmt", default="fancy_grid", help="tabulate format: github, grid, simple, psql, plain, fancy_grid, etc.")
    ap.add_argument("--eta", action="store_true", help="Calculate ETA from checkpoint timestamps")
    ap.add_argument("--watch", type=int, help="Refresh ETA every N seconds (implies --eta)")
    ap.add_argument("--smooth", type=int, default=2, help="Average ETA over last N checkpoints (default: 2)")
    ap.add_argument("--dt", type=float, default=0.002, help="Timestep in ps if not found in log (default: 0.002)")
    ap.add_argument("--no_summary", action="store_true", help="Do not print summary table")
    args = ap.parse_args()

    # Always try to parse summary once (unless suppressed)
    summary = parse_summary(args.logfile)

    # Decide printing mode
    if args.json:
        out = {"summary": summary}
        if args.eta or args.watch:
            try:
                eta = estimate_eta(args.logfile, smooth_n=args.smooth, dt_ps_arg=args.dt)
                out["eta"] = eta
            except Exception as e:
                out["eta_error"] = str(e)
        print(json.dumps(out, indent=2))
        return

    if not args.no_summary and not (args.eta or args.watch):
        print_summary(summary, args.tablefmt, args.logfile)

    # ETA once or watch loop
    if args.eta or args.watch:
        import time
        interval = args.watch if args.watch else 0
        while True:
            try:
                # Re-parse summary to reflect latest log content while running
                summary = parse_summary(args.logfile)
                eta = estimate_eta(args.logfile, smooth_n=args.smooth, dt_ps_arg=args.dt)
                print_summary(summary, args.tablefmt, args.logfile, eta=eta)
            except Exception as e:
                print(f"\nETA error: {e}")

            if not interval:
                break
            time.sleep(interval)

if __name__ == "__main__":
    main()