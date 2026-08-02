#!/usr/bin/env python3
"""
Extract open-beta-only lipid tail combinations (not carried into the 2025
refined M3-Lipid-Parameters release) from the old monolithic phospholipids
bundle, before that bundle gets replaced. See conversation log for how this
target list and extraction logic were derived/validated.
"""
import re

missing = {"DIPA","DIPC","DIPE","DIPG","DIPS","JFPG","JPPG","LPPA","LPPC","LPPE",
        "LPPG","LPPS","OPPG","PGPA","PGPC","PGPE","PGPG","PGPS","PRPA","PRPC",
        "PRPE","PRPG","PRPS","PUPA","PUPC","PUPE","PUPS"}

src = "/usr/local/share/gromacs/top/martini3.ff/martini_v3.0.0_phospholipids_v1.itp"
with open(src) as f:
    text = f.read()

# tolerate whitespace inside brackets: [moleculetype] or [ moleculetype ]
parts = re.split(r'(?=^\[\s*moleculetype\s*\]\s*$)', text, flags=re.M)
blocks = [p for p in parts if re.match(r'^\[\s*moleculetype\s*\]', p.strip())]

molname_re = re.compile(r'^\s*([A-Za-z0-9]{2,8})\s+\d+\s*$', re.M)
found = {}
for b in blocks:
    body = b.split("\n", 1)[1] if "\n" in b else ""
    m = molname_re.search(body)
    name = m.group(1) if m else "UNRESOLVED"
    found.setdefault(name, []).append(b)

target_found = {k: v for k, v in found.items() if k in missing}
still_missing = missing - target_found.keys()
if still_missing:
    raise SystemExit(f"FATAL: extraction incomplete, missing: {sorted(still_missing)}")
if any(len(v) > 1 for v in found.values()):
    raise SystemExit(f"FATAL: duplicate moleculetype names detected during extraction")

out = "/usr/local/share/gromacs/top/martini3.ff/martini_v3.0.0_phospholipids_openbeta_extra.itp"
with open(out, "w") as f:
    f.write("; Open-beta lipid tail combinations not carried into M3-Lipid-Parameters (2025)\n\n")
    for name, blist in target_found.items():
        for b in blist:
            f.write(b)
print(f"Extracted {len(target_found)} open-beta-only lipids to {out}")