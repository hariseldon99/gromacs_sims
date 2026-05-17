#! /usr/bin/env python3
import re

with open("membrane.gro") as f:
    lines = f.readlines()
title = lines[0]
natoms_mem = int(lines[1].strip())
mem_atoms = lines[2:2+natoms_mem]
box = lines[2+natoms_mem]

with open("pmb_placed.gro") as f:
    plines = f.readlines()
natoms_pmb = int(plines[1].strip())
pmb_atoms = plines[2:2+natoms_pmb]

with open("system.gro", "w") as out:
    out.write(title)
    out.write(f"{natoms_mem + natoms_pmb}\n")
    out.writelines(mem_atoms)
    out.writelines(pmb_atoms)
    out.write(box)

