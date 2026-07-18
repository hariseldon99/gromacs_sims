#!/usr/bin/env python
import COBY

sysname = "ku04_gim"

COBY.COBY(
    box=[14, 14, 14],
    # Method A: using POPG and DOPE shortcuts to bypass internal diacyl mismatch bugs
    membrane="lipid:DVPE:46 lipid:POPG:19 lipid:POPE:12 lipid:DOPE:8 lipid:TOCL:2 apl:0.61",
    protein="file:KU04AMP01_cg_shifted.pdb moleculetypes:KU04 center_protein:False", #Protein was pre-aligned, so do not center again.
    solvation="default",
    itp_input="file:top_for_COBY.itp",
    # Feed COBY the single cardiolipin coordinate file you just generated
    molecule_import="file:tocl_single.gro moleculetypes:TOCL",
    out_sys=sysname,
    out_top=sysname + ".top",
    out_log=sysname + ".log",
    sn=sysname,
)
