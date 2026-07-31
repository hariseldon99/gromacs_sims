#!/usr/bin/env python
import COBY

# 3. Build system in COBY
sysname = "gim_only"

COBY.COBY(
    box=[14, 14, 12],
    membrane="lipid:DVPE:46 lipid:POPG:19 lipid:POPE:25 lipid:DOPE:8 lipid:TOCL:2 apl:0.61",
    solvation="default",
    itp_input="file:top_for_COBY.itp",
    molecule_import="file:tocl_single.gro moleculetypes:TOCL",
    out_sys=sysname,
    out_top=sysname + ".top",
    out_log=sysname + ".log",
    sn=sysname,
)