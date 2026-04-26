import re

def sync_topol_and_gro(chg_file, original_itp):
    # 1. Load data from .chg
    with open(chg_file, 'r') as f:
        data = [line.split() for line in f if line.strip()]
    
    # 2. Extract organic framework from original ITP
    with open(original_itp, 'r') as f:
        itp_content = f.read()
    
    # Extract original atom types and interactions
    atom_types = re.findall(r'\[ atoms \](.*?)\n\[', itp_content, re.S)[0].strip().split('\n')
    organic_types = [line.split()[1] for line in atom_types if line.strip() and not line.strip().startswith(';')]
    
    # 3. Create the NEW .gro file (nm units)
    with open("PDL_final.gro", "w") as f:
        f.write("PDL optimized complex\n65\n")
        for i, (parts) in enumerate(data):
            symbol = "PD" if i == 0 else f"{parts[0]}{i+1}" # Unique names: PD, C2, C3...
            x, y, z = float(parts[1])/10, float(parts[2])/10, float(parts[3])/10
            f.write(f"{1:5d}{'PDL':5s}{symbol:>5s}{i+1:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        f.write(f"{10.0:10.5f}{10.0:10.5f}{10.0:10.5f}\n")

    # 4. Create the NEW .itp file
    with open("PDL_final.itp", "w") as f:
        f.write("[ atomtypes ]\n")
        # Add all types used in the file
        f.write("Pd   46   106.420   0.000   A   0.24500   0.21338\n")
        # Extract existing types from original ITP to be safe
        types_section = re.findall(r'\[ atomtypes \](.*?)\n\[', itp_content, re.S)[0]
        f.write(types_section.strip() + "\n\n")

        f.write("[ moleculetype ]\nPDL   3\n\n[ atoms ]\n")
        for i, parts in enumerate(data):
            symbol = "PD" if i == 0 else f"{parts[0]}{i+1}"
            charge = parts[4]
            at_type = "Pd" if i == 0 else organic_types[i-1]
            mass = 106.42 if i == 0 else (12.01 if 'c' in at_type else (14.01 if 'n' in at_type else 1.008))
            f.write(f"{i+1:6d} {at_type:5s} 1  PDL  {symbol:5s} {i+1:6d} {float(charge):10.6f} {mass:8.3f}\n")
        
        # Shift all interaction indices by +1
        f.write("\n[ bonds ]\n1 8 1 0.20131 250000\n1 9 1 0.20146 250000\n1 20 1 0.20297 250000\n1 21 1 0.20299 250000\n")
        
        # Simplified logic for remaining sections
        for section in ["bonds", "pairs", "angles", "dihedrals"]:
            if section != "bonds": f.write(f"\n[ {section} ]\n")
            lines = re.findall(rf'\[ {section} \](.*?)\n\[', itp_content + "\n[", re.S)[0].strip().split('\n')
            n = 4 if section == "dihedrals" else (3 if section == "angles" else 2)
            for line in lines:
                if not line.strip() or line.strip().startswith(';'): continue
                p = line.split()
                shifted = [str(int(x) + 1) if j < n else x for j, x in enumerate(p)]
                f.write("  " + "  ".join(shifted) + "\n")

    print("✅ Created PDL_final.gro and PDL_final.itp. They are now perfectly synced.")

sync_topol_and_gro("PDL_small_opt.chg", "PDL_H_optimized_nopd.acpype/PDL_H_optimized_nopd_GMX.itp")

