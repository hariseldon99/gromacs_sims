import sys

def log_to_pdb(log_file, pdb_output):
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: {log_file} not found.")
        return

    # Look for the LAST coordinate table (Standard or Input orientation)
    start_idx = -1
    for i in range(len(lines) - 1, -1, -1):
        if "Standard orientation:" in lines[i] or "Input orientation:" in lines[i]:
            start_idx = i + 5 
            break
            
    if start_idx == -1:
        print("Still couldn't find a coordinate table. Please check the log for 'Input orientation'.")
        return

    periodic_table = {46: "PD", 6: "C", 7: "N", 8: "O", 1: "H"}
    
    atoms = []
    for line in lines[start_idx:]:
        parts = line.split()
        # End of table is marked by a long dashed line
        if "---------------------------------------" in line or len(parts) < 6:
            break
        try:
            at_num = int(parts[1])
            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
            atoms.append((at_num, x, y, z))
        except (ValueError, IndexError):
            continue

    with open(pdb_output, 'w') as f:
        f.write(f"REMARK   Optimized coordinates extracted from {log_file}\n")
        for i, (at_num, x, y, z) in enumerate(atoms):
            symbol = periodic_table.get(at_num, "X")
            # PDB Columns: ATOM, Index, Name, Res, Chain, ResID, X, Y, Z, Occ, B, Element
            f.write(f"ATOM  {i+1:5d} {symbol:<4s} LIG A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {symbol:>2s}\n")
        f.write("END\n")
    
    print(f"✅ Success! {len(atoms)} atoms written to {pdb_output}")

# Update the path to your log file here
log_to_pdb("logs/PDL_small_opt2.log", "PDL_H_optimized.pdb")
