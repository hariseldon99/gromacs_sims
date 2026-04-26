# 1. Box and Centering
echo 0 | gmx editconf -f PDL_GMX.gro -bt octahedron -d 1.0 -c -princ -o boxed.gro

# 2. Energy Minimization (EM)
gmx grompp -f em.mdp -c boxed.gro -p PDL_GMX.top -o em.tpr -v
# Run EM using all available threads (-nt 0 uses all detected cores)
gmx mdrun -nt 6 -v -deffnm em

# 3. Production MD (GPU PME Offload)
gmx grompp -f md.mdp -c em.gro -p PDL_GMX.top -o md.tpr -r em.gro
# -nt 0: Use all CPU threads
# -nb gpu: Non-bonded interactions on GPU
# -pme gpu: PME calculations on GPU
# -pmefft gpu: PME FFT on GPU (if supported by your GMX version)
gmx mdrun -nt 6 -nb gpu -pme gpu -pmefft gpu -v -deffnm md
