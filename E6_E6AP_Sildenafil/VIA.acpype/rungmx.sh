
echo 0 | gmx editconf -f VIA_GMX.gro -bt octahedron -d 1 -c -princ
gmx grompp -f em.mdp -c out.gro -p VIA_GMX.top -o em.tpr -v
gmx mdrun -nt 12 -pin on -nb gpu -gpu_id 0 -v -deffnm em
gmx grompp -f md.mdp -c em.gro -p VIA_GMX.top -o md.tpr -r em.gro
gmx mdrun -nt 12 -pin on -nb gpu -gpu_id 0 -v -deffnm md
