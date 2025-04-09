
echo 0 | gmx editconf -f lfx_GMX.gro -bt octahedron -d 10 -c -princ
gmx grompp -f em.mdp -c out.gro -p lfx_GMX.top -o em.tpr -v
gmx mdrun -ntmpi 1 -v -deffnm em
gmx grompp -f md.mdp -c em.gro -p lfx_GMX.top -o md.tpr -r em.gro
gmx mdrun -nt 12 -v -deffnm md
