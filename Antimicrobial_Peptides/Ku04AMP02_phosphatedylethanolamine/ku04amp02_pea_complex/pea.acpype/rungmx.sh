
#!/bin/bash
echo 0 | gmx editconf -f pea_GMX.gro -bt octahedron -d 1 -c -princ
gmx grompp -f em.mdp -c out.gro -p pea_GMX.top -o em.tpr -v
gmx mdrun -nt 12 -v -deffnm em
gmx grompp -f md.mdp -c em.gro -p pea_GMX.top -o md.tpr -r em.gro
gmx mdrun -nt 12 -v -deffnm md
