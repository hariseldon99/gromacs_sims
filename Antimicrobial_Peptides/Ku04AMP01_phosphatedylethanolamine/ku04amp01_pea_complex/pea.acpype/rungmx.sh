#!/bin/bash
echo 0 | gmx editconf -f pea_GMX.gro -bt octahedron -d 10 -c -princ
gmx grompp -f em.mdp -c out.gro -p pea_GMX.top -o em.tpr -v
gmx mdrun -ntmpi 1 -v -deffnm em
gmx grompp -f md.mdp -c em.gro -p pea_GMX.top -o md.tpr -r em.gro
gmx mdrun -nb gpu -gpu_id 0 -pme gpu -pmefft gpu -pin on -v -nt 16 -deffnm md
