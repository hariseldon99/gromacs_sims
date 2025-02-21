<<<<<<< HEAD

=======
#!/bin/bash
>>>>>>> 754d05ef0c287f1187842b73cac4b171a4f1acd6
echo 0 | gmx editconf -f pgl_GMX.gro -bt octahedron -d 10 -c -princ
gmx grompp -f em.mdp -c out.gro -p pgl_GMX.top -o em.tpr -v
gmx mdrun -ntmpi 1 -v -deffnm em
gmx grompp -f md.mdp -c em.gro -p pgl_GMX.top -o md.tpr -r em.gro
<<<<<<< HEAD
gmx mdrun -ntmpi 1 -v -deffnm md
=======
gmx mdrun -nb gpu -gpu_id 0 -pme gpu -pmefft gpu -pin on -v -nt 16 -deffnm md
>>>>>>> 754d05ef0c287f1187842b73cac4b171a4f1acd6
