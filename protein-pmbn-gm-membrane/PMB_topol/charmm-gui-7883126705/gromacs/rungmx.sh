echo 0 | gmx editconf -f ../ligandrm.pdb -bt octahedron -d 1 -c -princ
gmx grompp -f em.mdp -c out.gro -p topol.top -o em.tpr -v -maxwarn 9
gmx mdrun -nt 12 -nb gpu -gpu_id 0 -v -deffnm em 
gmx grompp -f md.mdp -c em.gro -p topol.top -o md.tpr -r em.gro -maxwarn 9
gmx mdrun -nt 12 -nb gpu -gpu_id 0 -v -deffnm md 

# Post-processing gromacs, center trajectory and remove artefacts
# Create the ndx file for the trjconv command
echo "q" | gmx make_ndx -f md.tpr -o index.ndx
# Center on PMB (group 2), output System (group 0)
echo "2 0" | gmx trjconv -s md.tpr -f md.trr -o md_noPBC.trr -pbc mol -center -ur compact -n index.ndx
