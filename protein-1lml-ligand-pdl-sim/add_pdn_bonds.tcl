package require topotools

# wait until molecule is loaded
set mol [molinfo top]

# add Pd–N bonds
topo addbond 6898 6905
topo addbond 6898 6906
topo addbond 6898 6917
topo addbond 6898 6918

# optional: rebuild topology graph for display purposes
topo guessangles
topo guessdihedrals
