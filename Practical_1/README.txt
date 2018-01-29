Format of command line input:
huckel.py -l <linear length> -c <cyclic length> -f <filename>

output is given in command line and as a text file with the format "out_<inputfilename>.txt" and "out_<linear/cyclic>_polyene_c<#atoms>.txt"
-------------------------------------
The format of input files is of a unique alpha-numeric identifier for each atom in the molecule, seperated from adjacent atoms by a ",".  Non-adjacent atoms are seperated by a ";".  Each neighboring atom must be listed preceding or following the atom it neighbors at least once.

For example:

linear triatomic:

A,b,3

triangular geometry:

0,1,2,0

OR

0,1,2;0,2

rhombus geometry:

0,1,2,3,0;1,3

OR

A,B,2,def;A,def,B

-------------------------------------------------

Some example text files have been included for the following geometries:
cube
dodecahedron
icosahedron
tetrahedron
fullerene_c60 (buckminster fullerene)

