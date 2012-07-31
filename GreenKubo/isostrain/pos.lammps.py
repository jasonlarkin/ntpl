## ## ## Kevin Parrish - 7/31/12 - lammps pos script ## ## ##
import os
import sys
sys.path.append('/home/kevin/ntpl') # Needed to recognize ntplpy module
import ntplpy.lattice as lt

## Set lattice constants
lat20 = 1.5632
lat50 = 1.5885
lat80 = 1.6256
lat100 = lat80 # No lattice constant for liquid

lats = [lat20, lat50, lat80, lat100]

## Loop through creating position files
for index, lat in zip(range(4), lats):
	mylamp = lt.Lammps('in.pos.'+ str(index))
	mylat = lt.Block([lats[index], lats[index], lats[index]], 'fcc', [4, 4, 4], atom_mass = [1])
	mylamp.buildLammps([mylat])



