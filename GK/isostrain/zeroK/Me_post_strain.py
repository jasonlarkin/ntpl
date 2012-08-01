## ## ## Kevin Parrish - 6/25/2012	## ## ##
## ## ## Post Processing for GK			## ## ##

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sp
from scipy import integrate
from scipy import stats
from scipy.optimize import leastsq
import numpy as np
import os

## ## ## Standard Constants
con_kb = 1.3806E-23
con_hbar = 1.054E-34
con_s2ps = 1E-12	#seconds to picoseconds

## ## ## LJ Parameters
LJ_eps = 1.67E-21
LJ_sigma = 3.4E-10
LJ_mass = 6.6326E-26
LJ_tau = sp.sqrt((LJ_mass*(LJ_sigma**2))/LJ_eps)

## ## ## Mechanics Simulation Parameters
Me_scaleJ = (LJ_eps)/((LJ_sigma**2)*LJ_tau)	##includes factor of 1/V
Me_s = 1
Me_p = 1
Me_d = GK_p * GK_s
Me_dt_LJ = 0.002
Me_dt = GK_dt_LJ * LJ_tau
Me_sample_rate = 5
Me_total_steps = 10

Me_strains = [-.12, -.10, -.08, -.06, -.04, -.02, -.01, 0, .005, .01, .015, .02, .04, .06, .08, .10, .12]
##### END Me params

Me_num_strains = len(GK_strains)

Me_stress = np.zeros( (GK_num_strains), dtype=float)
Me_volume = np.zeros( (GK_num_strains), dtype=float )

## ## ## Average the volume
	str_read = 'out.lammps.vol.'+ str(istrain)
	dummy_array = np.loadtxt(str_read, comments='#')
	GK_volume[itemp,istrain] = GK_volume[itemp,istrain]+ dummy_array[dummy_array[:,0].size - 1,1]**3	## NEED to access last element and * by 3
##### END VOLUME LOOP

## ## ## Loop over Strains
for istrain in range(GK_num_strains):
	str_read = 'out.lammps.strain.'+ str(istrain)
	dummy_array = np.loadtxt(str_read, comments='#')
	Me_stress[istrain] = (-1) * dummy_array[dummy_array[:, 0] - 1,2]
##### END STRAIN LOOP

## ## ## Print lattice constants
tlattice = np.zeros( (Me_volume[:].size), dtype=float)
for istrain in range(GK_num_strains):
	tlattice[istrain] = np.power(GK_volume[istrain],(1/3))
	#print 'Lattice vector for strain = '+ str(GK_strains[istrain])+	' and temp = '+ str(GK_temps[itemp])+ 'K is ['+ str(tlattice[itemp,istrain])+ ', '+ str(tlattice[itemp,istrain])+ ', '+ str(tlattice[itemp,istrain])+ ']'

## ## ## Save lattice constants to file
np.save('../post/post.zero.lattice.npy', tlattice)

## ## ## Convert to real units
Me_volume[:] = Me_volume[:] * (LJ_sigma**3)
Me_stress[:] = Me_stress[:] * (LJ_eps / (LJ_sigma**3))

## ## ## Save stress to file
np.save('../post/post.zero.stress.npy', Me_stress)

## ## ## Save volume to file
np.save('../post/post.zero.volume.npy', Me_volume)

## ## ## Save kappa to file
np.save('post/post.kappa.npy', kappa)

## ## ## Find Young's Modulus
def hook(m, x):
	return m*x

def resid(m, x, y):
	return y - hook(m, x)

young = np.zeros( (GK_num_temps) ) #array to hold young's moduli

loweps = 7 #index of low strain inclusive
higheps = 10 #index of high strain non-inclusive
m0 = 100000 #initial guess
x = np.zeros( (higheps - loweps), dtype=float)
x = GK_strains[loweps:higheps] #0 to 2% inclusive

for itemp in range(GK_num_temps):
	y = Me_stress[itemp, loweps:higheps]
	young[itemp],cov,infodict,mesg,ier = leastsq(resid, m0, args=(x, y), full_output=True)
	print 'At T = '+ str(GK_temps[itemp])+ ' Young\'s modulus is '+ str(young[itemp])

## ## ## Save young's modulus to file
np.save('post/post.young.npy', young)
