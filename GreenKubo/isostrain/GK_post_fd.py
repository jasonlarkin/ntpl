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

## ## ## Green Kubo Simulation Parameters
GK_scaleJ = (LJ_eps)/((LJ_sigma**2)*LJ_tau)	##includes factor of 1/V
GK_s = 5
GK_p = 25000
GK_d = GK_p * GK_s
GK_dt_LJ = 0.002
GK_dt = GK_dt_LJ * LJ_tau
GK_sample_rate = 5
GK_total_steps = 1000000

GK_seeds = [22222, 33333, 44444, 55555, 66666, 77777, 88888, 99999, 10101010, 11111111]
GK_temps = [20, 50, 80, 100]
GK_strains = [-.12, -.10, -.08, -.06, -.04, -.02, -.01, 0, .005, .01, .015, .02, .04, .06, .08, .10, .12]

#GK_seeds = [55555, 66666]
#GK_temps = [20, 50]
#GK_strains = [-.32, 0]
##### END GK params

## ## ## Mechanics Simulation Parameters
Me_p = 125000
Me_total_steps = GK_total_steps
Me_d = Me_total_steps / Me_p
##### END Me params

GK_num_seeds = len(GK_seeds)
GK_num_temps = len(GK_temps)
GK_num_strains = len(GK_strains)

GK_JJ = np.zeros( (GK_p, 2, GK_num_temps, GK_num_strains), dtype=float )
GK_JJ[:,0,0,0] = np.arange(0, (GK_JJ[:,0,0,0].size), 1) * GK_dt * GK_sample_rate
GK_volume = np.zeros( (GK_num_temps, GK_num_strains) )


## ## ## Loop over Strains
for istrain in range(GK_num_strains):
	## ## ## Loop over Temperatures
	for itemp in range(GK_num_temps):
		## ## ## Loop over Seeds
		for iseed in range(GK_num_seeds):
			## ## ## Use Grep to extract 1000 steps into a readable file
			str_cmd = 'grep -A '+ str(GK_p)+ ' \"'+ str(GK_total_steps)+ ' '+ str(GK_p)+ '\" out.lammps.corr.'+ str(iseed+1)+ '.'+ str(itemp+1)+ '.'+ str(istrain+1)+ ' > out.lammps.corr.'+ str(iseed+1)+ '.'+ str(itemp+1)+ '.'+ str(istrain+1)+ '.grep'
			os.system(str_cmd)
			## ## ## Read the file and average
			str_read = 'out.lammps.corr.'+ str(iseed+1)+ '.'+ str(itemp+1)+ '.'+ str(istrain+1)+ '.grep'
			dummy_array = np.loadtxt(str_read, skiprows = 1)
			GK_JJ[:,1,itemp,istrain] = GK_JJ[:,1,itemp,istrain]+ ((dummy_array[:,3]+ dummy_array[:,4]+ dummy_array[:,5])/3) 
			## ## ## Average the volume
			str_read = 'out.lammps.vol.'+ str(iseed+1)+ '.'+ str(itemp+1)+ '.'+ str(istrain+1)
			dummy_array = np.loadtxt(str_read, comments='#')
			GK_volume[itemp,istrain] = GK_volume[itemp,istrain]+ dummy_array[dummy_array[:,0].size - 1,1]**3	## NEED to access last element and * by 3
			##### END SEED LOOP
		GK_JJ[:,1,itemp,istrain] = GK_JJ[:,1,itemp,istrain] / GK_num_seeds	##Normalize by # of seeds
		GK_volume[itemp,istrain] = GK_volume[itemp,istrain] / GK_num_seeds	##Normalize by # of seeds
		GK_JJ[:,1,itemp,istrain] = GK_JJ[:,1,itemp,istrain] / (GK_volume[itemp,istrain]**2)	## USE ONLY IF NOT DIVIDED BY VOL IN LAMMPS
		##### END TEMP LOOP
	##### END STRAIN LOOP

Me_stress = np.zeros( (GK_num_temps, GK_num_strains), dtype=float)
## ## ## Loop over Strains
for istrain in range(GK_num_strains):
	## ## ## Loop over Temperatures
	for itemp in range(GK_num_temps):
		## ## ## Loop over Seeds
		for iseed in range(GK_num_seeds):
			str_read = 'out.lammps.strain.'+ str(iseed+1)+ '.'+ str(itemp+1)+ '.'+ str(istrain+1)
			dummy_array = np.loadtxt(str_read, comments='#')
			Me_stress[itemp, istrain] = Me_stress[itemp, istrain] + dummy_array[Me_d - 1,2]
		Me_stress[itemp, istrain] = (-1) * Me_stress[itemp, istrain] / GK_num_seeds ##Normalize by # of seeds, convert to positive
		##### END TEMP LOOP
	##### END STRAIN LOOP

## ## ## Print lattice constants
tlattice = np.zeros( (GK_volume[:,0].size, GK_volume[0,:].size), dtype=float)
for istrain in range(GK_num_strains):
	for itemp in range(GK_num_temps):
		tlattice[itemp, istrain] = np.power(GK_volume[itemp,istrain],(1/3))
		#print 'Lattice vector for strain = '+ str(GK_strains[istrain])+	' and temp = '+ str(GK_temps[itemp])+ 'K is ['+ str(tlattice[itemp,istrain])+ ', '+ str(tlattice[itemp,istrain])+ ', '+ str(tlattice[itemp,istrain])+ ']'

## ## ## Save lattice constants to file
np.save('post/post.lattice.npy', tlattice)

## ## ## Convert to real units
GK_volume[:,:] = GK_volume[:,:] * (LJ_sigma**3)
GK_JJ[:,1,:,:] = GK_JJ[:,1,:,:] * (GK_scaleJ**2)
Me_stress[:,:] = Me_stress[:,:] * (LJ_eps / (LJ_sigma**3))

## ## ## Save stress to file
np.save('post/post.stress.npy', Me_stress)

## ## ## Save volume to file
np.save('post/post.volume.npy', GK_volume)

## ## ## First Dip Method
GK_intJJ = np.zeros( (GK_JJ[:,0,0,0].size, GK_num_temps, GK_num_strains), dtype=float)	## Create array for integration

for istrain in range(GK_num_strains):	## Integrate the function
	for itemp in range(GK_num_temps):
		GK_intJJ[1:,itemp,istrain] = integrate.cumtrapz(GK_JJ[:,1,itemp,istrain], x = GK_JJ[:,0,0,0]) * (GK_volume[itemp,istrain] / (con_kb * (GK_temps[itemp]**2)))

## ## ## Find kappa
kappa = np.zeros( (GK_num_temps, GK_num_strains) )
warn = np.zeros( (GK_num_temps, GK_num_strains) )

for istrain in range(GK_num_strains):
	for itemp in range(GK_num_temps):
		z = np.where(GK_JJ[:,1,itemp,istrain] < 0)[0]	##Find first instance of a negative number
		if z.size is 0:
			uplim = GK_p
			warn[itemp, istrain] = 1 
			print 'Warning: T = '+str(GK_temps[itemp])+' K, $\epsilon$ = '+str(GK_strains[istrain])+ ' exceeds data set'
		elif z[0]+5000 > GK_p:
			uplim = GK_p
			warn[itemp, istrain] = 1
			print 'Warning: T = '+str(GK_temps[itemp])+' K, $\epsilon$ = '+str(GK_strains[istrain])+ ' exceeds data set'
		else:
			uplim = z[0]+5000
		kappa[itemp,istrain] = sp.mean(GK_intJJ[0:uplim,itemp,istrain])
		#kappa[itemp,0] = sp.mean(GK_intJJ[:,itemp])	##Use whole integral
		#kappa[itemp,istrain,1] = stats.sem(GK_intJJ[0:uplim,itemp,istrain])	## Standard error

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

## ## ## Plot HCACF const strain
for istrain in range(GK_num_strains):
	for itemp in range(GK_num_temps):
		l1 = plt.plot(GK_JJ[:,0,0,0]/con_s2ps, GK_JJ[:,1,itemp,istrain]/GK_JJ[0,1,itemp,istrain], label= ('T = '+ str(GK_temps[itemp])+ ' K'))
	ll = plt.legend(loc='upper right')
	ly = plt.ylabel('<q(t) * q(0)> / <q(0) * q(0)>')
	lx = plt.xlabel('Time, ps')
	lxm = plt.xlim(0, 5)
	lt = plt.title('HCACF over temperature at strain of '+ str(GK_strains[istrain]))
	plt.savefig(('post/pic.hcacf.strain.'+ str(istrain+1)+ '.png'))
	plt.clf()

## ## ## Plot HCACF const temp
for itemp in range(GK_num_temps):
	for istrain in range(GK_num_strains):
		l1 = plt.plot(GK_JJ[:,0,0,0]/con_s2ps, GK_JJ[:,1,itemp,istrain]/GK_JJ[0,1,itemp,istrain], label= ('$\sigma$ = '+ str(GK_strains[istrain])))
	ll = plt.legend(loc='upper right')
	ly = plt.ylabel('<q(t) * q(0)> / <q(0) * q(0)>')
	lx = plt.xlabel('Time, ps')
	lxm = plt.xlim(0, 5)
	lt = plt.title('HCACF over strain at temperature of '+ str(GK_temps[itemp])+ 'K')
	plt.savefig(('post/pic.hcacf.temp.'+ str(itemp+1)+ '.png'))
	plt.clf()

## ## ## Plot Thermal Conductivity const strain
for istrain in range(GK_num_strains):
	for itemp, t in zip(range(GK_num_temps), ['r','g','c','m']):
		l1 = plt.plot(GK_JJ[:,0,0,0]/con_s2ps, GK_intJJ[:,itemp,istrain], c=t, label= ('T = '+ str(GK_temps[itemp])+ ' K'))
		l2 = plt.axhline(y=kappa[itemp,istrain], c=t, ls=':')
	ll = plt.legend(loc='lower right')
	ly = plt.ylabel('$\kappa$ (W/m-K)')
	lx = plt.xlabel('Time, ps')
	lxm = plt.xlim(0, 200)
	lt = plt.title('Thermal Conductivity over temperature at strain of '+ str(GK_strains[istrain]))
	plt.savefig(('post/pic.cond.strain.'+ str(istrain+1)+ '.png'))
	plt.clf()

## ## ## Plot Thermal Conductivity const temp
for itemp in range(GK_num_temps):
	for istrain in range(GK_num_strains):
		l1 = plt.plot(GK_JJ[:,0,0,0]/con_s2ps, GK_intJJ[:,itemp,istrain], label= ('$\sigma$ = '+ str(GK_strains[istrain])))
	ll = plt.legend(loc='lower right')
	ly = plt.ylabel('$\kappa$ (W/m-K)')
	lx = plt.xlabel('Time, ps')
	lxm = plt.xlim(0, 200)
	lt = plt.title('Thermal Conductivity over strain at temperature of '+ str(GK_temps[itemp])+ 'K')
	plt.savefig(('post/pic.cond.temp.'+ str(itemp+1)+ '.png'))
	plt.clf()

