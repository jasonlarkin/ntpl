## ## ## Kevin Parrish - 5/16/12 - in.GK.LJ editor ## ## ##
import os

strain = [-0.12, -0.10, -0.08, -0.06, -0.04, -0.02, -0.01, 0, 0.05, 0.01, 0.015, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12]

numStrains = len(strain)

for istrain in range(numStrains):
	## ## ## Change strain 
	orig = 'STRAIN_PARAM'
	new  = str(strain[istrain])
	cmd1 = '-e \'s/' + orig + '/' + new + '/g\' '

	## ## ## Change output file
	orig = 'OUT_VOLUME'
	new  = 'out.lammps.vol.' + str(istrain)
	cmd2 = '-e \'s/' + orig + '/' + new + '/g\' '

	## ## ## Change output file
	orig = 'OUT_STRAIN'
	new  = 'out.lammps.strain.' + str(istrain)
	cmd3 = '-e \'s/' + orig + '/' + new + '/g\' '

	## ## ## Change position file
	orig = 'OUT_POS'
	new  = 'out.lammps.pos.' + str(istrain) + '.xyz'
	cmd4 = '-e \'s/' + orig + '/' + new + '/g\' '

	## ## ## Change main log file
	orig = 'MAIN_LOG_FILE'
	new  = 'out.lammps.main.' + str(istrain)
	cmd5 = '-e \'s/' + orig + '/' + new + '/g\' '
	
	## ## ## Execute Command
	os.system('sed '+ cmd1+ cmd2+ cmd3+ cmd4+ cmd5+ 'in.lammps.zero.temp > ./in.lammps.zero.'+ str(istrain))
		#sed -e 's/\<POS_IN_STRUCTURE\>/in_structure_1.pos/g' -e 's/\<OUT_CORR\>/out.GK.LJ.corr.1/g' in.GK.LJ > ./in.GK.LJ.1
## ## ## END FILE
