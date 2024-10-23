import os
import sys
import argparse
import re

WD = os.getcwd()
MDPFILES = sys.path[0] + '/mdps'


def find_gmx():
	RC = os.system('gmx -h >/dev/null 2>&1')
	if RC == 0:
		return 'gmx'
	RC = os.system('gmx_mpi -h >/dev/null 2>&1')
	if RC == 0:
		return 'gmx_mpi'
	raise SystemExit("\nERROR! Not found gmx or gmx_mpi.\n")


def pull_code(mdpfile, state):
	with open(MDPFILES + '/' + mdpfile, 'r') as file:
		lines = file.readlines()

	with open(mdpfile, 'w') as f:
		for line in lines:
			if state == 'bound' and line.startswith("pull"):
				line = line.replace('pull', ';pull')
			if state == 'unbound' and line.startswith(";pull"):
				line = line.replace(';pull', 'pull')
			f.write(line)

def COM_distance(state): #Calculating min distance between center_of_masses.
	import MDAnalysis as mda
	from MDAnalysis.analysis import distances
	
	u =  mda.Universe(state + '_1/md.tpr',state + '_1/md.gro')
	
	host = u.select_atoms("resname HHH").center_of_mass()
	guest = u.select_atoms("resname GGG").center_of_mass()

	distances = abs(host[0]-guest[0])
	return "%.2f" % (float(distances) / 10)

def protocol(state, ID, threads):
	RC = os.system('cp -r {0}/{1} {0}/{1}_{2};'.format(WD, state, ID))
	if RC != 0:
		raise SystemExit("\nERROR! Can't copy the input, check your {} directory.\n".format(state))
	os.chdir('{}/{}_{}'.format(WD, state, ID))


	GMX = find_gmx()
	
	if ID == 0:
		#Minimization
		RC = os.system('''{0} grompp -f {1}/em.mdp -c solv.gro -p topol.top -r solv.gro -o em.tpr -maxwarn 2 > grompp-em.log 2>&1
				          {0} mdrun -nt {2} -deffnm em > mdrun-em.log 2>&1'''.format(GMX, MDPFILES, threads))
		if RC != 0:
			raise SystemExit('\nERROR! see the log files\n')
			
		# Equilibration NVT
		RC = os.system('''{0} grompp -f {1}/nvt.mdp -c em.gro -t em.trr -p topol.top -o nvt.tpr -r em.gro -n index.ndx -maxwarn 2 > grompp-nvt.log 2>&1
				          {0} mdrun -nt {2} -deffnm nvt > mdrun-nvt.log 2>&1'''.format(GMX, MDPFILES, threads))
		if RC != 0:
			raise SystemExit('\nERROR! see the log files\n')
		
		# Equilibration NPT
		RC = os.system('''{0} grompp -f {1}/npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -r nvt.gro -n index.ndx -maxwarn 2 > grompp-npt.log 2>&1
				          {0} mdrun -nt {2} -deffnm npt > mdrun-npt.log 2>&1'''.format(GMX, MDPFILES, threads))
		if RC != 0:
			raise SystemExit('\nERROR! see the log files\n')
	
	if ID != 0:
		pull_code('npt-norest.mdp',state)
		
		# Equilibration NPT without restraints
		RC = os.system('''{0} grompp -f npt-norest.mdp -c {2}/{3}_0/npt.gro -t {2}/{3}_0/npt.cpt -p topol.top -o npt-norest.tpr -n index.ndx -maxwarn 2 > grompp-npt-norest.log 2>&1
				          {0} mdrun -nt {4} -deffnm npt-norest > mdrun-npt-norest.log 2>&1'''.format(GMX, MDPFILES, WD, state, threads))
		if RC != 0:
			raise SystemExit('\nERROR! Check the log files\n')
		
		pull_code('md.mdp',state)
		
		# Production
		RC = os.system('''{0} grompp -f md.mdp -c npt-norest.gro -t npt-norest.cpt -p topol.top -o md.tpr -n index.ndx -maxwarn 2 > grompp-md.log 2>&1
				          {0} mdrun -nt {2} -deffnm md > mdrun-md.log 2>&1'''.format(GMX, MDPFILES, threads))
		if RC != 0:
			raise SystemExit('\nERROR! see the log files\n')


def ParserOptions():
    parser = argparse.ArgumentParser()

    """Parse command line arguments and start main script for the workflow."""
    parser.add_argument("--system", dest="system", help="System type; host-guest (HG), protein-ligand (PL), protein-peptide (PP)", required=True)   
    parser.add_argument("--index", dest="index", default=1, type=int, help="Simulation index")   
    parser.add_argument("--length", dest="length", default=100, type=int, help="The length of each simulation repeat") 
    parser.add_argument("--state", dest="state", help="bound or unbound", required=True)   
    parser.add_argument("--nt", dest="nt", default=1, type=int, help="Total number of threads to start the simulations")   

    
    args = parser.parse_args()
    return args


if __name__ == '__main__':

	args = ParserOptions()
	
	print(COM_distance('unbound'))
	exit()
	
	protocol(args.state, args.index, args.nt)




