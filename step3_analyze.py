import os
import sys
import pandas as pd
import pyblock
import re
import numpy as np
import argparse
from argparse import ArgumentParser
from utils.blocking import BlockAnalysis
import matplotlib.pyplot as plt

import gromacs


def ParserOptions():
    parser = ArgumentParser()

    """Parse command line arguments and start main script for the workflow."""
    parser.add_argument("--repeat", dest="repeat", default=1, type=int, help="Number of simulations") 
    args = parser.parse_args()
    return args

    
def mda_edr(edrfile, term):
	import MDAnalysis as mda
	from MDAnalysisTests.datafiles import AUX_EDR
	aux = mda.auxiliary.EDR.EDRReader(edrfile)
	return aux.get_data(term)

def calc(term, repeat):
	# Calculate the means and SEMs for 'bound'
	data_bound = []
	for j in range(1, repeat + 1):
		df = pd.DataFrame.from_dict(mda_edr('bound_{}/md.edr'.format(j), term))
		data_bound.append(df)
	data_bound = pd.concat(data_bound)
	block_bound = BlockAnalysis(data_bound.iloc[:, 1], multi=repeat)
	block_bound.SEM()
	mean_bound = block_bound.av
	sem_bound = block_bound.sem
	
	# Calculate the means and SEMs for 'unbound'
	data_unbound = []
	for j in range(1, repeat + 1):
		df = pd.DataFrame.from_dict(mda_edr('unbound_{}/md.edr'.format(j), term))
		data_unbound.append(df)
	data_unbound = pd.concat(data_unbound)
	block_unbound = BlockAnalysis(data_unbound.iloc[:, 1], multi=repeat)
	block_unbound.SEM()
	mean_unbound = block_unbound.av
	sem_unbound = block_unbound.sem
	
	# Calculate the difference between means and the standard error
	diff_means = mean_bound - mean_unbound
	sem_diff = ((sem_bound ** 2) + (sem_unbound ** 2)) ** 0.5
	return [diff_means,sem_diff];

if __name__ == '__main__':

	args = ParserOptions()
	
	# Calculate the means and SEMs for 'bound'
	data_bound = []
	for j in range(1, args.repeat + 1):
		df = pd.DataFrame.from_dict(mda_edr('bound_{}/md.edr'.format(j), 'Potential'))
		data_bound.append(df)
	data_bound = pd.concat(data_bound)
	block_bound = BlockAnalysis(data_bound.iloc[:, 1], multi=args.repeat)
	block_bound.SEM()
	mean_bound = block_bound.av
	sem_bound = block_bound.sem	
	with open('blocking_bound.txt', 'w') as f:
		f.write(str(pd.DataFrame(block_bound.stat, columns=['Block size', 'SEM', 'err(SEM)'])))
	plt.errorbar(block_bound.stat[...,0], block_bound.stat[...,1], block_bound.stat[...,2], fmt='', color='k', ecolor='0.5')
	plt.scatter(block_bound.bs, block_bound.sem,zorder=10,c='tab:red')
	plt.xlabel('Block size')
	plt.ylabel('SEM')
	plt.savefig("blocking_bound.svg")
	plt.close()

	reblock_data = pyblock.blocking.reblock(data_bound)
	opt = pyblock.blocking.find_optimal_block(len(data_bound), reblock_data)
	data = pd.Series(data_bound)
	(data_length, reblock_data, covariance) = pyblock.pd_utils.reblock(data_bound)
	data_length
	reblock_data
	covariance
	#pyblock.plot.plot_reblocking(reblock_data)
	print(pyblock.pd_utils.reblock_summary(reblock_data))

	# Calculate the means and SEMs for 'unbound'
	data_unbound = []
	for j in range(1, args.repeat + 1):
		df = pd.DataFrame.from_dict(mda_edr('unbound_{}/md.edr'.format(j), 'Potential'))
		data_unbound.append(df)
	data_unbound = pd.concat(data_unbound)
	block_unbound = BlockAnalysis(data_unbound.iloc[:, 1], multi=args.repeat)
	block_unbound.SEM()
	mean_unbound = block_unbound.av
	sem_unbound = block_unbound.sem	
	with open('blocking_unbound.txt', 'w') as f:
		f.write(str(pd.DataFrame(block_unbound.stat, columns=['Block size', 'SEM', 'err(SEM)'])))
	plt.errorbar(block_unbound.stat[...,0], block_unbound.stat[...,1], block_unbound.stat[...,2], fmt='', color='k', ecolor='0.5')
	plt.scatter(block_unbound.bs, block_unbound.sem,zorder=10,c='tab:red')
	plt.xlabel('Block size')
	plt.ylabel('SEM')
	plt.savefig("blocking_unbound.svg")
	plt.close()

	cumulative_sum_bound = [data_bound.iloc[:, 1][:int(len(data_bound) * pct/100)].sum() for pct in range(1, 101)]
	cumulative_mean_bound = [cumulative_sum_bound[i] / len(data_bound[:int(len(data_bound) * pct/100)]) for i, pct in enumerate(range(1, 101))]
	cumulative_sum_unbound = [data_unbound.iloc[:, 1][:int(len(data_unbound) * pct/100)].sum() for pct in range(1, 101)]
	cumulative_mean_unbound = [cumulative_sum_unbound[i] / len(data_unbound[:int(len(data_unbound) * pct/100)]) for i, pct in enumerate(range(1, 101))]
	plt.plot(range(1, 101), np.array(cumulative_mean_bound) - np.array(cumulative_mean_unbound))
	plt.xlabel('Percentage of data')
	plt.ylabel('ΔH (kJ/mol)')
	plt.savefig("Cumulative_ΔH.svg")
	plt.close()
	
	
	# Calculate the difference between means and the standard error
	diff_means = mean_bound - mean_unbound
	sem_diff = ((sem_bound ** 2) + (sem_unbound ** 2)) ** 0.5
	
	# Print Overall ΔH
	print('\nOverall ΔH: {:.2f} ± {:.2f} kJ/mol'.format(diff_means,sem_diff))
	
	print('\nPhysical Determinants of Binding Enthalpy')
	#Bond Angle Proper-Dih. Improper-Dih.
	valence_mean = 0
	valence_sem = 0
	for i in ['Bond','Angle','Proper Dih.','Improper Dih.']:
		dat = calc(i,args.repeat)
		valence_mean += dat[0]
		valence_sem += dat[1]**2
	print('Valance: \t{:.2f} ± {:.2f} kJ/mol'.format(valence_mean,(valence_sem)** 0.5))
	
	#LJ-14 LJ-(SR) Disper.-corr.
	LJ_mean = 0
	LJ_sem = 0
	for i in ['LJ-14','LJ (SR)','Disper. corr.']:
		dat = calc(i,args.repeat)
		LJ_mean += dat[0]
		LJ_sem += dat[1]**2
	print('Lennard-Jones: \t{:.2f} ± {:.2f} kJ/mol'.format(LJ_mean,(LJ_sem)** 0.5))
	
	#Coulomb-14 Coulomb-(SR) Coul.-recip.
	Coulomb_mean = 0
	Coulomb_sem = 0
	for i in ['Coulomb-14','Coulomb (SR)','Coul. recip.']:
		dat = calc(i,args.repeat)
		Coulomb_mean += dat[0]
		Coulomb_sem += dat[1]**2
	print('Electrostatics: \t{:.2f} ± {:.2f} kJ/mol'.format(Coulomb_mean,(Coulomb_sem)** 0.5))
		


