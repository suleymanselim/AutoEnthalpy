import os, shutil, glob
import sys
import re
import numpy as np
import argparse
from pathlib import Path
from rdkit.Chem.rdmolfiles import *

def find_gmx():
    RC = os.system('gmx -h >/dev/null 2>&1')
    if RC == 0:
        return 'gmx'
    RC = os.system('gmx_mpi -h >/dev/null 2>&1')
    if RC == 0:
        return 'gmx_mpi'
    logging.error('Not found gmx or gmx_mpi.')
    exit(1)

    
GMX = find_gmx()

water_models = {
  'tip3p-fb': 'tip3p', 'bind3p': 'tip3p', 'spce': 'tip3p',
  'tip3p': 'tip3p', 'opc3': 'tip3p', 'spc': 'tip3p', 'opc4': 'tip4p',
  'tip4p-fb': 'tip4p', 'tip4p': 'tip4p', 'tip4pew': 'tip4p',
  'tip5pe': 'tip5p', 'tip5p-2018': 'tip5p', 'tip5p': 'tip5p', 'tip7p': 'tip7p'
}

   
def ParserOptions():
    parser = argparse.ArgumentParser()

    """Parse command line arguments and start main script for the workflow."""
    parser.add_argument("--system", dest="system", nargs=3, help="System type; host-guest (HG), protein-ligand (PL), protein-peptide (PP)", required=True)   
    parser.add_argument("--space", dest="space", type=float, default=1, help="Distance between molecules")
    parser.add_argument("--water", dest="water", default='tip3p', help="Set water model") 
    parser.add_argument("--d", dest="d", type=float, default=1.2, help="Distance between the solute and the box")
    parser.add_argument("--FF", dest="FF", help="Generalized Force Fields GAFF, CGenFF, OpenFF", choices=['gaff', 'gaff2','openff', 'cgenff'])
    parser.add_argument("--proteinFF", dest="proteinFF", default='amber99sb-ildn', help="Force fields for protein")
    parser.add_argument("--userFF", dest="userFF", action='store_true', help="Enable user defined Force Fields for protein")
    args = parser.parse_args()
    
    return args



def convert_to_mol2(input_file, resname, output_file):
	"""Converting to mol2 file format, cgenff parametrization requires mol2 file format"""
	import openbabel
	openbabel.obErrorLog.StopLogging()
	
	# Determine the file format using file extension
	input_format = input_file.split('.')[-1].lower()
	if input_format not in ['mol2', 'pdb', 'sdf', 'mol']:
		raise SystemExit(f"Unsupported input format: {input_format}")

	mol2 = openbabel.OBConversion()
	mol2.SetInAndOutFormats(input_format, "mol2")
	mol_mol2 = openbabel.OBMol()
	mol2.ReadFile(mol_mol2, input_file) 

	for residue in openbabel.OBResidueIter(mol_mol2):
		residue.SetName(resname)
	mol_mol2.SetTitle(resname)
    
	mol2.WriteFile(mol_mol2, output_file)

def MolFromInput(mol_input): #Reading any file format with RDKit
	"""Reading molecule file from any format (mol, mol2, pdb, sdf) as RDKit molecule"""
	FILE_PARSERS = {
    "mol": MolFromMolFile,
    "mol2": MolFromMol2File,
    "pdb": MolFromPDBFile,
    "sdf": MolFromMolFile}

	if os.path.isfile(mol_input):
		content_reader = FILE_PARSERS
		mol_format = os.path.splitext(mol_input)[1][1:]
	if mol_format:
		try:
			reader = content_reader[mol_format.lower()]
		except KeyError:
			raise TypeError(
				f"Molecule format {mol_format} not supported. "
				f"Supported formats: " + ", ".join(FILE_PARSERS))
		return reader(mol_input)

	for reader in content_reader.values():
		try:
			mol = reader(mol_input)
		except RuntimeError:
			pass
		else:
			if mol is not None:
				return mol

	raise TypeError(
		f"Could not create an RDKit molecule from {mol_input}. "
        "Try passing in a `mol_format`. "
		f"Supported formats: " + ", ".join(FILE_PARSERS)) #The code was adapted from https://pypi.org/project/rdkit-utilities/


def get_gaff(molfile, molname, gaff):
    """Build a ligand topology and coordinate file from a ligand file using acpype for gaff"""
    from rdkit import Chem
   
    paras = {
        'molfile': molfile,
        'molname': molname,
        'gaff' : gaff,
        'net_charge': Chem.GetFormalCharge(MolFromInput(molfile))
            }
    
    RC = os.system('''acpype -i {molfile} -b {molname} -n {net_charge} -f -a {gaff} -o gmx >acpype.{molname}.log 2>&1 '''.format(**paras))
    if RC != 0:
        raise SystemExit('\nERROR!\nFailed to run the acpype. see the %s for details.'%os.path.abspath("acpype.{0}.log\n".format(molname)))

def get_cgenff(molecule, molname):
	"""Build a ligand topology and coordinate file from a ligand file using SILCSBIO for cgenff"""
    
	cmd = '''$SILCSBIODIR/cgenff/cgenff_to_gmx.sh mol={} > cgenff.log 2>&1;'''.format(molecule)
	RC = os.system(cmd)
	if RC != 0:
		raise SystemExit('\nERROR! see the log file for details %s'%os.path.abspath("cgenff.log\n"))
	else:
		os.system('''
		mv ./charmm36.ff/*_ffbonded.itp {0}_ffbonded.itp;
		mv {1}_gmx.top {0}.itp;
		mv posre.itp posre_{0}.itp;'''.format(molname, Path(molecule).stem))
		
		from openbabel import openbabel
		openbabel.obErrorLog.StopLogging()
		obConversion = openbabel.OBConversion()
		obConversion.SetInAndOutFormats("pdb", "gro")

		mol = openbabel.OBMol()
		obConversion.ReadFile(mol, '{}_gmx.pdb'.format(Path(molecule).stem))
		obConversion.WriteFile(mol, '''{}.gro'''.format(molname))

def get_openff(molecule, molname):
	"""Build a ligand topology and coordinate file from a ligand file using openff"""
	from acpype.topol import AbstractTopol, ACTopol, MolTopol, header
	from openff.toolkit.topology.molecule import Molecule
	from rdkit import Chem

	rdkit_mol2 = MolFromInput(molecule)
	
	# Use the OpenFF Toolkit to generate a molecule object from a SMILES pattern
	molecule = Molecule.from_rdkit(rdkit_mol2)

	# Create an OpenFF Topology object from the molecule
	from openff.toolkit.topology import Topology
	topology = Topology.from_molecules(molecule)

	# Load the latest OpenFF force field release: version 2.0.0, codename "Sage"
	from openff.toolkit.typing.engines.smirnoff import ForceField
	forcefield = ForceField('openff-2.1.0.offxml')

	# Create an Interchange object
	from openff.interchange import Interchange
	out = Interchange.from_smirnoff(force_field=forcefield, topology=topology)

	# write to GROMACS files
	out.to_gro("{}.gro".format(molname))
	out.to_top("{}.itp".format(molname))
	#AbstractTopol.writeGromacsTopolFiles(out)
	

def hg_topol(host_itp, guest_itp, ff='gaff'):
	
	if ff == 'gaff':
		defaults = """
; Include forcefield parameters
#include "amber99sb-ildn.ff/forcefield.itp"\n
"""
	if ff == 'openff':
		defaults = """
#include "/pub/scinarog/amber14sb.ff/forcefield.itp"\n
"""
	if ff == 'cgenff':
		defaults = """
; Include forcefield parameters
#include "../charmm36.ff/forcefield.itp"
;
; Include ligand specific parameters
# include "HHH_ffbonded.itp"
# include "GGG_ffbonded.itp"
"""			


	template = """{0}
{1}

#include "atomtypes_waters.itp"

#include "HHH.itp"
#ifdef POSRES
#include "posre_HHH.itp"
#endif\n
#include "GGG.itp"
#ifdef POSRES
#include "posre_GGG.itp"
#endif\n
#include "{2}.itp"
#ifdef POSRES_WATER
[ position_restraints ]
   1    1       1000       1000       1000
#endif\n
[ system ]
Host-Guest
[ molecules ]
HHH 1
GGG 1
""".format(defaults, combine_atomtypes([host_itp, guest_itp], ff), args.water)


	for output in ('bound', 'unbound'):
		try:
			os.mkdir(output)
		except FileExistsError:
			pass
		
		f = open("{0}/{1}".format(output, 'topol.top'), "w")
		f.write(template)
		f.close()
		
		mol_itp(host_itp, 'HHH', "{0}/HHH.itp".format(output), ff)
		mol_itp(guest_itp, 'GGG', "{0}/GGG.itp".format(output), ff)
		shutil.copyfile(os.path.dirname(os.path.abspath(__file__))+'/waters/itp/atomtypes_waters.itp',"{0}/atomtypes_waters.itp".format(output))
		shutil.copyfile(os.path.dirname(os.path.abspath(__file__))+'/waters/itp/{0}.itp'.format(args.water),"{0}/{1}.itp".format(output, args.water))
		if ff == 'gaff':
			shutil.copyfile('HHH.acpype/posre_HHH.itp', "{0}/posre_HHH.itp".format(output))
			shutil.copyfile('GGG.acpype/posre_GGG.itp', "{0}/posre_GGG.itp".format(output))
		elif ff == 'openff':
			for grofile in ['HHH','GGG']:
				with open(grofile+'.gro') as gro:
					file_lines = gro.readlines()
					posreFile = open("{}/posre_{}.itp".format(output,grofile), "w")
					posreFile.write("\n[ position_restraints ]\n; atom  type      fx      fy      fz\n")
					for line in file_lines[2:-2]:
						line_list = line.split()
						atom_name = str(line_list[1])
						if not atom_name.startswith("H"):
							posreFile.write("{0}     1 {1} {1} {1}\n".format(line_list[2], 1000))
					posreFile.close()
		else:
			for mol in ['HHH','GGG']:
				shutil.copyfile('posre_{}.itp'.format(mol), "{0}/posre_{1}.itp".format(output, mol))
				shutil.copyfile("{}_ffbonded.itp".format(mol), "{0}/{1}_ffbonded.itp".format(output, mol))
	
########## Combine atomtypes from itp files
def combine_atomtypes(itp_files, ff):
	'''Combine [ atomtypes ] from **itp_files** a list of itp files '''

	p = False
	atomtypes = []
	for f in itp_files:
		resname = Path(f).stem
		with open(f, 'r') as f:
			for line in f:
				if re.search('\[ atomtypes \]', line):
					p = True
					atomtypes.append(line)
					continue
				elif re.search('\[ moleculetype \]', line):
					p = False
					break
				if p:
					if line.strip() != '':
						if ff == 'openff':
							line_rn = line.replace("MOL0","")
							atomtypes.append(line_rn)
						else:
							atomtypes.append(line)

	return ''.join(map(str, list(dict.fromkeys(atomtypes))))


########## Getting itp file
def mol_itp(itp_file, resname, output_itp, ff):
	p = False
	rename = False
	with open(output_itp, 'w') as output:
		with open(itp_file, 'r') as file:
			for line in file:
				if re.search('\[ moleculetype \]', line):
					output.write(line)
					output.write('''{}              3\n'''.format(resname))
					continue
				elif line[0] == "#":
					p = False
				elif re.search('\[ system \]', line):
					p = False
					break
				elif re.search('\[ atoms \]', line):
					p = True
					rename = True
					output.write(line)
				elif re.search('\[ bonds \]', line) or re.search('\[ pairs \]', line):
					rename = False
					output.write(line)
				elif p == True and rename == False:
					output.write(line)
				elif p == True and rename == True:
					if line[0] != ";" and line[0] != "#" and line[0] != "[" and line.strip() != "":
						line_list = line.split()
						res_name = line_list[3]
						line_rn = line.replace(res_name,resname)
						output.write(line.replace(res_name,resname))


def pp_itp(itp_file, molname, output_itp):
	p = False
	with open(output_itp, 'w') as output:
		with open(itp_file, 'r') as file:
					for line in file:
						if re.search('\[ moleculetype \]', line):
							output.write(line)
							output.write(''';name            nrexcl 
{}              3\n'''.format(molname))
							continue
						elif re.search('\[ atoms \]', line):
							p = True
							output.write(line)
						elif line[0] == "#":
							p = False
						elif re.search('\[ system \]', line):
							p = False
							break
						elif p == True:
							output.write(line)


def pp_topol(mol1_itp, mol2_itp,  ff='amber99sb-ildn', userff=False):
	
	if userff:
		ff = '../'+ff
		
	template = """
; Include forcefield parameters
#include "{}.ff/forcefield.itp\n	
#include "mol1.itp"
#ifdef POSRES
#include "posre_mol1.itp"
#endif\n
#include "mol2.itp"
#ifdef POSRES
#include "posre_mol2.itp"
#endif\n
#include "{}.itp"
#ifdef POSRES_WATER
[ position_restraints ]
   1    1       1000       1000       1000
#endif\n
[ system ]
PP
[ molecules ]
mol1 1
mol2 1
""".format(ff,  args.water)
	
	for output in ('bound', 'unbound'):
		try:
			os.mkdir(output)
		except FileExistsError:
			pass
		
		f = open("{0}/{1}".format(output, 'topol.top'), "w")
		f.write(template)
		f.close()
		
		for mol in ('mol1','mol2'):
			pp_itp('{}.top'.format(mol), mol, "{0}/{1}.itp".format(output, mol))
			shutil.copyfile('posre_{}.itp'.format(mol), "{0}/posre_{1}.itp".format(output, mol))
			
		shutil.copyfile(os.path.dirname(os.path.abspath(__file__))+'/waters/itp/atomtypes_waters.itp',"{0}/atomtypes_waters.itp".format(output))
		shutil.copyfile(os.path.dirname(os.path.abspath(__file__))+'/waters/itp/{0}.itp'.format(args.water),"{0}/{1}.itp".format(output, args.water))

def pl_topol(pro_itp, lig_itp,  lig_gro, ff='gaff', proteinFF='amber99sb-ildn', userff=False):
	
	if userff:
		proteinFF = '../'+proteinFF
	
	if ff == 'gaff' or ff == 'openff':
		initial = """
; Include forcefield parameters
#include "{}.ff/forcefield.itp"\n
    """.format(proteinFF)
	elif ff == 'cgenff':
		if proteinFF == 'charmm36':
			initial = """
; Include forcefield parameters
#include "../charmm36.ff/forcefield.itp"
;
; Include ligand specific parameters
# include "LIG_ffbonded.itp"
    """
		else:
			initial = """
; Include forcefield parameters
#include "{}.ff/forcefield.itp

; Include ligand specific parameters
# include "LIG_ffbonded.itp""\n
    """.format(proteinFF)
	else:
    # Handle the case when ff is neither 'gaff' nor 'cgenff'
		raise ValueError("Invalid value for -ff")	


	template = """{0}
{1}

#include "atomtypes_waters.itp"\n
#include "mol1.itp"
#ifdef POSRES
#include "posre_mol1.itp"
#endif\n
#include "LIG.itp"
#ifdef POSRES
#include "posre_LIG.itp"
#endif\n
#include "{2}.itp"
#ifdef POSRES_WATER
[ position_restraints ]
   1    1       1000       1000       1000
#endif\n
[ system ]
Protein-Ligand\n
[ molecules ]
mol1 1
LIG 1
""".format(initial, combine_atomtypes([lig_itp], ff), args.water)

	for output in ('bound', 'unbound'):
		try:
			os.mkdir(output)
		except FileExistsError:
			pass
		
		f = open("{0}/{1}".format(output, 'topol.top'), "w")
		f.write(template)
		f.close()
		
		mol_itp(lig_itp, 'LIG', "{0}/LIG.itp".format(output), ff)
		pp_itp(pro_itp, 'mol1', "{0}/mol1.itp".format(output))
		
		shutil.copyfile('posre_mol1.itp', "{0}/posre_mol1.itp".format(output))
		shutil.copyfile(os.path.dirname(os.path.abspath(__file__))+'/waters/itp/atomtypes_waters.itp',"{0}/atomtypes_waters.itp".format(output))
		shutil.copyfile(os.path.dirname(os.path.abspath(__file__))+'/waters/itp/{0}.itp'.format(args.water),"{0}/{1}.itp".format(output, args.water))
		if ff == 'gaff':
			shutil.copyfile('LIG.acpype/posre_LIG.itp', "{0}/posre_LIG.itp".format(output))
		elif ff == 'openff':
			with open(lig_gro) as gro:
				file_lines = gro.readlines()
				posreFile = open("{}/posre_LIG.itp".format(output), "w")
				posreFile.write("\n[ position_restraints ]\n; atom  type      fx      fy      fz\n")
				for line in file_lines[2:-2]:
					line_list = line.split()
					atom_name = str(line_list[1])
					if not atom_name.startswith("H"):
						posreFile.write("{0}     1 {1} {1} {1}\n".format(line_list[2], 1000))
				posreFile.close()
		else:
			shutil.copyfile('posre_LIG.itp', "{0}/posre_LIG.itp".format(output))
			shutil.copyfile("LIG_ffbonded.itp", "{0}/LIG_ffbonded.itp".format(output))
			
def combine_gro_files(mode, file1_path, file2_path, output_path, box='1 1 1'):
	with open(file1_path) as file1:
		file1_lines = file1.readlines()
	with open(file2_path) as file2:
		file2_lines = file2.readlines()
    
	num_atoms1 = int(file1_lines[1])
	num_atoms2 = int(file2_lines[1])
	total_atoms = num_atoms1 + num_atoms2
    
    # Create the header for the combined file
	title = "Combined Gro File"
	header = f"{title}\n{total_atoms}\n"
    
    # Combine the atom information from both files
	atoms = []
	if mode == 'HG':
		for i in range(2, num_atoms1 + 2):
			atoms.append(file1_lines[i][:5] + "HHH".ljust(5) + file1_lines[i][10:] )
		for i in range(2, num_atoms2 + 2):
			atoms.append(file2_lines[i][:5] + "GGG".ljust(5) + file2_lines[i][10:] )
	if mode == 'PP':
		for i in range(2, num_atoms1 + 2):
			atoms.append(file1_lines[i])
		for i in range(2, num_atoms2 + 2):
			atoms.append(file2_lines[i])
	if mode == 'PL':
		for i in range(2, num_atoms1 + 2):
			atoms.append(file1_lines[i])
		for i in range(2, num_atoms2 + 2):
			atoms.append(file2_lines[i][:5] + "LIG".ljust(5) + file2_lines[i][10:] )
		# Write the combined file
	with open(output_path, "w") as output_file:
		output_file.write(header)
		output_file.writelines(atoms)
		output_file.write(box)


def sphere_radius(filename):
    """
    Reads a GRO file and returns the coordinates of each atom as a numpy array.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()[2:-1]
        coords = np.zeros((len(lines), 3))
        for i, line in enumerate(lines):
            x, y, z = map(float, line[20:].split()[:3])
            coords[i, :] = x, y, z
    """
    Calculates the minimum sphere radius that covers the given set of coordinates,
    and the center of the sphere.
    """
    center = np.mean(coords, axis=0)
    radius = np.max(np.linalg.norm(coords - center, axis=1))
    return radius


def water_number():

	with open("unbound/topol.top", "r") as f1, open("bound/topol.top", "a") as f2:
		last_line = None
		for line in f1:
			last_line = line
		if last_line is not None:
			f2.write(last_line)
        
	with open("unbound/topol.top", "r") as f:
	    lines = f.readlines()

	last_line = lines[-1]
	molecules_line = last_line.split()  

	molecules_num = molecules_line[-1]  
	return molecules_num
	

def gmx_pdb2gmx(pdbfile, outcoord='protein.gro', outtop='topol.top', forcefield='amber99sb-ildn', water='tip3p', ignh=False, posre='posre.itp'):
    """Build a protein topology and coordinate file from a PDB file"""
    paras = {
        'gmx':GMX,
        'pdbfile':pdbfile,
        'outfile': outcoord,
        'topolfile': outtop,
        'forcefield': forcefield,
        'water': water,
        'ignh':ignh,
        'posre':posre}
    cmd = '{gmx} pdb2gmx -f {pdbfile} -o {outfile} -p {topolfile} -ff {forcefield} -water {water} -ignh {ignh} -i {posre} > gromacs.log 2>&1'.format(**paras)
    RC = os.system(cmd)
    if RC != 0:
        raise SystemExit('\nERROR! see the log file for details %s'%os.path.abspath("gromacs.log\n"))

def gmx_editconf(grofile, outcoord, center_x, center_y, center_z):
    """Build a protein topology and coordinate file from a PDB file"""
    paras = {
        'gmx':GMX,
        'grofile':grofile,
        'outfile': outcoord,
        'center_x': center_x,
        'center_y': center_y,
        'center_z': center_z}
    cmd = '{gmx} editconf -f {grofile} -o {outfile} -center {center_x} {center_y} {center_z} > gromacs.log 2>&1'.format(**paras)
    RC = os.system(cmd)
    if RC != 0:
        raise SystemExit('\nERROR! see the log file for details %s'%os.path.abspath("gromacs.log\n"))

def gmx_solvate(grofile, outcoord, topolfile, water):
    """Build a protein topology and coordinate file from a PDB file"""
    paras = {
        'gmx':GMX,
        'grofile':grofile,
        'outfile': outcoord,
        'topolfile': topolfile,
        'water': water}
    cmd = '{gmx} solvate -cp {grofile} -o {outfile} -p {topolfile} -cs {water}  > gromacs.log 2>&1'.format(**paras)
    RC = os.system(cmd)
    if RC != 0:
        raise SystemExit('\nERROR! see the log file for details %s'%os.path.abspath("gromacs.log\n"))



if __name__ == '__main__':

	args = ParserOptions()
	mode, mol1, mol2 = args.system
	
	if not os.path.isfile(mol1) or not os.path.isfile(mol2):
		raise SystemExit('\nERROR! {} or {} file missing! \n'.format(mol1, mol2))
	
	water_dir = os.path.dirname(os.path.abspath(__file__)) + '/waters/gro/'

	if args.water in water_models:
		water_model = water_models[args.water]
		WATER_1 = water_dir + '{}-1.gro'.format(water_model)
		WATER_S = water_dir + '{}.gro'.format(water_model)
	else:
		raise SystemExit('\nERROR!\nUnsupported water model')
	
	if args.userFF:
		if not os.path.isdir(args.proteinFF+'.ff') or not os.path.isdir(args.proteinFF+'.ff'):
			raise SystemExit('\nERROR! {} missing! \n'.format(args.proteinFF+'.ff'))

	### Setup for host-guest systems ###
	if mode == 'HG':
		
		if args.FF == 'openff':
			get_openff(mol1, 'HHH')
			get_openff(mol2, 'GGG')
			hg_topol('HHH.itp', 'GGG.itp', ff='openff')
			combine_gro_files(mode, "HHH.gro", "GGG.gro", "bound/input.gro")

			r1 = sphere_radius('HHH.gro')
			r2 = sphere_radius('GGG.gro')
			if r1 > r2:
				gmx_editconf('HHH.gro', 'unbound/H.gro', r1 + args.d, r1 + args.d, r1 + args.d)
				gmx_editconf('GGG.gro', 'unbound/G.gro', r2 + args.d + args.space + 2*r1, r1 + args.d, r1 + args.d)			
				combine_gro_files(mode, "unbound/H.gro", "unbound/G.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r1 + 2*args.d, 2*r1 + 2*args.d))
			else:
				gmx_editconf('HHH.gro', 'unbound/H.gro', r1 + args.space + args.d + 2*r2, r2 + args.d, r2 + args.d)
				gmx_editconf('GGG.gro', 'unbound/G.gro', r2 + args.d, r2 + args.d, r2 + args.d)	
				combine_gro_files(mode, "unbound/H.gro", "unbound/G.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r2 + 2*args.d, 2*r2 + 2*args.d))
			gmx_solvate('unbound/complex.gro', 'unbound/solv.gro', 'unbound/topol.top', WATER_S)
			os.system('''echo q|gmx make_ndx -f unbound/solv.gro -o unbound/index.ndx > gromacs.log 2>&1''')			
			RC = os.system('{0} editconf -f bound/input.gro -o bound/complex.gro -c -box {1} {1} {1} > gromacs.log 2>&1'.format(GMX, 2*r2 + args.space + 2*args.d + 2*r1))
			if RC != 0:
				raise SystemExit('\nERROR! see the log file for details %s'%os.path.abspath("gromacs.log\n"))			
			os.system('''gmx insert-molecules -f bound/complex.gro -nmol {0} -ci {1} -o bound/solv.gro > gromacs.log 2>&1'''.format(water_number(), WATER_1))
			os.system('''echo q|gmx make_ndx -f bound/solv.gro -o bound/index.ndx > gromacs.log 2>&1''')
		
		if args.FF in ('gaff', 'gaff2'):
			get_gaff(mol1, 'HHH', args.FF)	
			get_gaff(mol2, 'GGG', args.FF)
			hg_topol('HHH.acpype/HHH_GMX.itp', 'GGG.acpype/GGG_GMX.itp', ff='gaff')
			combine_gro_files(mode, "HHH.acpype/HHH_GMX.gro", "GGG.acpype/GGG_GMX.gro", "bound/input.gro")		
			
			r1 = sphere_radius('HHH.acpype/HHH_GMX.gro')
			r2 = sphere_radius('GGG.acpype/GGG_GMX.gro')
			if r1 > r2:
				gmx_editconf('HHH.acpype/HHH_GMX.gro', 'unbound/H.gro', r1 + args.d, r1 + args.d, r1 + args.d)
				gmx_editconf('GGG.acpype/GGG_GMX.gro', 'unbound/G.gro', r2 + args.d + args.space + 2*r1, r1 + args.d, r1 + args.d)			
				combine_gro_files(mode, "unbound/H.gro", "unbound/G.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r1 + 2*args.d, 2*r1 + 2*args.d))
			else:
				gmx_editconf('HHH.acpype/HHH_GMX.gro', 'unbound/H.gro', r1 + args.space + args.d + 2*r2, r2 + args.d, r2 + args.d)
				gmx_editconf('GGG.acpype/GGG_GMX.gro', 'unbound/G.gro', r2 + args.d, r2 + args.d, r2 + args.d)	
				combine_gro_files(mode, "unbound/H.gro", "unbound/G.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r2 + 2*args.d, 2*r2 + 2*args.d))
			gmx_solvate('unbound/complex.gro', 'unbound/solv.gro', 'unbound/topol.top', WATER_S)
			os.system('''echo q|gmx make_ndx -f unbound/solv.gro -o unbound/index.ndx > gromacs.log 2>&1''')			
			RC = os.system('{0} editconf -f bound/input.gro -o bound/complex.gro -c -box {1} {1} {1} > gromacs.log 2>&1'.format(GMX, 2*r2 + args.space + 2*args.d + 2*r1))
			if RC != 0:
				raise SystemExit('\nERROR! see the log file for details %s'%os.path.abspath("gromacs.log\n"))			
			os.system('''gmx insert-molecules -f bound/complex.gro -nmol {0} -ci {1} -o bound/solv.gro > gromacs.log 2>&1'''.format(water_number(), WATER_1))
			os.system('''echo q|gmx make_ndx -f bound/solv.gro -o bound/index.ndx > gromacs.log 2>&1''')
						
		if args.FF == 'cgenff':
			if 'SILCSBIODIR' not in os.environ or not os.environ['SILCSBIODIR']:
				raise ValueError("SILCSBIODIR environment variable is not set.")
			else:
				for name, mfile in zip(['HHH', 'GGG'], [mol1, mol2]):
					input_format = mfile.split('.')[-1].lower()

					if input_format != 'mol2':
						convert_to_mol2(mfile, name, name + 'ligand.mol2')
						get_cgenff(name + 'ligand.mol2', name)
					else:
						get_cgenff(mfile, name)
					
				hg_topol('HHH.itp', 'GGG.itp', ff='cgenff')
				combine_gro_files(mode, "HHH.gro", "GGG.gro", "bound/input.gro")
	
				r1 = sphere_radius('HHH.gro')
				r2 = sphere_radius('GGG.gro')
				if r1 > r2:
					gmx_editconf('HHH.gro', 'unbound/H.gro', r1 + args.d, r1 + args.d, r1 + args.d)
					gmx_editconf('GGG.gro', 'unbound/G.gro', r2 + args.d + args.space + 2*r1, r1 + args.d, r1 + args.d)			
					combine_gro_files(mode, "unbound/H.gro", "unbound/G.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r1 + 2*args.d, 2*r1 + 2*args.d))
				else:
					gmx_editconf('HHH.gro', 'unbound/H.gro', r1 + args.space + args.d + 2*r2, r2 + args.d, r2 + args.d)
					gmx_editconf('GGG.gro', 'unbound/G.gro', r2 + args.d, r2 + args.d, r2 + args.d)	
					combine_gro_files(mode, "unbound/H.gro", "unbound/G.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r2 + 2*args.d, 2*r2 + 2*args.d))
				gmx_solvate('unbound/complex.gro', 'unbound/solv.gro', 'unbound/topol.top', WATER_S)
				os.system('''echo q|gmx make_ndx -f unbound/solv.gro -o unbound/index.ndx > gromacs.log 2>&1''')			
				RC = os.system('{0} editconf -f bound/input.gro -o bound/complex.gro -c -box {1} {1} {1} > gromacs.log 2>&1'.format(GMX, 2*r2 + args.space + 2*args.d + 2*r1))
				if RC != 0:
					raise SystemExit('\nERROR! see the log file for details %s'%os.path.abspath("gromacs.log\n"))			
				os.system('''gmx insert-molecules -f bound/complex.gro -nmol {0} -ci {1} -o bound/solv.gro > gromacs.log 2>&1'''.format(water_number(), WATER_1))
				os.system('''echo q|gmx make_ndx -f bound/solv.gro -o bound/index.ndx > gromacs.log 2>&1''')
	
	### Setup for protein-protein systems ###
	if mode == 'PP':
	
		gmx_pdb2gmx(mol1, outcoord='mol1.gro', outtop='mol1.top', forcefield=args.proteinFF, water=args.water, ignh=False, posre='posre_mol1.itp')
		gmx_pdb2gmx(mol2, outcoord='mol2.gro', outtop='mol2.top', forcefield=args.proteinFF, water=args.water, ignh=False, posre='posre_mol2.itp')
		pp_topol('mol1.top', 'mol2.top',  ff=args.proteinFF, userff=args.userFF)
		combine_gro_files(mode, 'mol1.gro', 'mol2.gro', "bound/input.gro", box='1 1 1')
		
		r1 = sphere_radius('mol1.gro')
		r2 = sphere_radius('mol2.gro')
		if r1 > r2:
			gmx_editconf('mol1.gro', 'unbound/mol1.gro', r1 + args.d, r1 + args.d, r1 + args.d)
			gmx_editconf('mol2.gro', 'unbound/mol2.gro', r2 + args.d + args.space + 2*r1, r1 + args.d, r1 + args.d)			
			combine_gro_files(mode, "unbound/mol1.gro", "unbound/mol2.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r1 + 2*args.d, 2*r1 + 2*args.d))
		else:
			gmx_editconf('mol1.gro', 'unbound/mol1.gro', r1 + args.space + args.d + 2*r2, r2 + args.d, r2 + args.d)
			gmx_editconf('mol2.gro', 'unbound/mol2.gro', r2 + args.d, r2 + args.d, r2 + args.d)	
			combine_gro_files(mode, "unbound/mol1.gro", "unbound/mol2.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r2 + 2*args.d, 2*r2 + 2*args.d))
		gmx_solvate('unbound/complex.gro', 'unbound/solv.gro', 'unbound/topol.top', WATER_S)
		os.system('''echo q|gmx make_ndx -f unbound/solv.gro -o unbound/index.ndx > gromacs.log 2>&1''')			
		RC = os.system('{0} editconf -f bound/input.gro -o bound/complex.gro -c -box {1} {1} {1} > gromacs.log 2>&1'.format(GMX, 2*r2 + args.space + 2*args.d + 2*r1))
		if RC != 0:
			raise SystemExit('\nERROR! see the log file for details %s'%os.path.abspath("gromacs.log\n"))			
		os.system('''gmx insert-molecules -f bound/complex.gro -nmol {0} -ci {1} -o bound/solv.gro > gromacs.log 2>&1'''.format(water_number(), WATER_1))
		os.system('''echo q|gmx make_ndx -f bound/solv.gro -o bound/index.ndx > gromacs.log 2>&1''')
		
	### Setup for protein-ligand systems ###
	if mode == 'PL':

		if args.FF == 'openff':
			gmx_pdb2gmx(mol1, outcoord='mol1.gro', outtop='mol1.top', forcefield=args.proteinFF, water=args.water, ignh=False, posre='posre_mol1.itp')
			get_openff(mol2, 'LIG')
			pl_topol('mol1.top', 'LIG.itp', 'LIG.gro', ff=args.FF, proteinFF=args.proteinFF, userff=args.userFF)
			combine_gro_files(mode, 'mol1.gro', 'LIG.gro', "bound/input.gro", box='1 1 1')
			r1 = sphere_radius('mol1.gro')
			r2 = sphere_radius('LIG.gro')
			if r1 > r2:
				gmx_editconf('mol1.gro', 'unbound/P.gro', r1 + args.d, r1 + args.d, r1 + args.d)
				gmx_editconf('LIG.gro', 'unbound/L.gro', r2 + args.d + args.space + 2*r1, r1 + args.d, r1 + args.d)			
				combine_gro_files(mode, "unbound/P.gro", "unbound/L.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r1 + 2*args.d, 2*r1 + 2*args.d))
			else:
				gmx_editconf('mol1.gro', 'unbound/P.gro', r1 + args.space + args.d + 2*r2, r2 + args.d, r2 + args.d)
				gmx_editconf('LIG.gro', 'unbound/L.gro', r2 + args.d, r2 + args.d, r2 + args.d)	
				combine_gro_files(mode, "unbound/P.gro", "unbound/L.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r2 + 2*args.d, 2*r2 + 2*args.d))
			gmx_solvate('unbound/complex.gro', 'unbound/solv.gro', 'unbound/topol.top', WATER_S)
			os.system('''echo q|gmx make_ndx -f unbound/solv.gro -o unbound/index.ndx > gromacs.log 2>&1''')			
			RC = os.system('{0} editconf -f bound/input.gro -o bound/complex.gro -c -box {1} {1} {1} > gromacs.log 2>&1'.format(GMX, 2*r2 + args.space + 2*args.d + 2*r1))
			if RC != 0:
				raise SystemExit('\nERROR! see the log file for details %s'%os.path.abspath("gromacs.log\n"))			
			os.system('''gmx insert-molecules -f bound/complex.gro -nmol {0} -ci {1} -o bound/solv.gro > gromacs.log 2>&1'''.format(water_number(), WATER_1))
			os.system('''echo q|gmx make_ndx -f bound/solv.gro -o bound/index.ndx > gromacs.log 2>&1''')
						
		
		if args.FF in ('gaff', 'gaff2'):
			gmx_pdb2gmx(mol1, outcoord='mol1.gro', outtop='mol1.top', forcefield=args.proteinFF, water=args.water, ignh=False, posre='posre_mol1.itp')
			get_gaff(mol2, 'LIG', args.FF)
			pl_topol('mol1.top', 'LIG.acpype/LIG_GMX.itp', 'LIG.acpype/LIG_GMX.gro', ff=args.FF, proteinFF=args.proteinFF, userff=args.userFF)
			combine_gro_files(mode, 'mol1.gro', 'LIG.acpype/LIG_GMX.gro', "bound/input.gro", box='1 1 1')
			r1 = sphere_radius('mol1.gro')
			r2 = sphere_radius('LIG.acpype/LIG_GMX.gro')
			if r1 > r2:
				gmx_editconf('mol1.gro', 'unbound/P.gro', r1 + args.d, r1 + args.d, r1 + args.d)
				gmx_editconf('LIG.acpype/LIG_GMX.gro', 'unbound/L.gro', r2 + args.d + args.space + 2*r1, r1 + args.d, r1 + args.d)			
				combine_gro_files(mode, "unbound/P.gro", "unbound/L.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r1 + 2*args.d, 2*r1 + 2*args.d))
			else:
				gmx_editconf('mol1.gro', 'unbound/P.gro', r1 + args.space + args.d + 2*r2, r2 + args.d, r2 + args.d)
				gmx_editconf('LIG.acpype/LIG_GMX.gro', 'unbound/L.gro', r2 + args.d, r2 + args.d, r2 + args.d)	
				combine_gro_files(mode, "unbound/P.gro", "unbound/L.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r2 + 2*args.d, 2*r2 + 2*args.d))
			gmx_solvate('unbound/complex.gro', 'unbound/solv.gro', 'unbound/topol.top', WATER_S)
			os.system('''echo q|gmx make_ndx -f unbound/solv.gro -o unbound/index.ndx > gromacs.log 2>&1''')			
			RC = os.system('{0} editconf -f bound/input.gro -o bound/complex.gro -c -box {1} {1} {1} > gromacs.log 2>&1'.format(GMX, 2*r2 + args.space + 2*args.d + 2*r1))
			if RC != 0:
				raise SystemExit('\nERROR! see the log file for details %s'%os.path.abspath("gromacs.log\n"))			
			os.system('''gmx insert-molecules -f bound/complex.gro -nmol {0} -ci {1} -o bound/solv.gro > gromacs.log 2>&1'''.format(water_number(), WATER_1))
			os.system('''echo q|gmx make_ndx -f bound/solv.gro -o bound/index.ndx > gromacs.log 2>&1''')

							
		if args.FF == 'cgenff':
			if 'SILCSBIODIR' not in os.environ or not os.environ['SILCSBIODIR']:
				raise ValueError("SILCSBIODIR environment variable is not set.")
			else:
				input_format = mol2.split('.')[-1].lower()

				if input_format != 'mol2':
					convert_to_mol2(mol2, 'LIG', 'ligand.mol2')
					get_cgenff('ligand.mol2', 'LIG')
				else:
					get_cgenff(mol2, 'LIG')

				if args.proteinFF == 'charmm36':
					if not os.path.exists('charmm36.ff'):
						shutil.copytree(os.environ.get('SILCSBIODIR') + '/data/gromacs/charmm36.ff', 'charmm36.ff')
					gmx_pdb2gmx(mol1, outcoord='mol1.gro', outtop='mol1.top', forcefield=args.proteinFF, water=args.water, ignh=False, posre='posre_mol1.itp')
					
				if args.proteinFF == 'charmm27':
					gmx_pdb2gmx(mol1, outcoord='mol1.gro', outtop='mol1.top', forcefield=args.proteinFF, water=args.water, ignh=False, posre='posre_mol1.itp')

				if args.proteinFF not in ['charmm27', 'charmm36']:
					print('\nWarning! {} is not compatible with {} !\n'.format(args.FF, args.proteinFF))
					gmx_pdb2gmx(mol1, outcoord='mol1.gro', outtop='mol1.top', forcefield=args.proteinFF, water=args.water, ignh=False, posre='posre_mol1.itp')
									
				pl_topol('mol1.top', 'LIG.itp', 'LIG.gro', ff=args.FF, proteinFF=args.proteinFF, userff=args.userFF)
				combine_gro_files(mode, 'mol1.gro', 'LIG.gro', "bound/input.gro", box='1 1 1')

				r1 = sphere_radius('mol1.gro')
				r2 = sphere_radius('LIG.gro')				
				if r1 > r2:
					gmx_editconf('mol1.gro', 'unbound/P.gro', r1 + args.d, r1 + args.d, r1 + args.d)
					gmx_editconf('LIG.gro', 'unbound/L.gro', r2 + args.d + args.space + 2*r1, r1 + args.d, r1 + args.d)			
					combine_gro_files(mode, "unbound/P.gro", "unbound/L.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r1 + 2*args.d, 2*r1 + 2*args.d))
				else:
					gmx_editconf('mol1.gro', 'unbound/P.gro', r1 + args.space + args.d + 2*r2, r2 + args.d, r2 + args.d)
					gmx_editconf('LIG.gro', 'unbound/L.gro', r2 + args.d, r2 + args.d, r2 + args.d)	
					combine_gro_files(mode, "unbound/P.gro", "unbound/L.gro", "unbound/complex.gro",'{0} {1} {2}'.format(2*r2 + args.space + 2*args.d + 2*r1, 2*r2 + 2*args.d, 2*r2 + 2*args.d))
				gmx_solvate('unbound/complex.gro', 'unbound/solv.gro', 'unbound/topol.top', WATER_S)
				os.system('''echo q|gmx make_ndx -f unbound/solv.gro -o unbound/index.ndx > gromacs.log 2>&1''')			
				RC = os.system('{0} editconf -f bound/input.gro -o bound/complex.gro -c -box {1} {1} {1} > gromacs.log 2>&1'.format(GMX, 2*r2 + args.space + 2*args.d + 2*r1))
				if RC != 0:
					raise SystemExit('\nERROR! see the log file for details %s'%os.path.abspath("gromacs.log\n"))			
				os.system('''gmx insert-molecules -f bound/complex.gro -nmol {0} -ci {1} -o bound/solv.gro > gromacs.log 2>&1'''.format(water_number(), WATER_1))
				os.system('''echo q|gmx make_ndx -f bound/solv.gro -o bound/index.ndx > gromacs.log 2>&1''')



