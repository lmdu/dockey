import os
import time
import math
from openbabel import openbabel

__all__ = ['AttrDict', 'draw_gridbox', 'convert_dimension_to_coordinates',
	'convert_coordinates_to_dimension', 'get_atom_types_from_pdbqt',
	'get_molecule_center_from_pdbqt', 'time_format', 'convert_pdbqt_to_pdb',
	'ligand_efficiency_assessment', 'get_molecule_information'
]

class AttrDict(dict):
	def __getattr__(self, attr):
		try:
			return self[attr]
		except KeyError:
			raise Exception("no attribute")

	def __setattr__(self, attr, val):
		self[attr] = val

def convert_dimension_to_coordinates(x, y, z, cx, cy, cz, spacing):
	spacing = float(spacing)
	
	x_size = x * spacing
	y_size = y * spacing
	z_size = z * spacing

	x_half = x_size/2
	y_half = y_size/2
	z_half = z_size/2

	min_points = [cx-x_half, cy-y_half, cz-z_half]
	max_points = [cx+x_half, cy+y_half, cz+z_half]

	return (min_points, max_points)

def convert_coordinates_to_dimension(points, spacing):
	[(xmin, ymin, zmin), (xmax, ymax, zmax)] = points
	x_size = xmax - xmin
	y_size = ymax - ymin
	z_size = zmax - zmin

	x = int(x_size/spacing)
	y = int(y_size/spacing)
	z = int(z_size/spacing)

	if x % 2 != 0:
		x -= 1

	if y % 2 != 0:
		y -= 1

	if z % 2 != 0:
		z -= 1

	cx = xmin + x_size/2
	cy = ymin + y_size/2
	cz = zmin + z_size/2

	return (x, y, z, cx, cy, cz)

def draw_gridbox(cmd, data):
	cmd.draw_box(
		points = convert_dimension_to_coordinates(
			data.x, data.y, data.z,
			data.cx, data.cy, data.cz,
			data.spacing
		),
		show_face = data.show_face,
		show_edge = data.show_edge,
		use_line = data.use_line,
		use_cylinder = data.use_cylinder,
		opacity = data.opacity,
		edge_width = data.edge_width,
		edge_color = data.edge_color,
		bg_x = data.bg_x,
		bg_y = data.bg_y,
		bg_z = data.bg_z
	)

def convert_pdbqt_to_pdb(pdbqt_str):
	obc = openbabel.OBConversion()
	obc.SetInAndOutFormats('pdbqt', 'pdb')
	mol = openbabel.OBMol()
	obc.ReadString(mol, pdbqt_str)
	return obc.WriteString(mol)

def convert_other_to_pdbqt(infile, informat, outfile):
	obc = openbabel.OBConversion()
	obc.SetInAndOutFormats(informat, "pdbqt")
	mol = openbabel.OBMol()
	obc.ReadFile(mol, infile)

	non_bond = 0

	for atom in openbabel.OBMolAtomIter(mol):
		flag = 0
		for bond in openbabel.OBAtomBondIter(atom):
			flag = 1
			break

		if flag == 0:
			non_bond += 1
			

	print(non_bond)


	#correct for ph
	mol.CorrectForPH()

	#add hydrogens
	mol.AddHydrogens()
	#mol.AddPolarHydrogens()

	cmodel = openbabel.OBChargeModel.FindType('gasteiger')
	cmodel.ComputeCharges(mol)

	#set to be rigid
	obc.AddOption('r')
	obc.AddOption('x')
	#obc.AddOption('h')
	obc.AddOption('b')
	obc.AddOption('p')

	obc.WriteFile(mol, outfile)

def get_atom_types_from_pdbqt(pdbqt_file):
	atom_types = set()
	with open(pdbqt_file) as fh:
		for line in fh:
			if not line.startswith(('ATOM', 'HETATM')):
				continue

			cols = line.strip().split()

			atom_types.add(cols[-1])

	return sorted(atom_types)

def get_molecule_center_from_pdbqt(pdbqt_file):
	c = 0
	x = 0
	y = 0
	z = 0
	with open(pdbqt_file) as fh:
		for line in fh:
			if not line.startswith(('ATOM', 'HETATM')):
				continue

			cols = line.strip().split()
			c += 1 
			x += float(cols[6])
			y += float(cols[7])
			z += float(cols[8])

	x = round(x/c, 3)
	y = round(y/c, 3)
	z = round(z/c, 3)

	return x, y, z

def time_format(seconds):
	if seconds:
		t = time.localtime(seconds)
		return time.strftime("%Y-%m-%d %H:%M:%S", t)

def convert_ki_to_log(ki, unit):
	scales = {
		'M': 0, 'mM': 3, 'µM': 6, 'nM': 9, 'pM': 12,
		'fM': 15, 'aM': 18, 'zM': 21, 'yM': 24
	}
	ki = ki/math.pow(10, scales[unit])
	return round(math.log10(ki), 3)

def convert_energy_to_ki(energy):
	R = 0.001987207
	T = 298.15
	ki = math.exp(energy/(R*T))
	return round(math.log10(ki), 3)

def get_molecule_information(mol_file):
	mol_name, mol_format = os.path.splitext(os.path.basename(mol_file))
	mol_format = mol_format.lstrip('.')

	obc = openbabel.OBConversion()
	obc.SetInAndOutFormats(mol_format, 'pdb')
	mol = openbabel.OBMol()
	obc.ReadFile(mol, mol_file)
	des = openbabel.OBDescriptor.FindType('logP')

	return AttrDict(
		name = mol_name,
		atoms = mol.NumAtoms(),
		bonds = mol.NumBonds(),
		hvyatoms = mol.NumHvyAtoms(),
		residues = mol.NumResidues(),
		rotors = mol.NumRotors(),
		formula = mol.GetFormula(),
		energy = mol.GetEnergy(),
		weight = mol.GetMolWt(),
		logp = des.Predict(mol),
		pdb = obc.WriteString(mol)
	)

def ligand_efficiency_assessment(pdb_str, energy, ki=0):
	if ki:
		logki = convert_ki_to_log(ki)
	else:
		logki = convert_energy_to_ki(energy)

	if pdb_str:
		obc = openbabel.OBConversion()
		obc.SetInFormat("pdb")
		mol = openbabel.OBMol()
		obc.ReadString(mol, pdb_str)
		ha = mol.NumHvyAtoms()
		mw = mol.GetMolWt()

		#calculate ligand efficiency
		le = round(energy/ha*-1, 3)

		#calculate logP
		des = openbabel.OBDescriptor.FindType('logP')
		logp = des.Predict(mol)

		#calculate ligand lipophilic efficiency
		lle = -1*logki - logp

		#calculate fit quality
		le_scale = 0.0715 + 7.5328/ha + 25.7079/math.pow(ha,2) - 361.4722/math.pow(ha, 3)
		fq = round(le/le_scale, 3)

		#calculate ligand efficiency lipophilic price
		lelp = round(logp/le, 3)

		return (logki, le, lle, fq, lelp)
	else:
		return (logki, None, None, None, None)

if __name__ == '__main__':
	import sys
	infile = sys.argv[1]
	outfile = infile.replace('pdb', 'pdbqt')
	convert_other_to_pdbqt(infile, 'pdb', outfile)
