import os
import time
from openbabel import openbabel

__all__ = ['AttrDict', 'draw_gridbox', 'convert_dimension_to_coordinates',
	'convert_coordinates_to_dimension', 'get_atom_types_from_pdbqt',
	'get_molecule_center_from_pdbqt', 'time_format', 'convert_pdbqt_to_pdb'
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


if __name__ == '__main__':
	import sys
	infile = sys.argv[1]
	outfile = infile.replace('pdb', 'pdbqt')
	convert_other_to_pdbqt(infile, 'pdb', outfile)
