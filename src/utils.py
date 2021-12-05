import os
from openbabel import openbabel

__all__ = ['AttrDict', 'draw_gridbox', 'convert_dimension_to_coordinates',
	'convert_coordinates_to_dimension', 'convert_other_to_pdbqt'
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

if __name__ == '__main__':
	import sys
	infile = sys.argv[1]
	outfile = infile.replace('pdb', 'pdbqt')
	convert_other_to_pdbqt(infile, 'pdb', outfile)
