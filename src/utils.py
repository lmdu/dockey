import os
from openbabel import openbabel

__all__ = ['AttrDict', 'draw_gridbox', 'convert_dimension_to_coordinates',
	'convert_coordinates_to_dimension'
]


class AttrDict(dict):
	def __getattr__(self, attr):
		try:
			return self[attr]
		except KeyError:
			raise Exception("no attribute")

	def __setattr__(self, attr, val):
		self[attr] = val

def convert_dimension_to_coordinates(x_size, y_size, z_size, x_center, y_center, z_center):
	x_half = x_size/2
	y_half = y_size/2
	z_half = z_size/2

	min_point = [x_center-x_half, y_center-y_half, z_center-z_half]
	max_point = [x_center+x_half, y_center+y_half, z_center+z_half]

	return (min_point, max_point)

def convert_coordinates_to_dimension(points):
	[(xmin, ymin, zmin), (xmax, ymax, zmax)] = points
	x_size = xmax - xmin
	y_size = ymax - ymin
	z_size = zmax - zmin

	x_center = xmin + x_size/2
	y_center = ymin + y_size/2
	z_center = zmin + z_size/2

	return (x_size, y_size, z_size, x_center, y_center, z_center)

def draw_gridbox(cmd, data):
	cmd.draw_box(
		points = convert_dimension_to_coordinates(
			data.x, data.y, data.z,
			data.cx, data.cy, data.cz
		),
		show_face = data.show_face,
		show_edge = data.show_edge,
		use_line = data.use_line,
		use_cylinder = data.use_cylinder,
		padding = data.padding,
		opacity = data.opacity,
		edge_width = data.edge_width,
		edge_color = data.edge_color,
		bg_x = data.bg_x,
		bg_y = data.bg_y,
		bg_z = data.bg_z
	)

def get_atom_types_from_pdbqt(pdbqt_file):
	atom_types = set()
	with open(pdbqt_file) as fh:
		for line in fh:
			if not line.startswith(('ATOM', 'HETATM')):
				continue

			atom_types.add(line[78:80].strip())

	return sorted(atom_types)

def convert_other_to_pdbqt(other_file):
	fname, fext = os.path.splitext(other_file)
	obc = openbabel.OBConversion()
	obc.SetInAndOutFormats(fext, "pdbqt")
	new_file = ''.join([fname, '.pdbqt'])

	mol = openbabel.OBMol()
	obc.ReadFile(mol, other_file)
	obc.WriteFile(mol, new_file)
