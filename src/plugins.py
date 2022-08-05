from pymol import cmd
from pymol.cgo import *

def register_plugin(func):
	name = func.__name__

	if hasattr(cmd, name):
		func.__doc__ = getattr(cmd, name).__doc__

	setattr(cmd, name, func)
	cmd.extend(func)
	return func

@register_plugin
def draw_box(points=([0,0,0],[5,5,5]), show_face=True, bg_x=(1.0,0.0,0.0), bg_y=(0.0,1.0,0.0),
			 bg_z=(0.0,0.0,1.0), show_edge=True, use_line=True, use_cylinder=False, radius=0.2,
			 edge_width=2, edge_color=(1.0,1.0,1.0), opacity=50.0):

	([xmin, ymin, zmin], [xmax, ymax, zmax]) = points

	cmd.delete('gridbox')
	view = cmd.get_view()

	opacity = float(opacity)/100.0
	edge_width = float(edge_width)

	xmin = float(xmin)
	ymin = float(ymin)
	zmin = float(zmin)
	xmax = float(xmax)
	ymax = float(ymax)
	zmax = float(zmax)

	obj = []

	if show_face:
		obj.extend([
			ALPHA, opacity,
			BEGIN, TRIANGLE_STRIP,
			COLOR, bg_x[0], bg_x[1], bg_x[2],
			VERTEX, xmin, ymin, zmin, #1
			VERTEX, xmin, ymin, zmax, #2
			VERTEX, xmin, ymax, zmin, #3
			VERTEX, xmin, ymax, zmax, #4
			END,

			BEGIN, TRIANGLE_STRIP,
			COLOR, bg_x[0], bg_x[1], bg_x[2],
			VERTEX, xmax, ymin, zmin, #5
			VERTEX, xmax, ymin, zmax, #6
			VERTEX, xmax, ymax, zmin, #7
			VERTEX, xmax, ymax, zmax, #8
			END,

			BEGIN, TRIANGLE_STRIP,
			COLOR, bg_y[0], bg_y[1], bg_y[2],
			VERTEX, xmin, ymin, zmin, #1
			VERTEX, xmin, ymin, zmax, #2
			VERTEX, xmax, ymin, zmin, #5
			VERTEX, xmax, ymin, zmax, #6 
			END,

			BEGIN, TRIANGLE_STRIP,
			COLOR, bg_y[0], bg_y[1], bg_y[2],
			VERTEX, xmin, ymax, zmin, #3
			VERTEX, xmin, ymax, zmax, #4
			VERTEX, xmax, ymax, zmin, #7
			VERTEX, xmax, ymax, zmax, #8
			END,

			BEGIN, TRIANGLE_STRIP,
			COLOR, bg_z[0], bg_z[1], bg_z[2],
			VERTEX, xmin, ymin, zmin, #1
			VERTEX, xmin, ymax, zmin, #3
			VERTEX, xmax, ymin, zmin, #5
			VERTEX, xmax, ymax, zmin, #7
			END,

			BEGIN, TRIANGLE_STRIP,
			COLOR, bg_z[0], bg_z[1], bg_z[2],
			VERTEX, xmin, ymin, zmax, #2
			VERTEX, xmin, ymax, zmax, #4
			VERTEX, xmax, ymin, zmax, #6
			VERTEX, xmax, ymax, zmax, #8
			END,
		])

	if show_edge and use_line:
		obj.extend([
			ALPHA, 1,
			LINEWIDTH, edge_width,

			BEGIN, LINES,
			COLOR, edge_color[0], edge_color[1], edge_color[2],

			VERTEX, xmin, ymin, zmin, #1
			VERTEX, xmin, ymin, zmax, #2

			VERTEX, xmin, ymax, zmin, #3
			VERTEX, xmin, ymax, zmax, #4

			VERTEX, xmax, ymin, zmin, #5
			VERTEX, xmax, ymin, zmax, #6

			VERTEX, xmax, ymax, zmin, #7
			VERTEX, xmax, ymax, zmax, #8

			VERTEX, xmin, ymin, zmin, #1
			VERTEX, xmax, ymin, zmin, #5

			VERTEX, xmin, ymax, zmin, #3
			VERTEX, xmax, ymax, zmin, #7

			VERTEX, xmin, ymax, zmax, #4
			VERTEX, xmax, ymax, zmax, #8

			VERTEX, xmin, ymin, zmax, #2
			VERTEX, xmax, ymin, zmax, #6

			VERTEX, xmin, ymin, zmin, #1
			VERTEX, xmin, ymax, zmin, #3

			VERTEX, xmax, ymin, zmin, #5
			VERTEX, xmax, ymax, zmin, #7

			VERTEX, xmin, ymin, zmax, #2
			VERTEX, xmin, ymax, zmax, #4

			VERTEX, xmax, ymin, zmax, #6
			VERTEX, xmax, ymax, zmax, #8

			END,
		])

	if show_edge and use_cylinder:
		c = edge_color
		radius = float(edge_width)/10

		obj.extend([
			CYLINDER, xmin, ymin, zmin, xmin, ymin, zmax, radius, c[0], c[1], c[2], c[0], c[1], c[2],
			CYLINDER, xmin, ymax, zmin, xmin, ymax, zmax, radius, c[0], c[1], c[2], c[0], c[1], c[2],
			CYLINDER, xmax, ymin, zmin, xmax, ymin, zmax, radius, c[0], c[1], c[2], c[0], c[1], c[2],
			CYLINDER, xmax, ymax, zmin, xmax, ymax, zmax, radius, c[0], c[1], c[2], c[0], c[1], c[2],
			CYLINDER, xmin, ymin, zmin, xmax, ymin, zmin, radius, c[0], c[1], c[2], c[0], c[1], c[2],
			CYLINDER, xmin, ymax, zmin, xmax, ymax, zmin, radius, c[0], c[1], c[2], c[0], c[1], c[2],
			CYLINDER, xmin, ymax, zmax, xmax, ymax, zmax, radius, c[0], c[1], c[2], c[0], c[1], c[2],
			CYLINDER, xmin, ymin, zmax, xmax, ymin, zmax, radius, c[0], c[1], c[2], c[0], c[1], c[2],
			CYLINDER, xmin, ymin, zmin, xmin, ymax, zmin, radius, c[0], c[1], c[2], c[0], c[1], c[2],
			CYLINDER, xmax, ymin, zmin, xmax, ymax, zmin, radius, c[0], c[1], c[2], c[0], c[1], c[2],
			CYLINDER, xmin, ymin, zmax, xmin, ymax, zmax, radius, c[0], c[1], c[2], c[0], c[1], c[2],
			CYLINDER, xmax, ymin, zmax, xmax, ymax, zmax, radius, c[0], c[1], c[2], c[0], c[1], c[2]
		])

	#auto_zoom = cmd.get('auto_zoom')
	#cmd.set('auto_zoom', 0, quiet=1)
	cmd.load_cgo(obj, 'gridbox')
	#cmd.set('auto_zoom', auto_zoom, quiet=1)
	cmd.set("cgo_line_width", edge_width)
	cmd.set_view(view)

	#return 'gridbox'

@register_plugin
def get_box(name="gridbox"):
	((xmin, ymin, zmin),(xmax, ymax, zmax)) = cmd.get_extent(name)



#pymol.cmd.extend('draw_box', draw_box)
#pymol.cmd.extend('get_box', get_box)
#from random import randint

@register_plugin
def draw_bounding_box(selection="(all)"):
	draw_box(cmd.get_extent(selection))

@register_plugin
def initialize():
	val = cmd.get('internal_gui')
	cmd.reinitialize()
	cmd.set('internal_gui', val)
       