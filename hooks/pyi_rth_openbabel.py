import os
import sys

d = sys._MEIPASS

if sys.platform.startswith('win32'):
	os.environ['BABEL_DATADIR'] = os.path.join(d, 'openbabel', 'data')
	os.environ['BABEL_LIBDIR'] = os.path.join(d, 'openbabel')
else:
	os.environ['BABEL_DATADIR'] = os.path.join(d, 'openbabel', 'data')
	os.environ['BABEL_LIBDIR'] = os.path.join(d, 'openbabel', 'lib')
