import sys

d = sys._MEIPASS

if sys.platform.startswith('linux'):
	os.environ['BABEL_DATADIR'] = os.path.join(d, 'openbabel', 'data')
	os.environ['BABEL_LIBDIR'] = os.path.join(d, 'openbabel', 'lib')
