import os
import sys

d = sys._MEIPASS
is_win = sys.platform.startswith('win')

if is_win:
	os.environ['BABEL_DATADIR'] = os.path.join(d, 'data')
	os.environ['BABEL_LIBDIR'] = os.path.join(d)
else:
	os.environ['BABEL_DATADIR'] = os.path.join(d, 'openbabel', 'data')
	os.environ['BABEL_LIBDIR'] = os.path.join(d, 'openbabel', 'plugin')
