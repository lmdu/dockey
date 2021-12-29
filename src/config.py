import sys
import pymol
import PySide6
import openbabel

__all__ = ['DOCKEY_VERSION', 'DOCKEY_BUILD', 'DOCKEY_ABOUT']

DOCKEY_VERSION = "0.1.0"

DOCKEY_BUILD = "211118"

DOCKEY_ABOUT = """
<p>Dockey - Molecular Docking and Virtual Screening</p>
<p><b>Version</b> v{version} <b>Build</b> {build}</p>
<p>Dockey is a morden tool for molecular docking</p>
<p><b>Acknowledgements:</b></p>
<table cellpadding="0" cellspacing="5">
	<tr>
		<td>Python</td>
		<td>v{python}</td>
		<td>
			<a href="https://www.python.org/">
				https://www.python.org
			</a>
		</td>
	</tr>
	<tr>
		<td>PySide</td>
		<td>v{pyside}</td>
		<td>
			<a href="https://doc.qt.io/qtforpython/">
				https://doc.qt.io/qtforpython
			</a>
		</td>
	</tr>
	<tr>
		<td>Pymol</td>
		<td>v{pymol}</td>
		<td>
			<a href="https://pymol.org">
				https://pymol.org
			</a>
		</td>
	</tr>
	<tr>
		<td>OpenBabel</td>
		<td>v{babel}</td>
		<td>
			<a href="http://openbabel.org">
				http://openbabel.org
			</a>
		</td>
	</tr>
</table>
<p><b>Requirements:</b></p>
<table cellpadding="0" cellspacing="5">
	<tr>
		<td>MGLTools</td>
		<td>v1.5.7</td>
		<td>
			<a href="https://ccsb.scripps.edu/mgltools/">
				https://ccsb.scripps.edu/mgltools
			</a>
		</td>
	</tr>
	<tr>
		<td>AutoDock</td>
		<td>v4.2.6</td>
		<td>
			<a href="https://autodock.scripps.edu/">
				https://autodock.scripps.edu
			</a>
		</td>
	</tr>
	<tr>
		<td>AutoDock Vina</td>
		<td>v1.2.3</td>
		<td>
			<a href="https://github.com/ccsb-scripps/AutoDock-Vina/">
				https://github.com/ccsb-scripps/AutoDock-Vina
			</a>
		</td>
	</tr>
</table>
""".format(
	version = DOCKEY_VERSION,
	build = DOCKEY_BUILD,
	python = sys.version.split()[0],
	pyside = PySide6.__version__,
	pymol = pymol.cmd.get_version()[0],
	babel = openbabel.__version__
)
