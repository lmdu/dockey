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
<table cellspacing="5">
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
		<td>PySide6</td>
		<td>v{pyside}</td>
		<td>
			<a href="https://doc.qt.io/qtforpython/">
				https://doc.qt.io/qtforpython
			</a>
		<td>
	</tr>
	<tr>
		<td>Pymol</td>
		<td>v{pymol}</td>
		<td>
			<a href="https://pymol.org">
				https://pymol.org
			</a>
		<td>
	</tr>
	<tr>
		<td>OpenBabel</td>
		<td>v{babel}</td>
		<td>
			<a href="http://openbabel.org">
				http://openbabel.org
			</a>
		<td>
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
