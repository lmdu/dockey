import sys
import apsw
import numpy
import rdkit
import pymol
import meeko
import openmm
import psutil
import pdb2pqr
import requests
import pdbfixer
import openbabel
from PySide6.QtCore import __version__ as pyside_version
from plip.basic.config import __version__ as plip_version

__all__ = ['DOCKEY_VERSION', 'DOCKEY_BUILD', 'DOCKEY_ABOUT',
			'DOCKEY_THANKS']

DOCKEY_VERSION = "0.10.0"

DOCKEY_BUILD = "231110"

DOCKEY_ABOUT = """
<p>Dockey - Molecular Docking and Virtual Screening</p>
<p><b>Version</b> v{version} <b>Build</b> {build}</p>
<p>Dockey is a morden graphical user interface tool for
molecular docking. Dockey seamlessly integrates serveral
external tools to implement a complete streamlined docking
pipeline including molecular preprocessing, molecular preparation,
paralleled docking execution, interaction detection and
conformation visualization.
<p><b>Citation:</b><br>
Du L, Geng C, Zeng Q et al. Dockey: a modern integrated tool for
large-scale molecular docking and virtual screening. <i>Briefings in
Bioinformatics</i>. 2023, 24(2):bbad047.
(doi:<a href="https://doi.org/10.1093/bib/bbad047">10.1093/bib/bbad047</a>)
</p>
""".format(
	version = DOCKEY_VERSION,
	build = DOCKEY_BUILD
)

DOCKEY_THANKS = """
<table cellpadding="0" cellspacing="15">
	<tr>
		<td>Python</td>
		<td>v{python}</td>
		<td>
			<a href="https://www.python.org/">
				https://www.python.org
			</a>
		</td>
		<td>PSFL</td>
	</tr>
	<tr>
		<td>APSW</td>
		<td>v{apsw}</td>
		<td>
			<a href="https://github.com/rogerbinns/apsw/">
				https://github.com/rogerbinns/apsw
			</a>
		</td>
		<td>OSI Approved</td>
	</tr>
	<tr>
		<td>PySide6</td>
		<td>v{pyside}</td>
		<td>
			<a href="https://doc.qt.io/qtforpython-6/">
				https://doc.qt.io/qtforpython-6
			</a>
		</td>
		<td>LGPL</td>
	</tr>
	<tr>
		<td>Pymol</td>
		<td>v{pymol}</td>
		<td>
			<a href="https://pymol.org">
				https://pymol.org
			</a>
		</td>
		<td>BSD-like</td>
	</tr>
	<tr>
		<td>OpenBabel</td>
		<td>v{babel}</td>
		<td>
			<a href="http://openbabel.org">
				http://openbabel.org
			</a>
		</td>
		<td>GPL-2.0</td>
	</tr>
		<tr>
		<td>RDKit</td>
		<td>v{rdkit}</td>
		<td>
			<a href="https://www.rdkit.org/">
				https://www.rdkit.org
			</a>
		</td>
		<td>BSD 3-Clause</td>
	</tr>
	<tr>
		<td>Meeko</td>
		<td>v{meeko}</td>
		<td>
			<a href="https://github.com/forlilab/Meeko">
				https://github.com/forlilab/Meeko
			</a>
		</td>
		<td>LGPL-2.1</td>
	</tr>
	<tr>
		<td>Plip</td>
		<td>v{plip}</td>
		<td>
			<a href="https://github.com/pharmai/plip">
				https://github.com/pharmai/plip
			</a>
		</td>
		<td>GPL-2.0</td>
	</tr>
	<tr>
		<td>Psutil</td>
		<td>v{psutil}</td>
		<td>
			<a href="https://github.com/giampaolo/psutil">
				https://github.com/giampaolo/psutil
			</a>
		</td>
		<td>BSD 3-Clause</td>
	</tr>
	<tr>
		<td>NumPy</td>
		<td>v{numpy}</td>
		<td>
			<a href="https://numpy.org/">
				https://numpy.org
			</a>
		</td>
		<td>BSD 3-Clause</td>
	</tr>
	<tr>
		<td>OpenMM</td>
		<td>v{openmm}</td>
		<td>
			<a href="https://openmm.org/">
				https://openmm.org
			</a>
		</td>
		<td>MIT/LGPL</td>
	</tr>
	<tr>
		<td>PDBFixer</td>
		<td>v{pdbfixer}</td>
		<td>
			<a href="https://github.com/openmm/pdbfixer">
				https://github.com/openmm/pdbfixer
			</a>
		</td>
		<td>MIT</td>
	</tr>
	<tr>
		<td>PDB2PQR</td>
		<td>v{pdb2pqr}</td>
		<td>
			<a href="https://www.poissonboltzmann.org/">
				https://www.poissonboltzmann.org
			</a>
		</td>
		<td>BSD</td>
	</tr>
	<tr>
		<td>Requests</td>
		<td>v{requests}</td>
		<td>
			<a href="https://github.com/psf/requests">
				https://github.com/psf/requests
			</a>
		</td>
		<td>Apache-2.0</td>
	</tr>
	<tr>
		<td>Icons</td>
		<td>v{icon}</td>
		<td>
			<a href="https://icons.getbootstrap.com/">
				https://icons.getbootstrap.com
			</a>
		</td>
		<td>MIT</td>
	</tr>
	<tr>
		<td>AutoDock</td>
		<td>v4.2.6</td>
		<td>
			<a href="https://autodock.scripps.edu/">
				https://autodock.scripps.edu
			</a>
		</td>
		<td>GPL</td>
	</tr>
	<tr>
		<td>Vina</td>
		<td>v1.2.5</td>
		<td>
			<a href="https://github.com/ccsb-scripps/AutoDock-Vina/">
				https://github.com/ccsb-scripps/AutoDock-Vina
			</a>
		</td>
		<td>Apache-2.0</td>
	</tr>
	<tr>
		<td>QuickVina-W</td>
		<td>v1.1</td>
		<td>
			<a href="https://qvina.github.io/">
				https://qvina.github.io
			</a>
		</td>
		<td>Apache-2.0</td>
	</tr>
</table>
""".format(
	python = sys.version.split()[0],
	apsw = apsw.apswversion(),
	pyside = pyside_version,
	pymol = pymol.cmd.get_version()[0],
	babel = openbabel.__version__,
	rdkit = rdkit.__version__,
	meeko = meeko.__version__,
	plip = plip_version,
	psutil = psutil.__version__,
	openmm = openmm.__version__,
	pdbfixer = pdbfixer.pdbfixer.__version__,
	pdb2pqr = pdb2pqr.__version__,
	numpy = numpy.__version__,
	requests = requests.__version__,
	icon = '1.7.0'
)
