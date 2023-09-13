Import molecules
================

After creating or opening project file, you will be allowed to import receptors and ligands. The Dockey supports various formats that can be read by `OpenBabel <http://openbabel.org/docs/current/FileFormats/Overview.html>`_.

.. note::

	The Dockey will identify molecule file format through extension name. Please make sure that the molecule file extension is identical with its content.

Import Receptors
----------------

Import from local file
^^^^^^^^^^^^^^^^^^^^^^

Go **File** menu -> **Import Receptors** to select receptor files and click **Open** to import into Dockey.

Import from PDB database
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Go **File** menu -> **Import Receptor from Database** -> **PDB** to open a dialog:

.. rst-class:: wy-text-center

	|pdb|

Then input PDB IDs using comma to separate multiple ones and click ``OK`` button. The Dockey will automatically download the molecule from `RCSB PDB <https://www.rcsb.org/>`_ database and import it.

Import Ligands
--------------

Import from local file
^^^^^^^^^^^^^^^^^^^^^^

Go **File** menu -> **Import Ligands** to select ligand files and import to Dockey.

Import from SDF file
^^^^^^^^^^^^^^^^^^^^

If you want to import all ligands from SDF file database. You can go to **File** menu -> **Import Ligand from Database** -> **SDF** to select the sdf file with multiple molecules to import.

Import from Zinc database
^^^^^^^^^^^^^^^^^^^^^^^^^

Go **File** menu -> **Import Ligand from Database** -> **Zinc** to open a dialog:

.. rst-class:: wy-text-center

	|zinc|

Then input Zinc IDs using comma to separate multiple ones and click ``OK`` button. The Dockey will automatically download the molecule from `Zinc <https://zinc.docking.org/>`_ database and import it.

Import from PubChem database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Go **File** menu -> **Import Ligand from Database** -> **PubChem** to open a dialog:

.. rst-class:: wy-text-center

	|pubchem|

Then input PubChem IDs using comma to separate multiple ones and click ``OK`` button. The Dockey will automatically download the molecule from `PubChem <https://pubchem.ncbi.nlm.nih.gov/>`_ database and import it.

Import from ChEMBL database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Go **File** menu -> **Import Ligand from Database** -> **ChEMBL** to open a dialog:

.. rst-class:: wy-text-center

	|chembl|

Then input ChEMBL IDs using comma to separate multiple ones and click ``OK`` button. The Dockey will automatically download the molecule from `ChEMBL <https://www.ebi.ac.uk/chembl/>`_ database and import it.

Molecular List
--------------

The imported ligands and receptors will be separately displayed in molecular list.

.. rst-class:: wy-text-center

	|mol|

The molecular list has right-click menu. You can use **Import Receptors** and **Import Ligands** in menu list to import receptor and ligand files.

.. rst-class:: wy-text-center

	|molmenu|

You also allowed to use **Delete** to remove current selected molecule and **Delete All** to remove all molecules from Dockey.

You can use **View Details** to obtain detailed information of molecule including number of atoms, bonds, heavy atoms, residues and rotors, formula, molecular weight as well as calculated *logp*.

.. rst-class:: wy-text-center

	|molinfo|

.. |pdb| image:: _static/pdb.png
	:width: 400
.. |zinc| image:: _static/zinc.png
	:width: 400
.. |pubchem| image:: _static/pubchem.png
	:width: 400
.. |chembl| image:: _static/chembl.png
	:width: 400
.. |mol| image:: _static/molecules.png
	:width: 400
.. |molmenu| image:: _static/molmenus.png
	:width: 400
.. |molinfo| image:: _static/molinfo.png
	:width: 400
