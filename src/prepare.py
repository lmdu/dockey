#Modified from MGLTools v1.5.7 prepare_ligand4.py and prepare_receptor4.py
import os
import MolKit.molecule
import MolKit.protein

from rdkit import Chem
from MolKit import Read
from PySide6.QtCore import QSettings
from meeko.preparation import MoleculePreparation
#from meeko.utils.obutils import OBMolSupplier
from AutoDockTools.MoleculePreparation import AD4LigandPreparation
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation

from utils import load_molecule_from_file

__all__ = ['prepare_autodock_ligand', 'prepare_autodock_receptor',
	'prepare_meeko_ligand', 'prepare_ligand'
]

def prepare_autodock_ligand(ligand_file, ligand_pdbqt):
	settings = QSettings()
	settings.beginGroup('Ligand')
	repairs = settings.value('repairs', '')
	charges_to_add = settings.value('charges_to_add', 'gasteiger')
	#preserve_charge_types = settings.value('preserve_charge_types', '')
	cleanup = settings.value('cleanup', 'nphs_lps')
	allowed_bonds = settings.value('allowed_bonds', 'backbone')
	check_for_fragments = settings.value('check_for_fragments', False, bool)
	inactivate_all_torsions = settings.value('inactivate_all_torsions', False, bool)
	attach_nonbonded_fragments = settings.value('attach_nonbonded_fragments', False, bool)
	attach_singletons = settings.value('attach_singletons', False, bool)
	settings.endGroup()

	if charges_to_add == 'None':
		charges_to_add = None

	if attach_singletons:
		attach_nonbonded_fragments = True

	mols = Read(ligand_file)
	#if more than one moleucle, use the one molecule with the most atoms
	mol = max(mols, key=lambda x: len(x.allAtoms))
	coord_dict = {a: a.coords for a in mol.allAtoms}

	mol.buildBondsByDistance()

	'''
	if charges_to_add is not None:
		preserved = {}
		preserved_types = preserve_charge_types.split(',')

		for t in preserved_types:
			if not t:
				continue

			ats = mol.allAtoms.get(lambda x: x.autodock_element == t)

			for a in ats:
				if a.chargeSet is not None:
					preserved[a] = [a.chargeSet, a.charge]
	'''

	AD4LigandPreparation(mol,
		repairs = repairs,
		charges_to_add = charges_to_add,
		cleanup = cleanup,
		allowed_bonds = allowed_bonds,
		root = 'auto',
		outputfilename = ligand_pdbqt,
		check_for_fragments = check_for_fragments,
		bonds_to_inactivate = [],
		inactivate_all_torsions = inactivate_all_torsions,
		attach_nonbonded_fragments = attach_nonbonded_fragments,
		attach_singletons = attach_singletons
	)

	'''
	if charges_to_add is not None:
		for atom, chargeList in preserved.items():
			atom._charges[chargeList[0]] = chargeList[1]
			atom.chargeSet = chargeList[0]
	'''

	if mol.returnCode != 0:
		raise Exception(mol.returnMsg)

#modified from prepare_receptor4.py
def prepare_autodock_receptor(receptor_file, receptor_pdbqt):
	settings = QSettings()
	settings.beginGroup('Receptor')
	repairs = settings.value('repairs', '')
	charges_to_add = settings.value('charges_to_add', 'gasteiger')
	#preserve_charge_types = settings.value('preserve_charge_types', '')
	cleanup = settings.value('cleanup', 'nphs_lps_waters_nonstdres')
	delete_single_nonstd_residues = settings.value('delete_single_nonstd_residues', False, bool)
	unique_atom_names = settings.value('unique_atom_names', False, bool)
	settings.endGroup()

	if charges_to_add == 'None':
		charges_to_add = None

	mols = Read(receptor_file)
	#mol = max(mols, key=lambda x: len(x.allAtoms))
	mol = mols[0]

	if unique_atom_names:
		for at in mol.allAtoms:
			if mol.allAtoms.get(at.name) > 1:
				at.name = "{}{}".format(at.name, at._uniqIndex+1)

	'''
	preserved = {}
	has_autodock_element = False

	if charges_to_add and preserve_charge_types:
		if hasattr(mol, 'allAtoms') and not hasattr(mol.allAtoms[0], 'autodock_element'):
			_, ext = os.path.splitext(receptor_infile)

			if ext == '.pdbqt':
				has_autodock_element = True

		if preserve_charge_types and not has_autodock_element:
			raise Exception("receptor file does not have autodock_element, unable to preserve charges on {}".format(preserve_charge_types))

		preserved_types = preserve_charge_types.split(',')

		for t in preserved_types:
			if not t:
				continue

			ats = mol.allAtoms.get(lambda x: x.autodock_element==t)

			for a in :
				if a.chargeSet is not None:
					preserved[a] = [a.chargeSet, a.charge]
	'''

	if len(mols) > 1:
		mol = max(mols, key=lambda x: len(x.allAtoms))

	mol.buildBondsByDistance()

	AD4ReceptorPreparation(mol,
		repairs = repairs,
		charges_to_add=charges_to_add,
		cleanup = cleanup, 
		outputfilename = receptor_pdbqt,
		#preserved = preserved,
		delete_single_nonstd_residues=delete_single_nonstd_residues
	)

	'''
	if charges_to_add:
		for atom, chargeList in preserved.items():
			atom._charges[chargeList[0]] = chargeList[1]
			atom.chargeSet = chargeList[0]
	'''

	if mol.returnCode != 0:
		raise Exception(mol.returnMsg)

def prepare_meeko_ligand(ligand_file, ligand_pdbqt):
	settings = QSettings()
	settings.beginGroup('Meeko')
	rigid_macrocycles = settings.value('rigid_macrocycles', False, bool)
	keep_chorded_rings = settings.value('keep_chorded_rings', False, bool)
	keep_equivalent_rings = settings.value('keep_equivalent_rings', False, bool)
	hydrate = settings.value('hydrate', False, bool)
	keep_nonpolar_hydrogens = settings.value('keep_nonpolar_hydrogens', False, bool)
	flexible_amides = settings.value('flexible_amides', False, bool)
	add_index_map = settings.value('add_index_map', False, bool)
	remove_smiles = settings.value('remove_smiles', False, bool)
	rigidify_bonds_smarts = settings.value('rigidify_bonds_smarts', '')
	rigidify_bonds_indices = settings.value('rigidify_bonds_indices', '')
	atom_type_smarts = settings.value('atom_type_smarts', '')
	double_bond_penalty = settings.value('double_bond_penalty', 50, int)
	settings.endGroup()

	if rigidify_bonds_smarts.strip():
		rigidify_bonds_smarts = rigidify_bonds_smarts.strip().split(',')
	else:
		rigidify_bonds_smarts = []

	rbi_list = []
	if rigidify_bonds_indices.strip():
		for item in rigidify_bonds_indices.strip().split(','):
			for i, j in item.split():
				rbi_list.append([int(i)-1, int(j)-1])
	rigidify_bonds_indices = rbi_list

	if atom_type_smarts.strip():
		atom_type_smarts = json.load(atom_type_smarts.strip())
	else:
		atom_type_smarts = {}

	#mols = OBMolSupplier(ligand_file, 'pdb')
	#mols = iter(mols)
	#mol = next(mols)
	#mol = load_molecule_from_file(ligand_file, 'pdb')
	mol = Chem.MolFromPDBFile(ligand_file)

	preparator = MoleculePreparation(
		keep_nonpolar_hydrogens = keep_nonpolar_hydrogens,
		hydrate = hydrate,
		flexible_amides = flexible_amides,
		rigid_macrocycles = rigid_macrocycles,
		keep_chorded_rings = keep_chorded_rings,
		keep_equivalent_rings = keep_equivalent_rings,
		rigidify_bonds_smarts = rigidify_bonds_smarts,
		rigidify_bonds_indices = rigidify_bonds_indices,
		double_bond_penalty = double_bond_penalty,
		atom_type_smarts = atom_type_smarts,
		add_index_map = add_index_map,
		remove_smiles = remove_smiles
	)

	preparator.prepare(mol)
	preparator.write_pdbqt_file(ligand_pdbqt)

def prepare_ligand(ligand_file, ligand_pdbqt):
	settings = QSettings()
	tool = settings.value('Ligand/prepare_tool', 'prepare_ligand4')

	if tool == 'prepare_ligand4':
		prepare_autodock_ligand(ligand_file, ligand_pdbqt)
	elif tool == 'meeko':
		prepare_meeko_ligand(ligand_file, ligand_pdbqt)
