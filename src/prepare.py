#Modified from MGLTools v1.5.7 prepare_ligand4.py and prepare_receptor4.py
import MolKit.molecule
import MolKit.protein

from rdkit import Chem
from MolKit import Read
from PySide6.QtCore import QSettings
from meeko import MoleculePreparation
from AutoDockTools.MoleculePreparation import AD4LigandPreparation
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation

__all__ = ['prepare_autodock_ligand', 'prepare_autodock_receptor']

def prepare_autodock_ligand(ligand_file, ligand_pdbqt, repaires='bonds_hydrogens', add_charges='gasteiger',
	preserve_charges=[], cleanup_type='nphs_lps', rotate_bonds='backbone', root_index='auto',
	nonbonded_fragments=False, mode='automatic', inactivate_bonds='', inactivate_torsions=False,
	attach_fragments=False, attach_singletons=False):
	
	add_bonds = False

	if attach_singletons:
		attach_fragments = True

	mols = Read(ligand_file)
	#if more than one moleucle, use the one molecule with the most atoms
	mol = max(mols, key=lambda x: len(x.allAtoms))
	#coord_dict = {a: a.coords for a in mol.allAtoms}
	mol.buildBondsByDistance()

	AD4LigandPreparation(mol, mode, repaires, add_charges, cleanup_type, rotate_bonds,
		root_index, ligand_pdbqt, check_for_fragments=nonbonded_fragments,
		bonds_to_inactivate=inactivate_bonds, inactivate_all_torsions=inactivate_torsions,
		attach_nonbonded_fragments=attach_fragments, attach_singletons=attach_singletons
	)

	'''
	if add_charges:
		atoms = mol.allAtoms.get(lambda x: x.autodock_element in preserve_charges)
		preserved = {a: [a.chargeSet, a.charge] for a in atoms if a.chargeSet}

		for atom in preserved:
			charge_set, charge = preserved[atom]
			atom._charges[charge_set] = charge
	'''

	if mol.returnCode != 0:
		raise Exception(mol.returnMsg)

def prepare_autodock_receptor(receptor_infile, receptor_outfile, repaires='bonds_hydrogens',
	add_charges='gasteiger', preserve_charges=[], cleanup_type='nphs_lps_waters_nonstdres',
	delete_residues=None, mode='automatic', dictionary=None, unique_name=False):
	mols = Read(receptor_infile)
	mol = max(mols, key=lambda x: len(x.allAtoms))

	if unique_name:
		for atom in mol.allAtoms:
			if mol.allAtoms.get(atom.name) > 1:
				atom.name = "{}{}".format(atom.name, atom._uniqIndex+1)

	preserved = {}
	has_autodock_element = False

	if add_charges and preserve_charges:
		if hasattr(mol, 'allAtoms') and not hasattr(mol.allAtoms[0], 'autodock_element'):
			_, ext = os.path.splitext(receptor_infile)

			if ext == '.pdbqt':
				has_autodock_element = True

		if preserve_charges and not has_autodock_element:
			raise Exception("Unable to preserve charges")

		atoms = mol.allAtoms.get(lambda x: x.autodock_element in preserve_charges)
		preserved = {a: [a.chargeSet, a.charge] for a in atoms if a.chargeSet}

	mol.buildBondsByDistance()

	AD4ReceptorPreparation(mol, mode, repaires, add_charges, cleanup_type, receptor_outfile,
		preserved = preserved, delete_single_nonstd_residues=delete_residues, dict_file=dictionary
	)

	if add_charges:
		for atom in preserved:
			charge_set, charge = preserved[atom]
			atom._charges[charge_set] = charge

	if mol.returnCode != 0:
		raise Exception(mol.returnMsg)
