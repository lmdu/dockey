#Modified from MGLTools v1.5.7 prepare_ligand4.py and prepare_receptor4.py
import os
import MolKit.molecule
import MolKit.protein

from rdkit import Chem
from MolKit import Read
from MolKit.molecule import BondSet
from MolKit.protein import ProteinSet, ResidueSet, AtomSet
from meeko.preparation import MoleculePreparation
from AutoDockTools.MoleculePreparation import AD4LigandPreparation
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation
from AutoDockTools.MoleculePreparation import AD4FlexibleReceptorPreparation

from utils import load_molecule_from_file

__all__ = ['prepare_autodock_ligand', 'prepare_autodock_receptor',
	'prepare_meeko_ligand', 'prepare_ligand'
]

def prepare_autodock_ligand(ligand_file, ligand_pdbqt, params):
	repairs = params['repairs']
	charges_to_add = params['charges_to_add']
	cleanup = params['cleanup']
	allowed_bonds = params['allowed_bonds']
	check_for_fragments = params['check_for_fragments']
	inactivate_all_torsions = params['inactivate_all_torsions']
	attach_nonbonded_fragments = params['attach_nonbonded_fragments']
	attach_singletons = params['attach_singletons']

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
		#mode = 'automatic',
		repairs = repairs,
		charges_to_add = charges_to_add,
		cleanup = cleanup,
		allowed_bonds = allowed_bonds,
		#root = 'auto',
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
def prepare_autodock_receptor(receptor_file, receptor_pdbqt, params):
	repairs = params['repairs']
	charges_to_add = params['charges_to_add']
	cleanup = params['cleanup']
	delete_single_nonstd_residues = params['delete_single_nonstd_residues']
	unique_atom_names = params['unique_atom_names']

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

#modified from prepare_flexreceptor4.py
def prepare_flex_receptor(receptor_pdbqt, res_dict):
	#mol.chains.get(lambda x: x.id == chainID)
	rigid_file = receptor_pdbqt.replace('.pdbqt', '_rigid.pdbqt')
	flex_file = receptor_pdbqt.replace('.pdbqt', '_flex.pdbqt')

	rec = Read(receptor_pdbqt)[0]
	bnds = rec.buildBondsByDistance()
	all_res = ResidueSet()

	for chain_id, res_names in res_dict.items():
		chains = rec.chains.get(lambda x: x.id == chain_id)
		res = ResidueSet()
		for res_name in res_names:
			matched = chains.residues.get(lambda x: x.name == res_name)

			if matched.__class__ == AtomSet:
				res = matched.parent.uniq()

			if len(res):
				all_res += res

	#check for duplicates
	d = {res: 1 for res in all_res}
	all_res = list(d.keys())
	all_res = ResidueSet(all_res).uniq()
	all_res.sort()

	all_bnds = BondSet()

	AD4FlexibleReceptorPreparation(rec,
		residues = all_res,
		rigid_filename = rigid_file, 
		flexres_filename = flex_file,
		non_rotatable_bonds = all_bnds,
		mode = 'automatic'
	)

def prepare_meeko_ligand(ligand_file, ligand_pdbqt, params):
	rigid_macrocycles = params['rigid_macrocycles']
	keep_chorded_rings = params['keep_chorded_rings']
	keep_equivalent_rings = params['keep_equivalent_rings']
	hydrate = params['hydrate']
	keep_nonpolar_hydrogens = params['keep_nonpolar_hydrogens']
	flexible_amides = params['flexible_amides']
	add_index_map = params['add_index_map']
	remove_smiles = params['remove_smiles']
	rigidify_bonds_smarts = params['rigidify_bonds_smarts']
	rigidify_bonds_indices = params['rigidify_bonds_indices']
	atom_type_smarts = params['atom_type_smarts']
	double_bond_penalty = params['double_bond_penalty']

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

	ligand_format = os.path.splitext(ligand_file)[1]

	if ligand_format == '.pdb':
		mol = Chem.MolFromPDBFile(ligand_file, removeHs=False)

	elif ligand_format == '.mol2':
		mol = Chem.MolFromMol2File(ligand_file, removeHs=False)

	elif ligand_format == '.mol':
		mol = Chem.MolFromMolFile(ligand_file, removeHs=False)

	elif ligand_format == '.sdf':
		#only get the first mol in sdf file
		for mol in Chem.SDMolSupplier(ligand_file, removeHs=False):
			break

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

def prepare_ligand(ligand_file, ligand_pdbqt, params):
	tool = params['tool']

	if tool == 'prepare_ligand4':
		prepare_autodock_ligand(ligand_file, ligand_pdbqt, params)
	elif tool == 'meeko':
		prepare_meeko_ligand(ligand_file, ligand_pdbqt, params)

if __name__ == '__main__':
	prepare_flex_receptor('las2.pdbqt', {'A': ['LYS16', 'GLY126']})
