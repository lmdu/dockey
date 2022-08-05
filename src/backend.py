import apsw

__all__ = ['DB']

DB_TABLES = {
	'molecular': [
		('id', 'INTEGER PRIMARY KEY'),
		('name', 'TEXT'),
		('type', 'INTEGER'),
		('pdb', 'TEXT'),
		('atoms', 'INTEGER'),
		('bonds', 'INTEGER'),
		('hvyatoms', 'INTEGER'),
		('residues', 'INTEGER'),
		('rotors', 'INTEGER'),
		('formula', 'TEXT'),
		('energy', 'REAL'),
		('weight', 'REAL'),
		('logp', 'REAL')
	],
	'grid': [
		('id', 'INTEGER PRIMARY KEY'),
		('rid', 'INTEGER'),
		('x', 'INTEGER'),
		('y', 'INTEGER'),
		('z', 'INTEGER'),
		('cx', 'REAL'),
		('cy', 'REAL'),
		('cz', 'REAL'),
		('spacing', 'REAL')
	],
	'jobs': [
		('id', 'INTEGER PRIMARY KEY'),
		('rid', 'INTEGER'),
		('lid', 'INTEGER'),
		('status', 'INTEGER'),
		('progress', 'REAL'),
		('started', 'INTEGER'),
		('finished', 'INTEGER'),
		('message', 'TEXT')
	],
	'pose': [
		('id', 'INTEGER PRIMARY KEY'),
		('jid', 'INTEGER'),
		('run', 'INTEGER'),
		('energy', 'REAL'),
		('rmsd1', 'REAL'),
		('rmsd2', 'REAL'),
		('logki', 'REAL'),
		('le', 'REAL'),
		('sile', 'REAL'),
		('fq', 'REAL'),
		('lle', 'REAL'),
		('lelp', 'REAL'),
		('ki', 'TEXT'),
		('mode', 'TEXT'),
		('complex', 'TEXT')
	],
	'best': [
		('id', 'INTEGER PRIMARY KEY'),
		('pid', 'INTEGER')
	],
	'logs': [
		('id', 'INTEGER PRIMARY KEY'),
		('jid', 'INTEGER'),
		('name', 'TEXT'),
		('content', 'TEXT')
	],
	'option': [
		('id', 'INTEGER PRIMARY KEY'),
		('name', 'TEXT'),
		('value', 'TEXT')
	],
	'binding_site': [
		('id', 'INTEGER PRIMARY KEY'),
		('pid', 'INTEGER'),
		('site', 'TEXT')
	],
	'hydrogen_bond': [
		('id', 'INTEGER PRIMARY KEY'),
		('bid', 'INTEGER'),
		('chain', 'TEXT'),
		('residue', 'INTEGER'),
		('amino_acid', 'TEXT'),
		('distance_ha', 'REAL'),
		('distance_da', 'REAL'),
		('donor_angle', 'REAL'),
		('protein_donor', 'TEXT'),
		('side_chain', 'TEXT'),
		('donor_atom', 'TEXT'),
		('acceptor_atom', 'TEXT')
	],
	'halogen_bond': [
		('id', 'INTEGER PRIMARY KEY'),
		('bid', 'INTEGER'),
		('chain', 'TEXT'),
		('residue', 'INTEGER'),
		('amino_acid', 'TEXT'),
		('distance', 'REAL'),
		('donor_angle', 'REAL'),
		('acceptor_angle', 'REAL'),
		('donor_atom', 'TEXT'),
		('acceptor_atom', 'TEXT')
	],
	'hydrophobic_interaction': [
		('id', 'INTEGER PRIMARY KEY'),
		('bid', 'INTEGER'),
		('chain', 'TEXT'),
		('residue', 'INTEGER'),
		('amino_acid', 'TEXT'),
		('distance', 'REAL'),
		('ligand_atom', 'INTEGER'),
		('protein_atom', 'INTEGER')
	],
	'salt_bridge': [
		('id', 'INTEGER PRIMARY KEY'),
		('bid', 'INTEGER'),
		('chain', 'TEXT'),
		('residue', 'INTEGER'),
		('amino_acid', 'TEXT'),
		('distance', 'REAL'),
		('protein_positive', 'INTEGER'),
		('ligand_group', 'TEXT'),
		('ligand_atoms', 'TEXT')
	],
	'water_bridge': [
		('id', 'INTEGER PRIMARY KEY'),
		('bid', 'INTEGER'),
		('chain', 'TEXT'),
		('residue', 'INTEGER'),
		('amino_acid', 'TEXT'),
		('distance_aw', 'REAL'),
		('distance_dw', 'REAL'),
		('donor_angle', 'REAL'),
		('water_angle', 'REAL'),
		('protein_donor', 'INTEGER'),
		('donor_atom', 'REAL'),
		('acceptor_atom', 'REAL'),
		('water_atom', 'INTEGER')
	],
	'pi_stacking': [
		('id', 'INTEGER PRIMARY KEY'),
		('bid', 'INTEGER'),
		('chain', 'TEXT'),
		('residue', 'INTEGER'),
		('amino_acid', 'TEXT'),
		('distance', 'REAL'),
		('angle', 'REAL'),
		('offset', 'REAL'),
		('type', 'TEXT'),
		('ligand_atoms', 'TEXT')
	],
	'pi_cation': [
		('id', 'INTEGER PRIMARY KEY'),
		('bid', 'INTEGER'),
		('chain', 'TEXT'),
		('residue', 'INTEGER'),
		('amino_acid', 'TEXT'),
		('distance', 'REAL'),
		('offset', 'REAL'),
		('protein_charged', 'INTEGER'),
		('ligand_group', 'TEXT'),
		('ligand_atoms', 'TEXT')
	],
	'metal_complex': [
		('id', 'INTEGER PRIMARY KEY'),
		('bid', 'INTEGER'),
		('chain', 'TEXT'),
		('residue', 'INTEGER'),
		('amino_acid', 'TEXT'),
		('metal', 'INTEGER'),
		('target', 'INTEGER'),
		('distance', 'REAL'),
		('location', 'TEXT')
	]
}

def get_fields(table):
	return [field[0] for field in DB_TABLES[table]]

class DataRow(dict):
	def __getattr__(self, attr):
		return self[attr]

def row_factory(cursor, row):
	fields = [name for name, _ in cursor.getdescription()]
	return DataRow(zip(fields, row))

class DataBackend:
	conn = None

	def __init__(self):
		#self.connect()
		pass

	def __del__(self):
		if self.conn:
			self.conn.close()

	@property
	def cursor(self):
		return self.conn.cursor()

	def query(self, sql, paras=None):
		if paras is None:
			return self.cursor.execute(sql)
		else:
			return self.cursor.execute(sql, paras)

	def _create_tables(self):
		for table, fields in DB_TABLES.items():
			columns = ','.join(["{} {}".format(*field) for field in fields])
			sql = "CREATE TABLE IF NOT EXISTS {} ({})".format(table, columns)
			self.query(sql)

	def connect(self, db_file=':memory:'):
		self.close()
		self.conn = apsw.Connection(db_file)
		#self.conn.setrowtrace(row_factory)
		self._create_tables()

	def close(self):
		if self.conn:
			self.conn.close()
			self.conn = None

	def active(self):
		return self.conn is not None

	def save(self, db_file):
		target = apsw.Connection(db_file)
		with target.backup('main', self.conn, 'main') as b:
			b.step()

	def clear_table(self, table):
		self.query("DELETE FROM {}".format(table))

	def insert_rows(self, sql, rows):
		self.cursor.executemany(sql, rows)

	def table_exists(self, table):
		try:
			self.query("SELECT 1 FROM {}".format(table))
			return True
		except:
			return False

	def get_exists(self, sql, paras=None):
		for _ in self.query(sql, paras):
			return True
		return False

	def get_one(self, sql, paras=None):
		for row in self.query(sql, paras):
			return row[0]

	def get_row(self, sql, paras=None):
		for row in self.query(sql, paras):
			return row

	def get_dict(self, sql, paras=None):
		cur = self.query(sql, paras)
		for row in cur:
			fields = [col[0] for col in cur.getdescription()]
			return DataRow(zip(fields, row))

	def get_column(self, sql, paras=None):
		return [row[0] for row in self.query(sql, paras)]

	def get_set(self, sql, paras=None):
		return {row[0] for row in self.query(sql, paras)}

	def get_count(self, table):
		if self.table_exists(table):
			return self.get_one("SELECT COUNT(1) FROM {}".format(table))
		else:
			return 0

	def get_field(self, table):
		return [row[1] for row in self.query("PRAGMA table_info({})".format(table))]

	def get_field_type(self, table):
		return [row[2] for row in self.query("PRAGMA table_info({})".format(table))]

	def get_tables(self):
		sql = "SELECT name FROM sqlite_master WHERE type=?"
		return self.get_column(sql, ('table',))

	def get_option(self, name):
		sql = "SELECT value FROM option WHERE name=? LIMIT 1"
		return self.get_one(sql, (name,))

	def set_option(self, name, val):
		if self.get_option(name) is not None:
			sql = "UPDATE option SET value=? WHERE name=?"
			self.query(sql, (val, name))
		else:
			sql = "INSERT INTO option VALUES (?,?,?)"
			self.query(sql, (None, name, val))

DB = DataBackend()
