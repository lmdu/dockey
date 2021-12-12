import apsw

__all__ = ['DB']

DB_TABLES = {
	'molecular': [
		('id', 'INTEGER PRIMARY KEY'),
		('name', 'TEXT'),
		('type', 'INTEGER'),
		('format', 'TEXT'),
		('path', 'TEXT')
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
	'option': [
		('id', 'INTEGER PRIMARY KEY'),
		('name', 'TEXT'),
		('value', 'TEXT')
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
		self.connect()

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
		if self.conn:
			self.conn.close()

		self.conn = apsw.Connection(db_file)
		#self.conn.setrowtrace(row_factory)
		self._create_tables()

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
		table = sql.split()[3]
		fields = get_fields(table)
		row = cur.fetchone()
		if row:
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


DB = DataBackend()
