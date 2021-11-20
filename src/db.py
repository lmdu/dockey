from PySide6.QtSql import *

__all__ = ['DockeyDatabase']

DB_TABLES = {
	'molecular': [
		('id', 'INTEGER PRIMARY KEY'),
		('name', 'TEXT'),
		('type', 'INTEGER'),
		('content', 'TEXT')
	],
	'tasks': [
		('id', 'INTEGER PRIMARY KEY'),
		('receptor', 'TEXT'),
		('ligand', 'TEXT'),
		('status', 'TEXT'),
		('created', 'INTEGER'),
		('finished', 'INTEGER')
	]
}

class DockeyDatabase:
	_conn = None

	def __init__(self):
		self.query = None

	def __del__(self):
		self._conn.commit()

	def __initialize_db(self):
		self.query = QSqlQuery()

		for table, fields in DB_TABLES.items():
			columns = ','.join(["{} {}".format(*field) for field in fields])
			sql = "CREATE TABLE IF NOT EXISTS {} ({})".format(table, columns)
			self.query.exec(sql)

	def connect(self, db_file):
		self._conn = QSqlDatabase.addDatabase("QSQLITE")
		self._conn.setDatabaseName(db_file)
		if self._conn.open():
			self.__initialize_db()
			return True
		else:
			return False

	def prepare(self, table):
		self.query = QSqlQuery()
		signs = ','.join(['?'] * len(DB_TABLES[table]))
		sql = "INSERT INTO {} VALUES ({})".format(table, signs)
		self.query.prepare(sql)

	def insert(self, *args):
		self.query.addBindValue(None)
		for arg in args:
			self.query.addBindValue(arg)
		self.query.exec()

		return self.query.lastInsertId()

	def get_content(self, name):
		query = QSqlQuery()
		query.prepare("SELECT content FROM molecular WHERE name=?")
		query.addBindValue(name)
		query.exec()

		while query.next():
			return query.value('content')

