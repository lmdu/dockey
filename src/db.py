from PySide6.QtSql import *

from utils import AttrDict

__all__ = ['DB']

DB_TABLES = {
	'molecular': [
		('id', 'INTEGER PRIMARY KEY'),
		('name', 'TEXT'),
		('type', 'INTEGER'),
		('content', 'TEXT'),
		('format', 'TEXT')
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
		('status', 'TEXT')
	]
}

class DB:

	@staticmethod
	def connect(db_file):
		db = QSqlDatabase.addDatabase("QSQLITE")
		db.setDatabaseName(db_file)
		if db.open():
			query = QSqlQuery()

			for table, fields in DB_TABLES.items():
				columns = ','.join(["{} {}".format(*field) for field in fields])
				sql = "CREATE TABLE IF NOT EXISTS {} ({})".format(table, columns)
				query.exec(sql)

			return True
		else:
			return False

	@staticmethod
	def prepare(table):
		query = QSqlQuery()
		signs = ','.join(['?'] * len(DB_TABLES[table]))
		sql = "INSERT INTO {} VALUES ({})".format(table, signs)
		query.prepare(sql)
		return query

	@staticmethod
	def insert(query, *args):
		query.addBindValue(None)
		for arg in args:
			query.addBindValue(arg)
		query.exec()

		return query.lastInsertId()

	@staticmethod
	def get_receptors():
		query = QSqlQuery()
		query.exec("SELECT * FROM molecular WHERE type=1")
		while query.next():
			yield AttrDict({
				'id': query.value('id'),
				'name': query.value('name')
			})

	@staticmethod
	def get_ligands():
		query = QSqlQuery()
		query.exec("SELECT * FROM molecular WHERE type=2")
		while query.next():
			yield AttrDict({
				'id': query.value('id'),
				'name': query.value('name')
			})

	@staticmethod
	def get_molecular_by_id(_id):
		query = QSqlQuery()
		query.prepare("SELECT * FROM molecular WHERE id=? LIMIT 1")
		query.addBindValue(_id)
		query.exec()

		while query.next():
			return AttrDict({
				'id': query.value('id'),
				'name': query.value('name'),
				'type': query.value('type'),
				'content': query.value('content'),
				'format': query.value('format')
			})

	@staticmethod
	def get_molecular_by_name(name):
		query = QSqlQuery()
		query.prepare("SELECT * FROM molecular WHERE name=? LIMIT 1")
		query.addBindValue(name)
		query.exec()

		while query.next():
			return AttrDict({
				'id': query.value('id'),
				'name': name,
				'type': query.value('type'),
				'content': query.value('content'),
				'format': query.value('format')
			})

	@staticmethod
	def get_grid_box(_id):
		query = QSqlQuery()
		query.prepare("SELECT * FROM grid WHERE id=? LIMIT 1")
		query.addBindValue(_id)
		query.exec()

		while query.next():
			return AttrDict({
				'x': query.value('x'),
				'y': query.value('y'),
				'z': query.value('z'),
				'cx': query.value('cx'),
				'cy': query.value('cy'),
				'cz': query.value('cz'),
				'spacing': query.value('spacing')
			})
		else:
			return None

	@staticmethod
	def update_grid_box(*args):
		query = QSqlQuery()
		query.prepare("SELECT 1 FROM grid WHERE id=?")
		query.addBindValue(args[-1])

		if query.first():
			query.prepare(
				"UPDATE grid SET x=?, y=?, z=?, "
				"cx=?, cy=?, cz=?, spacing=? "
				"WHERE rid=?"
			)

			for v in args[1:]:
				query.addBindValue(v)

			query.addBindValue(args[0])
			query.exec()

		else:
			query.prepare("INSERT INTO grid VALUES (?,?,?,?,?,?,?,?,?)")
			query.addBindValue(None)

			for v in args:
				query.addBindValue(v)

			query.exec()

	@staticmethod
	def delete_grid_box(_id):
		query = QSqlQuery()
		query.prepare("DELETE FROM grid WHERE rid=?")
		query.addBindValue(_id)
		query.exec()
