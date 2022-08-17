import meeko

with open(meeko.__file__, 'w') as fw:
	fw.write('')

with open(meeko.__file__) as fh:
	print(fh.read())
