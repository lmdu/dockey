import psutil

meeko_file = psutil.__file__.replace('psutil', 'meeko')

with open(meeko_file) as fh:
	for line in fh:
		if line.startswith('__version__'):
			break
		else:
			line = ''

with open(meeko_file, 'w') as fw:
	fw.write(line)
