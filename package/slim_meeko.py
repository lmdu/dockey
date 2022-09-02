import psutil

meeko_file = psutil.__file__.replace('psutil', 'meeko')
with open(meeko_file, 'w') as fw:
	fw.write('')
