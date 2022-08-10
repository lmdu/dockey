import sys
import cairosvg

for img in sys.argv[1:]:
	out = img.replace('svg', 'png')
	cairosvg.svg2png(url=img, write_to=out, output_width=256, output_height=256)
