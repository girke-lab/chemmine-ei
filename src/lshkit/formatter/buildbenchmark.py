"""build benchmark data from eucsearch or chemsearch format"""
K = 100
import sys

def convert(fo, ofo=sys.stdout):
	fo.next()
	fo.next()
	for lid, line in enumerate(fo):
		for i, pair in enumerate(line.split()[:K]):
			idx, dist = pair.split(':')
			dist = float(dist)
			idx = int(idx) - 1
			if i == 0:
				assert dist == 0
				ofo.write("%d\t%s\t" % (idx, K))
			ofo.write("\t%d\t%s" % (idx, dist))
		ofo.write("\n")
		sys.stderr.write("%d     \r" % lid)
	sys.stderr.write("\n")
	
if __name__ == '__main__':
	if len(sys.argv) == 3:
		convert(file(sys.argv[1]), file(sys.argv[2], 'w'))
	elif len(sys.argv) == 2:
		convert(file(sys.argv[1]))
	elif len(sys.argv) == 1:
		sys.stderr.write("Reading from STDIN...\n")
		convert(sys.stdin)
	else:
		sys.stderr.write("Usage: %s input [output]\n" % sys.argv[0])
		sys.exit(1)

