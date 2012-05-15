from struct import pack

def guess_dim(f):
	i = 0
	dim = 0
	for line in f:
		if not dim: dim = len(line.split())
		i += 1
	sys.stderr.write("%d x %d\n" % (dim, i))
	return (dim, i)

def pack_header(n, dim):
	return pack('III', 4, n, dim)

def pack_line(line, dim=None):
	rst = ''
	if dim: assert dim == len(line.split())
	for i in line.split():
		rst += pack('f', float(i))
	return rst

if __name__ == '__main__':
	import sys
	try:
		inp = sys.argv[1]
		outp = sys.argv[2]
		f = file(inp)
		of = file(outp, 'w')
	except:
		sys.stderr.write("Usage: mk_record_file.py input.txt output\n")
		sys.exit(1)
	
	dim, n = guess_dim(f)
	f.seek(0)
	
	of.write(pack_header(n, dim))
	for line in f:
		of.write(pack_line(line, dim))
		n -= 1
	assert n == 0
	of.close()
	f.close()
