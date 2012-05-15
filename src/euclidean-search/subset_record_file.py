from struct import pack, unpack

def pack_header(n, dim):
	return pack('III', 4, n, dim)

if __name__ == '__main__':
	import sys
	try:
		inp = sys.argv[1]
		indp = sys.argv[2]
		outp = sys.argv[3]
		f = file(inp)
		indf = file(indp)
		of = file(outp, 'w')
	except:
		sys.stderr.write("Usage: subset_record_file.py matrix index.file output\n")
		sys.exit(1)
	
	x = f.read(4 * 3);
	_, n, dim = unpack('III', x)
	assert _ == 4
	sys.stderr.write("%d x %d\n" % (dim, n))

	indices = [int(i) for i in indf]
	m = len(indices)
	sys.stderr.write("subsetting %d lines\n" % m)
	
	of.write(pack_header(m, dim))
	for i in range(n):
		x = f.read(4 * dim)
		if (i + 1) in indices:
			of.write(x)
			m -= 1
	assert m == 0
	of.close()
	f.close()
	indf.close()
