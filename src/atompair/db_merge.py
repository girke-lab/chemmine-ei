HEADER = "DESCDB" + "\x00" * 10
import sys
from struct import pack, unpack

def init_db(fo):
	"""init the database. fo is a file object"""
	fo.write(HEADER)
	fo.write(pack('B', 4))

def copy_db(fo, db):
	"""copy content of db at 'db' into fo"""
	f = file(db)
	x = f.read(17)
	assert unpack('B', x[-1])[0] == 4
	while True:
		x = f.read(1024)
		if not x: break
		fo.write(x)
	f.close()

if __name__ == '__main__':
	if len(sys.argv) < 2:
		sys.stderr.write("db_merge out.cdb in.cdb ...\n")
		sys.exit(1)

	f = file(sys.argv[1], 'w')
	init_db(f)
	for i in sys.argv[2:]:
		copy_db(f, i)
	f.close()
