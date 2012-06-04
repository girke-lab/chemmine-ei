#!/usr/bin/env python
"""
build the fingerprint library from PubChem SDF: read the fingerprint and 
write to a binary database
"""
HEADER = "DESCDB" + "\x00" * 10
import sys
from sdfiterator import sdf_iter
from base64 import decodestring
from struct import unpack, pack
from pubchemfp import getFpcalc

def init_db(fo):
	"""init the database. fo is a file object"""
	fo.write(HEADER)
	fo.write(pack('B', 4))

def close_db(fo):
	"""finishing writing to db"""
	fo.close()

def write_db(fo, data):
	fo.write(pack('i', len(data) / 4))	# every unsigned int has 4 bytes
	fo.write(data)

def read_fp(sdf):
	"""read the fingerprint from the sdf string of one molecule"""
	ready = False
	fp = None
	for i in sdf.splitlines():
		if ready:
			fp = i.strip()
			break
		if i.strip() == "> <PUBCHEM_CACTVS_SUBSKEYS>":
			ready = True
	
	if fp is None:
		fp = getFpcalc().calc(sdf)
	fp = decodestring(fp) + '\x00'
	assert unpack('>l112B', fp)[0] == 881
	return fp[4:]

def create_db(fp, db_fp, log_names=True, first=False, callback=None):
	"""create databse from a file path <fp>"""
	db = file(db_fp, 'w')
	if log_names: db_cids = file(db_fp + '.names', 'w')
	init_db(db)
	for i, sdf in enumerate(sdf_iter(fp)):
		if callable(callback):
			cb = callback(i)
			if cb == 'continue':
				continue
			if cb == 'break':
				break
			if cb == False:
				continue
			assert cb == True
		if log_names: db_cids.write(sdf.splitlines()[0] + '\n')
		fp = read_fp(sdf)
		write_db(db, fp)
		if first: break
	close_db(db)
	if log_names: db_cids.close()

if __name__ == '__main__':
	if len(sys.argv) != 3:
		sys.stderr.write("Usage: create_db.py mol.sdf db.cdb\n")
		sys.exit(1)

	create_db(sys.argv[1], sys.argv[2])
