#!/usr/bin/env python
"""
build the fingerprint library from PubChem SDF in parallel: read the fingerprint
and write to a binary database
"""
import sys
from fpcdb import create_db

def id_ok(per_job, this_id, i):
	if i < per_job * this_id:
		return 'continue'
	if i >= per_job * (this_id + 1):
		return 'break'
	return True

if __name__ == '__main__':
	if len(sys.argv) != 5:
		sys.stderr.write("Usage: create_db.py mol.sdf db.cdb per_job job_id\n")
		sys.exit(1)
	
	per_job, job_id = int(sys.argv[3]), int(sys.argv[4])
	callback = lambda x: id_ok(per_job, job_id, x)
	out_fp = sys.argv[2].rsplit('.', 1)
	out_fp.insert(1, '%d' % job_id)
	out_fp = '.'.join(out_fp)

	create_db(sys.argv[1], out_fp, callback=callback)
