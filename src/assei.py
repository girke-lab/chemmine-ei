#!/usr/bin/env python
"""
the ASSEI demo module
"""
import os
import sys
sys.path.append('eutils')
sys.path.append('eutils/fp')
#sys.path.append('/home/ycao/dev/fp')
#sys.path.append('/home/ycao/dev/optim.2')
#sys.path.append('/home/ycao/dev/lsh-2/build/bin')
#sys.path.append('/home/ycao/dev/descriptor')
from fpcdb import create_db
from coord import CoordinateSolver
from fpdbcompare import DBComparer
#_dir_ = os.path.dirname(os.path.abspath(__file__))
#solver = CoordinateSolver(os.path.join(_dir_, 'puzzle'))
solver = CoordinateSolver('puzzle')
import socket
import sys
from base64 import b64encode
from Queue import Queue
queue = Queue()
from time import time
from stat import ST_SIZE
from traceback import print_exc
from lshsearch import LSHSearcher

import signal

# CONFIG SECTION
DISTMAT_BIN = '/home/ycao/dev/descriptor/fp_compare db-db-distance'
#REF_CDB = os.path.join(_dir_, 'ref.cdb')
REF_CDB = 'ref.cdb'
LSH_SERVER = ('himem', 50009)
REFINE_SERVER = ('mem3g', 50010)
COUNTER_FILE = 'counter'
# CONFIG ENDS

sys.stderr.write("loading names...\n")
#names_fp = os.path.join(_dir_, 'names')
names_fp =  'names'
names_f = file(names_fp)
names = [i.strip() for i in names_f]
names_f.close()
sys.stderr.write("done\n")

comparer = DBComparer(REF_CDB)


def lsh(query):
	""" lshsearch local """
	lshsearcher = LSHSearcher("matrix")
	try:
		start = time()
		result = lshsearcher.search(query)
		elapsed = time() - start
	except:
		print_exc()
		elapsed = 0
		result = '\n'
	return "/t:"+str(elapsed) + " " + result
	

def lshRemote(query):
	""" lshsearch TCP client """
	HOST, PORT = LSH_SERVER
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	sock.connect((HOST, PORT))
	hello = sock.recv(3)
	assert hello == 'OK\n'
	if not query.endswith('\n'):
		query += '\n'
	sock.send(query)
	data = ''

	while True:
		ret = sock.recv(1024*16)
		if not ret: break
		data += ret
	sock.close()
	return data

def refine(query_cdb, candidates):
	""" refine through TCP service """
	f = file(query_cdb)
	c = f.read()
	f.close()
	HOST, PORT = REFINE_SERVER
	sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	sock.connect((HOST, PORT))
	hello = sock.recv(3)
	assert hello == 'OK\n'
	query = b64encode(c) + ' ' + candidates
	if not query.endswith('\n'):
		query += '\n'
	sock.send(query)
	data = ''

	while True:
		ret = sock.recv(1024*16)
		if not ret: break
		data += ret
	sock.close()
	return data

def update_status(dir, status):
	status_fp = os.path.join(dir, 'status')
	status_f = file(status_fp, 'w')
	status_f.write(status)
	status_f.close()

def make(dir):
	""" process in <dir>. A file named 'query.sdf' is expected there """
	update_status(dir, 'start')

	t = time()
	assert os.path.isfile(os.path.join(dir, 'query.sdf'))
	inp = os.path.join(dir, 'query.sdf')
	outp = os.path.join(dir, 'query.cdb')
	cdb = outp
	error_fp = os.path.join(dir, 'error')
	try:
		create_db(inp, outp, log_names=False, first=True)
		assert os.stat(outp)[ST_SIZE] != 17
	except:
		print_exc()
		f = file(error_fp, 'w')
		f.write("Error in parsing SDF\n")
		f.close()
		update_status(dir, 'failed')
		return False
	parsing_time = time() - t

	# create distance matrix
	inp = outp
	outp = os.path.join(dir, 'query.dist')
	try: os.unlink(outp)
	except: pass
	t = time()
	failed = False
	try:
		comparer.compare(inp, outp)
	except:
		failed = True
	if failed or os.stat(outp)[ST_SIZE] == 0:
		f = file(error_fp, 'w')
		f.write("Error in comparing to references\n")
		f.close()
		update_status(dir, 'failed')
		return False

	# solve the coordinate
	inp = outp
	outp = os.path.join(dir, 'query.coord')
	try:
		ret = solver.solve(inp).strip()
		assert ret
	except:
		print_exc()
		f = file(error_fp, 'w')
		f.write("Error in embedding\n")
		f.close()
		update_status(dir, 'failed')
		return False
	f = file(outp, 'w')
	f.write(ret)
	f.close()
	coord_time = time() - t

	# perform lsh search
	t = time()
	try:
		ret = lsh(ret).strip()
		assert ret
		print ret
		if ret.startswith('/t:'):
			_, ret = ret.split(None, 1)
			lsh_time_ = float(_.split(':')[1])
		else:
			lsh_time_ = 0
	except:
		print_exc()
		f = file(error_fp, 'w')
		f.write("Error in performing LSH query\n")
		f.close()
		update_status(dir, 'failed')
		return False
	lsh_time = time() - t

	# refine
	t = time()
	try:
		ret = refine(cdb, ret)
		assert ret
		if ret.startswith('/t:'):
			_, ret = ret.split(None, 1)
			refine_time_ = float(_.split(':')[1])
		else:
			refine_time_ = 0
	except:
		print_exc()
		f = file(error_fp, 'w')
		f.write("Error in performing refinement\n")
		f.close()
		update_status(dir, 'failed')
		return False
	refine_time = time() - t

	communication_time = lsh_time - lsh_time_ + refine_time - refine_time_

	sys.stderr.write('timing: parsing=%s embedding=%s lsh=%s refine=%s communication=%s\n' %
		(parsing_time, coord_time, lsh_time_, refine_time_, communication_time))
	f = file(os.path.join(dir, 'out'), 'w')
	f.write('# %s %s %s %s %s\n' % 
		(parsing_time, coord_time, lsh_time_, refine_time_, communication_time))
	for pair in ret.split():
		seq_id, dist = pair.split(':')
		cid = names[int(seq_id) - 1]
		f.write('%s %s\n' % (cid, dist))
	f.close()

	if os.path.exists(error_fp):
		os.unlink(error_fp)
	
	update_status(dir, 'done')
	return True

def start_server():
	import SocketServer
	import threading

	class SearchRequestHandler(SocketServer.StreamRequestHandler):
		def __init__(self, *args, **kargs):
			SocketServer.StreamRequestHandler.__init__(self, *args,
				**kargs)

		def setup(self):
			sys.stderr.write( "Now Serving %s\n" % str(self.client_address))
			self.request.send("OK\n")
			SocketServer.StreamRequestHandler.setup(self)

		def handle(self):
			data = self.rfile.readline().strip()
			queue.put(data)
			s = queue.qsize()
			self.wfile.write("%d" % s)
	class ThreadedTCPServer(SocketServer.ThreadingMixIn, SocketServer.TCPServer):
		pass
	

	server = ThreadedTCPServer(('127.0.0.1', 50008), SearchRequestHandler)
	sys.stderr.write( "Now Serving\n" )
	t = threading.Thread(target=server.serve_forever)
	t.setDaemon(True)
	t.start()
	
if __name__ == '__main__':
	start_server()
	try:
		counter = int(file(COUNTER_FILE).read())
	except:
		counter = 0
	while True:
		counter += 1
		f = file(COUNTER_FILE, 'w')
		f.write('%d\n' % counter)
		f.close()
		dir = queue.get()
		if make(dir):
			sys.stderr.write("Successfully finished one query: %s\n" % dir)
		else:
			sys.stderr.write("Failed one query: %s\n" % dir)
