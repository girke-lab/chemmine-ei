"""
perform chemical search
"""

import os
import sys
from pexpect import spawn, TIMEOUT
bin = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fp_search_server')
db = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'db.fp_cdb')
q_db = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'query.fp_cdb')
candidates_f = os.path.join(os.path.dirname(os.path.abspath(__file__)),
	'candidates.data')
from base64 import b64decode
from time import time
import signal
from traceback import print_exc

class RefinerStartupError(Exception):
	pass

class RefinerError(Exception):
	pass

class Refiner(object):
	"""The fingerprint calculator wrapper class"""
	def __init__(self, cdb, k):
		"""start the process"""
		self.cmd = "%s %s %s" % (bin, cdb, k)
		self.start()
	
	def start(self):
		try:
			self.child = spawn(self.cmd)
			self.child.logfile_read = sys.stderr
			self.child.expect_exact('>>', timeout=1000)
		except TIMEOUT:
			try: self.child.terminate(force=True)
			except: pass
			sys.stderr.write("Error: Cannot start fp_db_isearch program\n")
			raise RefinerStartupError
	
	def close(self):
		self.child.close(force=True)

	def refine(self, line):
		if not self.child.isalive():
			self.start()
		if not line.endswith('\n'):
			line += '\n'
		start = time()
		#self.child.send(line)
		self.child.kill(signal.SIGINT)
		sys.stderr.write("sending took %f seconds\n" % (time() - start))
		start = time()
		index = self.child.expect_exact(['OK:', 'Cannot open', TIMEOUT])
		sys.stderr.write("expecting took %f seconds\n" % (time() - start))
		if index == 0:
			fp = self.child.readline()
			return fp
		elif index == 1:
			raise RefinerError
		elif index == 2:
			self.close()
			self.start()
			raise RefinerError
		else:
			self.close()
			self.start()
			raise RefinerError

if __name__ == '__main__':
	sys.stderr.write( "loading database..." )
	refiner = Refiner(db, 200)
	sys.stderr.write("done\n")

	import SocketServer

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
			try:
				c, candidates = data.split(None, 1)
			except ValueError:
				c = data.strip()
				candidates = ''
			f = file(q_db, 'w')
			f.write(b64decode(c))
			f.close()
			if candidates:
				f = file(candidates_f, 'w')
				f.write(candidates)
				line = "%s %s" % (q_db, candidates_f)
			else:
				line = q_db
			try:
				start = time()
				ret = refiner.refine(line)
				elapsed = time() - start
				sys.stderr.write("search took %f seconds\n" % elapsed)
			except:
				print_exc()
				elapsed = 0
				ret = '\n'
			self.wfile.write('/t:%f %s' % (elapsed, ret))

	server = SocketServer.TCPServer(('0.0.0.0', 50010), SearchRequestHandler)
	sys.stderr.write( "Now Serving\n" )
	server.serve_forever()


