"""
Performing LSH search (without refinement)
"""

import os
import sys
from pexpect import spawn, TIMEOUT
bin = 'ei-search_server'
matrix =  'matrix'
exchange_fp = 'coord.in'
import signal
from time import time
from traceback import print_exc

class LSHSearcherStartupError(Exception):
	pass

class LSHSearchError(Exception):
	pass

class LSHSearcher(object):
	"""The fingerprint calculator wrapper class"""
	def __init__(self, matrix, w=0.65664, m=16, l=10, k=1000, t=80): # l=30
		"""start the process"""
		sys.stderr.write("Using parameters:\n"
			" -W %s -M %s -L %s -K %s -T %s -D %s\n" % (w, m, l, k, t, matrix)
			)
		self.cmd = bin + " -W %s -M %s -L %s -K %s -T %s -D %s" % (
			w, m, l, k, t, matrix)
		self.start()
	
	def tell(self, msg):
		f = file(exchange_fp, 'w')
		if not msg.endswith('\n'):
			msg += '\n'
		f.write(msg)
		f.close()
		self.child.kill(signal.SIGINT)
		
	def start(self):
		try:
			self.child = spawn(self.cmd)
			self.child.logfile_read = sys.stderr
			self.child.expect_exact('>>', timeout=4000)
		except TIMEOUT:
			try: self.child.terminate(force=True)
			except: pass
			sys.stderr.write("Error: Cannot start lsh program\n") 
			raise LSHSearcherStartupError
	
	def close(self):
		self.child.tell('/Q\n')

	def search(self, line):
		self.tell(line)
		index = self.child.expect_exact(['OK:', TIMEOUT], timeout=10)
		if index == 0:
			fp = self.child.readline()
			return fp
		elif index == 1:
			raise LSHSearchError
		else:
			self.child.terminate(force=True)
			self.start()
			raise LSHSearchError

if __name__ == '__main__':
	sys.stderr.write( "creating index..." )
	lshsearcher = LSHSearcher(matrix)
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
				start = time()
				ret = lshsearcher.search(data)
				elapsed = time() - start
			except:
				print_exc()
				elapsed = 0
				ret = '\n'
			self.wfile.write('/t:%f %s' % (elapsed, ret))

	server = SocketServer.TCPServer(('0.0.0.0', 50009), SearchRequestHandler)
	sys.stderr.write( "Now Serving\n" )
	server.serve_forever()

