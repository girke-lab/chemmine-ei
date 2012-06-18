"""
Calculating coordinate from distances to references
"""

import os
import sys
from pexpect import spawn, TIMEOUT
bin =  'ei-coord_server'
exchange_fp = "/tmp/__coord_server.in"
from time import time
from signal import SIGINT

class CoordinateSolverStartupError(Exception):
	pass

class SolverError(Exception):
	pass

class CoordinateSolver(object):
	"""The fingerprint calculator wrapper class"""
	def __init__(self, reference):
		"""start the process"""
		self.cmd = bin + " " + reference
		self.start()
	
	def start(self):
		try:
			self.child = spawn(self.cmd)
			#self.child.logfile_read = sys.stderr
			self.child.expect_exact('>>', timeout=40)
		except TIMEOUT:
			try: self.child.terminate(force=True)
			except: pass
			sys.stderr.write("Error: Cannot start icoord\n")
			raise CoordinateSolverStartupError
	
	def close(self):
		self.child.close(force=True)
	
	def tell(self, line):
		f = file(exchange_fp, 'w')
		if not line.endswith('\n'):
			line += '\n'
		f.write(line)
		self.child.kill(SIGINT)
		f.close()

	def solve(self, line):
		if not self.child.isalive():
			self.start()

		start = time()
		self.tell(line)
		start = time()
		index = self.child.expect_exact(['OK:', 'Input: File I/O failed when reading first line', TIMEOUT])
		if index == 0:
			fp = self.child.readline()
			return fp
		elif index == 1:
			print "before: "+self.child.before
			print "after: "+self.child.after
			raise SolverError
		elif index == 2:
			self.close()
			self.start()
			raise SolverError
		else:
			self.close()
			self.start()
			raise SolverError

