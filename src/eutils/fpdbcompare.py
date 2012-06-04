"""
Calculating distances to all reference compounds
"""

import os
import sys
from pexpect import spawn, TIMEOUT
bin = 'ei-fp_db_compare_server'
exchange_fp = "/tmp/__db_compare.in"
from time import time
from signal import SIGINT

class DBComparerStartupError(Exception):
	pass

class ComparerError(Exception):
	pass

class DBComparer(object):
	"""The db comparer wrapper class"""
	def __init__(self, reference):
		"""start the process"""
		self.cmd = bin + " " + reference
		self.start()
	
	def start(self):
		try:
			self.child = spawn(self.cmd)
			self.child.logfile_read = sys.stderr
			self.child.expect_exact('ready', timeout=4)
		except TIMEOUT:
			try: self.child.terminate(force=True)
			except: pass
			sys.stderr.write("Error: Cannot start fp_db_compare_server\n")
			raise DBComparerStartupError
	
	def close(self):
		self.child.close(force=True)
	
	def tell(self, line):
		f = file(exchange_fp, 'w')
		if not line.endswith('\n'):
			line += '\n'
		f.write(line)
		f.close()
		self.child.kill(SIGINT)

	def compare(self, db_in, out):
		if not self.child.isalive():
			self.start()

		start = time()
		self.tell('%s\n%s\n' % (db_in, out))
		start = time()
		index = self.child.expect_exact(['OK:', 'Input:', TIMEOUT])
		if index == 0:
			fp = self.child.readline()
			return fp
		elif index == 1:
			raise ComparerError
		elif index == 2:
			self.close()
			self.start()
			raise DBComparerError
		else:
			self.close()
			self.start()
			raise DBComparerError

