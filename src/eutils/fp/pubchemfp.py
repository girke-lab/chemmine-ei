"""
Calculating Pubchem Fingerprint/Keys using a Java program. This is a wraper.
"""

import os
from pexpect import spawn, TIMEOUT
lib_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lib')
bin_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'bin')
java_class = "PubChemFp"
lib_jar = os.path.join(lib_dir, 'jchem.jar')

class FPCalculatorStartupError(Exception):
	pass

class FPCalculationError(Exception):
	pass

class FPCalculator(object):
	"""The fingerprint calculator wrapper class"""
	def __init__(self):
		"""start the process"""
		self.cmd = "java -cp %s:. %s" % (lib_jar, java_class)
		self.start()
	
	def start(self):
		try:
			self.child = spawn(self.cmd, cwd=bin_dir)
			self.child.expect_exact('>>>', timeout=10)
		except TIMEOUT:
			try: self.child.terminate(force=True)
			except: pass
			sys.stderr.write("Error: Cannot start java program to calculate "
				"fingerprint.\n")
			raise FPCalculatorStartupError

	def close(self):
		self.child.close(force=True)

	def calc(self, sdf):
		sdf = sdf.replace('\r', '')
		sdf_end = sdf.find('$$$$')
		if sdf_end != -1:
			sdf = sdf[:sdf_end]
		sdf += '\n$$$$\n'

		self.child.send(sdf)
		index = self.child.expect_exact(['OK:', '**', TIMEOUT])
		if index == 0:
			fp = self.child.readline()
			return fp
		elif index == 1:
			raise FPCalculationError
		else:
			self.close()
			self.start()
			raise FPCalculationError

fpcalc = None
def getFpcalc():
	global fpcalc
	if fpcalc is None:
		fpcalc = FPCalculator()
	return fpcalc
