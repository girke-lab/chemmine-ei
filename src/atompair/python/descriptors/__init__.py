import CDescriptors

class UnknownFormat(Exception):
	pass
class UnknownMode(Exception):
	pass
class ParserError(Exception):
	pass
class DatabaseIOError(Exception):
	pass
class ModeError(Exception):
	pass

class Descriptors:
	def __init__(self, input, format="sdfile"):
		if format == 'raw':
			assert isinstance(input, CDescriptors.Descriptors)
			self.content = input
			return
		self.content = CDescriptors.Descriptors()
		if format == 'sdfile':
			status = self.content.parse_sdfile(input)
		elif format == 'sdf':
			status = self.content.parse_sdf(input)
		elif format == 'smiles':
			status = self.content.parse_smiles(input)
		else:
			raise UnknownFormat
		if status == 0:
			raise ParserError
	
	def __len__(self):
		return self.content.get_len()
	
	def __getitem__(self, i):
		if i >= self.content.get_len():
			return None
		return self.content.get_descriptor(i)
	
	def __iter__(self):
		def _next():
			offset = 0
			while offset < len(self):
				yield self[offset]
				offset += 1
		return _next()
	
	def similarity_to(self, other):
		assert isinstance(other, Descriptors)
		return CDescriptors.similarity(self.content, other.content)

class Database:
	def __init__(self, filepath, mode="in"):
		self.mode = "closed"
		self.db = CDescriptors.Database()
		if mode == "in":
			status = self.db.open(filepath, 'r')
		elif mode == "out":
			status = self.db.open(filepath, 'w')
		elif mode == "append":
			status = self.db.open(filepath, 'a')
		elif mode == "preloaded":
			status = self.db.open(filepath, 'm')
		elif mode == "indexed":
			status = self.db.open(filepath, 'i')
		else:
			raise UnknownMode
		
		if status == 0:
			raise DatabaseIOError
		self.mode = mode
	
	def __del__(self):
		self.close()

	def close(self):
		self.db.close()
		self.mode = "closed"

	def next(self):
		if self.mode not in ["in", "preloaded", "indexed"]: raise ModeError
		d = self.db.next()
		if d.get_len() == 0:
			raise StopIteration
		_d = Descriptors(d, "raw")
		return _d

	def rewind(self):
		if self.mode not in ["in", "preloaded", "indexed"]: raise ModeError
		self.db.rewind()
		
	def __getitem__(self, i):
		if self.mode not in ["preloaded", "indexed"]: raise ModeError
		d = self.db.get(i)
		if d.get_len() == 0:
			raise IndexError
		_d = Descriptors(d, "raw")
		return _d

	def __iter__(self):
		if self.mode not in ["in", "preloaded", "indexed"]: raise ModeError
		return self
	
	def append(self, d):
		if self.mode == "closed": raise ModeError
		assert isinstance(d, Descriptors)
		if self.mode not in ["append", "out"]:
			raise ModeError
		self.db.store(d.content)

__all__ = ['Descriptors', 'Database']
