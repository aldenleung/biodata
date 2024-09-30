import os

from .baseio import BaseIReader
from genomictools import GenomicPos

try:
	import tabix
	_SUPPORT_TABIX = True
except:
	_SUPPORT_TABIX = False

class TabixIReader(BaseIReader):
	def __init__(self, arg, tbi=None):
		if not _SUPPORT_TABIX:
			raise Exception("TabDelimitedGenomicIReader is not supported without pytabix")
		if isinstance(arg, str):
			if tbi is None:
				tbi = arg + ".tbi"
		if tbi is None:
			raise Exception("Cannot auto-determine tbi file")
		if not os.path.exists(arg) or not os.path.exists(tbi):
			raise Exception("Path not exists")
		self.tb = tabix.open(arg, tbi)
		
	def _parse_raw_entry(self, entry):
		return entry
	def entries_iterator(self, r):
		r = GenomicPos(r)
		for entry in self.tb.query(r.name, r.zstart, r.ostop):
			yield self._parse_raw_entry(entry)
	def entries(self, r):
		return list(self.entries_iterator(r))
	def __getitem__(self, key):
		return self.entries[key]
	
	def close(self):
		pass
