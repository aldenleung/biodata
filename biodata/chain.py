
from genomictools import StrandedGenomicPos
from .baseio import BaseReader
class Chain():
	'''
	The basic BED3 format, with only chrom, chromStart and chromEnd 
	'''
	__slots__ = 'score', 'tName', 'tSize', 'tStrand', 'tStart', 'tEnd', 'qName', 'qSize', 'qStrand', 'qStart', 'qEnd', 'id', 'datalist'
	def __init__(self, score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, id, datalist):
		self.score = score
		self.tName = tName
		self.tSize = tSize
		self.tStrand = tStrand
		self.tStart = tStart
		self.tEnd = tEnd
		self.qName = qName
		self.qSize = qSize
		self.qStrand = qStrand
		self.qStart = qStart
		self.qEnd = qEnd
		self.id = id
		self.datalist = datalist
	@property
	def t_stranded_genomic_pos(self):
		return StrandedGenomicPos(self.tName, self.tStart + 1, self.tEnd)
	@property
	def q_stranded_genomic_pos(self):
		return StrandedGenomicPos(self.qName, self.qStart + 1, self.qEnd)
	
class ChainAlignmentData():
	__slots__ = 'size', 'dt', 'dq'
	def __init__(self, size, dt, dq):
		self.size = size
		self.dt = dt
		self.dq = dq



class ChainReader(BaseReader):
	'''
	'''
	def _proceed_next_line(self):
		while True:
			line = self.f.readline()
			if line == '':
				self._line = None
				break
			line = line.rstrip("\r\n") # Auto-stripping
			if line != '': 
				self._line = line
				break
	
	def _read(self):
		if self._line is None:
			return None
		tmp = [func(d) for d, func in zip(self._line.split()[1:], [int, str, int, str, int, int, str, int, str, int, int, str])]
		
		self._proceed_next_line() 
		datalist = []
		while len(self._line.split()) != 1:
			datalist.append(ChainAlignmentData(*list(map(int, self._line.split()))))
			self._proceed_next_line() 
		datalist.append(ChainAlignmentData(int(self._line), None, None))
		self._proceed_next_line()
		
		return Chain(*tmp, datalist)
