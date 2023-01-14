'''
Created on Sep 21, 2021

@author: kl945
'''

import re
from genomictools import StrandedGenomicPos
from ..baseio import BaseReader

		
class FastMap():
	def __init__(self, qname, qlen, matches):
		self.qname = qname
		self.qlen = qlen
		self.matches = matches
		
class FastMapEM():
	'''
	qstart is 0-based
	qstop is 1-based
	'''
	def __init__(self, qstart, qstop, nmatch, submatches):
		self.qstart = qstart
		self.qstop = qstop
		self.nmatch = nmatch
		self.submatches = submatches
	@property
	def alignlen(self):
		return self.qstop - self.qstart

class FastMapReader(BaseReader):
	'''
	'''
	def _proceed_next_line(self):
		while True:
			line = self.f.readline()
			if line == '':
				self._line = None
				break
			line = line.rstrip("\r\n") # Auto-stripping
			line = line.strip()
			if line != '': # Auto comment removal. Blank lines are skipped by default
				self._line = line
				break
	
	def _read(self):
		line = self._line
		if line is None:
			return None
		while not line.startswith("SQ"):
			self._proceed_next_line() 
			line = self._line
			if line is None:
				return None
		
		_, qname, qlen = line.split("\t")
		qlen = int(qlen)
		self._proceed_next_line()
		line = self._line
		
		matches = []
		while not line.startswith("//"):
			qstart, qstop, nmatch, tmp_match_str = re.match(r"^EM\s([0-9]+)\s([0-9]+)\s([0-9]+)\s(.+)$", line).groups()
			qstart = int(qstart)
			qstop = int(qstop)
			nmatch = int(nmatch)
			qalignlen = qstop - qstart
			if tmp_match_str == "*":
				submatches = None
			else:
				submatches = []
				for i in tmp_match_str.split("\t"):
					rname, strand, rstart = re.match(r"^(.+):(\+|-)([0-9]+)$", i).groups()
					rstart = int(rstart)
					submatches.append(StrandedGenomicPos(rname, rstart, rstart + qalignlen - 1, strand))
			matches.append(FastMapEM(qstart, qstop, nmatch, submatches))
			self._proceed_next_line()
			line = self._line
		self._proceed_next_line()		
		return FastMap(qname, qlen, matches)
	
	
