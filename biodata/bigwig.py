try:
	import pyBigWig
	_SUPPORT_BIGWIG = True
except:
	_SUPPORT_BIGWIG = False


from .baseio import BaseIReader
from genomictools import GenomicCollection, GenomicPos

class BigWigIReader(BaseIReader):
	def __init__(self, f, normalization=None, missing_data_mode="zero"):
		if not _SUPPORT_BIGWIG:
			raise Exception("bigwig is not supported without pyBigWig.")
		# missing
		self.bw = pyBigWig.open(f)
		
		if normalization is None:
			self.normalization_factor = 1
		elif isinstance(normalization, str):
			if normalization == "rpm":
				self. normalization_factor = 1000000/abs(self.bw.header()["sumData"])
			else:
				raise Exception("Unknown normalization method")
		elif isinstance(normalization, int) or isinstance(normalization, float):
			self.normalization_factor = normalization
		else:
			raise Exception()
		if missing_data_mode == "zero":
			self.missing_data_value = 0
		elif missing_data_mode == "nan":
			self.missing_data_value = float("nan")
		
	def value(self, r, method="sum"):
		'''
		Return a value processed across a selected region. It could be sum, max, abssum
		Some bigwig data also contain negative data and thus abs could be useful
		The value is always calculated exactly
		
		'''
		if r is None:
			raise Exception()
		elif isinstance(r, GenomicCollection):
			raise Exception()
		else:
			r = GenomicPos(r)
			if r.name not in self.bw.chroms():
				return 0
			zstart = r.start - 1
			ostop = r.stop
			intervals = self.bw.intervals(r.name, zstart, ostop)
			if method == "sum":
				if intervals is None:
					return 0
				return sum((min(i_ostop, ostop) - max(i_zstart, zstart)) * v for i_zstart, i_ostop, v in intervals) * self.normalization_factor
			elif method == "abssum":
				if intervals is None:
					return 0
				func = lambda vs: sum(abs(v) for v in vs)
				return func((min(i_ostop, ostop) - max(i_zstart, zstart)) * v for i_zstart, i_ostop, v in intervals) * self.normalization_factor
			elif method == "max":
				if intervals is None:
					return 0
				return func(v for i_zstart, i_ostop, v in intervals) * self.normalization_factor 
			elif method == "absmax":
				if intervals is None:
					return 0
				func = lambda vs: max(abs(v) for v in vs)
				return func(v for i_zstart, i_ostop, v in intervals) * self.normalization_factor 
			else:
				raise Exception()
	def values(self, r):
		'''
		Return a list of values of size length of r. Inefficient if the region is very large 
		'''
		r = r.genomic_pos
		d = self.values_dict(r)
		return [d[i+1] if i+1 in d else self.missing_data_value for i in range(r.zstart, r.ostop)]

	def values_dict(self, r):
		'''
		Return a dict of values. The key of dict is 1-based coordinate. Missing data is not put in the dictionary
		'''
		r = GenomicPos(r)
		if r.name not in self.bw.chroms():
			return {}
		zstart = r.zstart
		ostop = r.ostop
		intervals = self.bw.intervals(r.name, zstart, ostop)
		if intervals is None:
			return {}
		return {p+1 : v * self.normalization_factor for i_zstart, i_ostop, v in intervals for p in range(max(i_zstart, zstart), min(i_ostop, ostop))}
		
	
	def close(self):
		self.bw.close()