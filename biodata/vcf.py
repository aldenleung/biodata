'''
Created on Dec 27, 2020

@author: Alden
'''
from builtins import int
'''
Created on Sep 28, 2018

@author: Alden
'''

'''
Created on May 18, 2018

@author: Alden
'''
from collections import OrderedDict

from genomictools import GenomicAnnotation, GenomicPos
from .baseio import BaseReader, BaseWriter
from .tabix import TabixIReader
import csv
import io
def _csv_split(s, delimiter):
	return next(iter(csv.reader(io.StringIO(s), delimiter=delimiter, quotechar='"')))


def _comma_split(value):
	open_quote = False
	ps = [-1]
	for p, c in enumerate(value):
		if c == '"':
			open_quote = not open_quote
		elif c == ",":
			if not open_quote:
				ps.append(p)
	ps.append(len(value))
	split_values = [value[ps[i] + 1:ps[i + 1]] for i in range(len(ps) - 1)]
	return split_values

class VCF(GenomicAnnotation):
	__slots__ = "chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "genotypes"	
	def __init__(self, chrom, pos, id, ref, alt, qual, filter, info, format=None, genotypes=None):
		self.chrom = chrom
		self.pos = pos
		self.id = id
		self.ref = ref
		self.alt = alt
		self.qual = qual
		self.filter = filter
		self.info = info
		self.format = format
		self.genotypes = genotypes
	@property
	def genomic_pos(self):
		return GenomicPos(self.chrom, self.pos)
		
class VCFInfo():
	__slots__ = "ID", "Number", "Type", "Description", "Source", "Version", "_typefunc"
	def __init__(self, ID=None, Number=None, Type=None, Description=None, Source=None, Version=None):
		self.ID = ID
		self.Number = Number
		self.Type = Type
		self.Description = Description
		self.Source = Source
		self.Version = Version
		self._compile_type_func()
		
	def _compile_type_func(self):
		if self.Type == "Integer":
			ntype = int
		elif self.Type == "Float":
			ntype = lambda a: float('NaN') if a == "." else float(a)
		else:
			ntype = str
		
		if self.Number == "A":
			self._typefunc = lambda s: list(map(ntype, s.split(",")))
		elif self.Number == "R":
			self._typefunc = lambda s: list(map(ntype, s.split(",")))
		elif self.Number == "G":
			self._typefunc = lambda s: list(map(ntype, s.split(",")))
		elif self.Number == ".":
			self._typefunc = lambda s: list(map(ntype, s.split(",")))
		else:
			n = int(self.Number)
			if n == 1:
				self._typefunc = ntype
			else:
				self._typefunc = lambda s: list(map(ntype, s.split(",")))
		 
	@property
	def typefunc(self):
		return self._typefunc
		
		
def _parse_words_array_VCF(vr, words_array):
	
	chrom = words_array[0]
	pos = int(words_array[1])
	id = words_array[2]
	ref = words_array[3]
	alt = words_array[4]
	qual = words_array[5]
	filter = words_array[6]
	info = OrderedDict()
	for mark in words_array[7].split(";"):
		if "=" in mark:
			tk, tv = mark.split("=", maxsplit=1)
			if tk in vr.metainfo["INFO"]:
				try:
					info[tk] = vr.metainfo["INFO"][tk].typefunc(tv)
				except:
					print(tk, words_array)
			else:
				info[tk] = tv
		else:
			info[mark] = None # Flag
	if vr.is_genotype_info_available:
		format = words_array[8]
		sample_fields = words_array[9:]
		genotypes = OrderedDict(zip(vr.samples, sample_fields))
	else:
		format = None
		genotypes = None
	return VCF(chrom, pos, id, ref, alt, qual, filter, info, format, genotypes)		
class VCFReader(BaseReader):	
	__slots__ = "metainfo", "samples", "is_genotype_info_available"	
	def _metainfo_reader(self):
		self._proceed_next_line()
		self.is_genotype_info_available = False
		self.metainfo = {}
		self.metainfo["INFO"] = {}
		while self._line is not None and self._line.startswith("##"):
			line = self._line[2:]
			idx = line.index("=")
			key = line[:idx]
			value = line[idx + 1:]
			if key == "INFO":
				# Split comma not in quote
				value = value[1:-1]
				open_quote = False
				ps = [-1]
				for p, c in enumerate(value):
					if c == '"':
						open_quote = not open_quote
					elif c == ",":
						if not open_quote:
							ps.append(p)
				ps.append(len(value))
				split_values = [value[ps[i] + 1:ps[i + 1]] for i in range(len(ps) - 1)]
				vcfinfo = VCFInfo(**OrderedDict([_csv_split(i, "=") for i in split_values]))
				self.metainfo[key][vcfinfo.ID] = vcfinfo 
				# Type Integer, Float, Flag, Character, and String
				pass
			elif key == "FILTER":				
				pass
			elif key == "FORMAT":
				pass
			elif key == "contig":
				value = value[1:-1]
# 				d = dict(OrderedDict([_csv_split(i, "=") for i in _comma_split(value)]))
# 				self.metainfo[key][d["ID"]] 
# 			self.metainfo[key] = value
			self._proceed_next_line()
		if self._line is not None and self._line.startswith("#"):
			
			header_fields = self._line.split('\t')
			if len(header_fields) > 8:
				# GENOTYPE NAME
				self.samples = header_fields[9:]
				self.is_genotype_info_available = True
			else:
				self.samples = None
				self.is_genotype_info_available = False
			self._proceed_next_line()
		
	def _read(self):
		line = self._line
		if line is None:
			return None
		self._proceed_next_line() 
		words_array = line.split('\t')
		return _parse_words_array_VCF(self, words_array)
	
# class VCFIReader(TabixIReader):
# 	def __init__(self, arg, tbi=None, ):
# 		super(VCFIReader, self).__init__(arg, tbi)
# 		with VCFReader(arg) as vr:
# 			self.metainfo = vr.metainfo
# 			self.samples = vr.samples
# 			self.is_genotype_info_available = vr.is_genotype_info_available
# 	def _parse_raw_entry(self, entry):
# 		return _parse_words_array_VCF(self, entry)	
# 
# 	
# class VCFWriter(BaseWriter):
# 	__slots__ = "metainfo"
# 	def __init__(self, arg, metainfo=None):
# 		raise Exception("Incomplete function")
# 		self.metainfo = metainfo
# 		super(VCFWriter, self).__init__(arg)
# 	def _initialize_header(self):
# 		pass
# 	def write(self, vcf):
# 		pass
# 	
