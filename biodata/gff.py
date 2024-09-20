import csv
from collections import OrderedDict
import logging


from genomictools import StrandedGenomicAnnotation, StrandedGenomicPos
from .baseio import BaseReader, BaseWriter
from .tabix import TabixIReader
def _quote_split(s, delimiter):
	return next(csv.reader([s],delimiter=delimiter))
		
	

__all__ = ["GFF3", "GFF3Reader", "GFF3Writer", "GTF", "GTFReader", "GTFWriter"]
class GFF3(StrandedGenomicAnnotation):
	def __init__(self, seqid, source, feature, start, end, score, strand, phase, attributes):
		super(GFF3, self).__init__()
		self.seqid = seqid
		self.source = source
		self.feature = feature
		self.start = start
		self.end = end
		self.score = score
		self.strand = strand
		self.phase = phase		
		self.attributes = attributes
	
	# For support to old GTF field name
	@property
	def seqname(self):
		return self.seqid	
	@property
	def attribute(self):
		return self.attributes
	@property
	def frame(self):
		return self.phase
	
	@property
	def stranded_genomic_pos(self):
		return StrandedGenomicPos(self.seqname, self.start, self.end, self.strand)
	
	
def _parse_words_array_GFF3(words_array):
	seqid = words_array[0]
	source = words_array[1]
	feature = words_array[2]
	start = int(words_array[3])
	end = int(words_array[4])
	score = float(words_array[5]) if words_array[5] != "." else None
	strand = words_array[6]
	phase = words_array[7]
	attributes = OrderedDict([_quote_split(item.strip(), "=") for item in _quote_split(words_array[8], ";")])
	return GFF3(seqid, source, feature, start, end, score, strand, phase, attributes)
class GFF3Reader(BaseReader):
	def __init__(self, arg):
		if isinstance(arg, str):
			super(GFF3Reader, self).__init__(arg)
	def _metainfo_reader(self):
		line = self.f.readline()
		line = line.rstrip("\r\n") # Auto-stripping
		if line is not None:
			if line != "##gff-version 3":
				logging.warning("Missing GFF3-version header. This file may not be a proper GFF3 file.")
				# This warning may be modified in the future
			self._proceed_next_line()
		else:
			raise Exception("The first line must be ##gff-version 3.")
		
	def _proceed_next_line(self):
		'''
		Attempts to assign the next line to self._line. 
		'''
		while True:
			line = self.f.readline()
			if line == '':
				self._line = None
				break
			line = line.rstrip("\r\n") # Auto-stripping
			if line != '' and not line.startswith("#"): # Auto comment removal. Blank lines are skipped by default
				self._line = line
				break
		
	def _read(self):
		line = self._line
		if line is None:
			return None
		self._proceed_next_line()
		words_array = line.split('\t')
		return _parse_words_array_GFF3(words_array)
	
	@staticmethod
	def read_all_by_features(features, read_all_func, *args, **kwargs):
		'''
		A short cut method to read GFFs with selected features
		'''
		return GFF3Reader.read_all(lambda gff3_generator: read_all_func(filter(lambda gff3: gff3.feature in features, gff3_generator)), *args, **kwargs)
class GFF3IReader(TabixIReader):
	def _parse_raw_entry(self, entry):
		return _parse_words_array_GFF3(entry)	

class GFF3Writer(BaseWriter):
	def __init__(self, arg):
		if isinstance(arg, str):
			super(GFF3Writer, self).__init__(arg)
	def _initialize_header(self):
		super(GFF3Writer, self).write("##gff-version 3\n")
	def write(self, gff3):
		'''
		Output the GFF
		Note that the start position is stored in 1-based
		'''
		super(GFF3Writer, self).write(("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n").format(
			gff3.seqid,
			gff3.source,
			gff3.feature,
			gff3.start,
			gff3.end,
			gff3.score if gff3.score is not None else ".",
			gff3.strand,
			gff3.phase,	
			";".join([key + "=" + value for key, value in gff3.attributes.items()])))
			


class GTF(StrandedGenomicAnnotation):
	# Equivalent to GFF2
	def __init__(self, seqname, source, feature, start, end, score, strand, frame, attribute):
		self.seqname = seqname
		self.source = source
		self.feature = feature
		self.start = start
		self.end = end
		self.score = score
		self.strand = strand
		self.frame = frame		
		self.attribute = attribute

	# For support to GFF3 field name
	@property
	def seqid(self):
		return self.seqname	
	@property
	def attributes(self):
		return self.attribute
	@property
	def phase(self):
		return self.frame
	
	@property
	def stranded_genomic_pos(self):
		return StrandedGenomicPos(self.seqname, self.start, self.end, self.strand)
	
def _parse_words_array_GTF(words_array):
	seqname = words_array[0]
	source = words_array[1]
	feature = words_array[2]
	start = int(words_array[3])
	end = int(words_array[4])
	score = float(words_array[5]) if words_array[5] != "." and words_array[5] != "None" else None
	strand = words_array[6]
	frame = words_array[7]
	attribute = OrderedDict([_quote_split(item.strip(), " ") for item in _quote_split(words_array[8], ";") if len(item.strip()) != 0]) # Defective
	return GTF(seqname, source, feature, start, end, score, strand, frame, attribute)

class GTFReader(BaseReader):
	def _proceed_next_line(self):
		'''
		Attempts to assign the next line to self._line. 
		'''
		while True:
			line = self.f.readline()
			if line == '':
				self._line = None
				break
			line = line.rstrip("\r\n") # Auto-stripping
			if line != '' and not line.startswith("#"): # Auto comment removal. Blank lines are skipped by default
				self._line = line
				break
	
	def _read(self):
		line = self._line
		if line is None:
			return None
		self._proceed_next_line()
		words_array = line.split('\t')
		return _parse_words_array_GTF(words_array)
class GTFIReader(TabixIReader):
	def _parse_raw_entry(self, entry):
		return _parse_words_array_GTF(entry)	

class GTFWriter(BaseWriter):
	def __init__(self, arg):
		super(GTFWriter, self).__init__(arg)
			
	def write(self, gtf):
		'''
		Output the GTFNode
		'''
		super(GTFWriter, self).write(("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n").format(
			gtf.seqname,
			gtf.source,
			gtf.feature,
			gtf.start,
			gtf.end,
			gtf.score if gtf.score is not None else ".",
			gtf.strand,
			gtf.frame,	
			"; ".join([key + " " + "\"" + value + "\"" for key, value in gtf.attribute.items()])))
