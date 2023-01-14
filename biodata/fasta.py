'''
Created on Dec 29, 2020

@author: Alden
'''
import textwrap
from .baseio import BaseReader, BaseWriter, BaseIReader
from genomictools import GenomicAnnotation, StrandedGenomicAnnotation

_reverse_complment_map = {
"A":"T",
"G":"C",
"C":"G",
"T":"A",
"Y":"R",
"R":"Y",
"W":"W",
"S":"S",
"K":"M",
"M":"K",
"D":"H",
"V":"B",
"H":"D",
"B":"V",
"X":"X",
"N":"N",
"a":"t",
"g":"c",
"c":"g",
"t":"a",
"y":"r",
"r":"y",
"w":"w",
"s":"s",
"k":"m",
"m":"k",
"d":"h",
"v":"b",
"h":"d",
"b":"v",
"x":"x",
"n":"n",
"-":"-"
}

def _reverse_complement(s):
	return "".join([_reverse_complment_map[c] for c in s][::-1])

class FASTA(object):
	'''
	The basic FASTA format
	'''
	
	def __init__(self, name, seq):
		self.name = name
		self.seq = seq
	
	def __len__(self):
		return len(self.seq)
	
class FASTAReader(BaseReader):
	'''
	'''
	def _read(self):
		if self._line is None:
			return None
		words_array = self._line.split('>')
		name = words_array[1]
		self._proceed_next_line()
		
		seqarray = []
		while self._line is not None and not self._line.startswith(">"):
			seqarray.append(self._line)
			self._proceed_next_line()	
		seq = "".join(seqarray)
				
		return FASTA(name, seq)
	
	
class Faidx():
	'''
	File format specification from htslib
	
	NAME	Name of this reference sequence
	LENGTH	Total length of this reference sequence, in bases
	OFFSET	Offset in the FASTA/FASTQ file of this sequence's first base
	LINEBASES	The number of bases on each line
	LINEWIDTH	The number of bytes in each line, including the newline
	QUALOFFSET	Offset of sequence's first quality within the FASTQ file
	'''
	def __init__(self, name, length, offset, linebases, linewidth, qualoffset):
		self.name = name
		self.length = length
		self.offset = offset
		self.linebases = linebases
		self.linewidth = linewidth
		self.qualoffset = qualoffset

class FaidxReader(BaseReader):
	def _read(self):
		if self._line is None:
			return None
		words_array = self._line.split()
		self._proceed_next_line()
		
		name = words_array[0]
		length = int(words_array[1])
		offset = int(words_array[2])
		linebases = int(words_array[3])
		linewidth = int(words_array[4])
		if len(words_array) > 5:
			qualoffset = int(words_array[5])
		else:
			qualoffset = None
		
		return Faidx(name, length, offset, linebases, linewidth, qualoffset)
	
class FASTAIReader(BaseIReader):
	def __init__(self, f, fai=None):
		if isinstance(f, str):
			if fai is None:
				fai = f + ".fai"
		if fai is None:
			raise Exception("Cannot auto-determine fai file")
		with FaidxReader(fai) as fr:
			self.faidx_dict = {idx.name:idx for idx in fr}
		self.f = open(f, "rb")
	def __getitem__(self, key):
		if isinstance(key, StrandedGenomicAnnotation):
			r = key.stranded_genomic_pos
			faidx = self.faidx_dict[r.name]
			if r.name not in self.faidx_dict:
				raise Exception("Seq name not found")
			faidx = self.faidx_dict[r.name]
			if r.zstart < 0 or r.ostop > faidx.length:
				raise Exception("Invalid range")
			start_offset = faidx.offset + r.zstart // faidx.linebases * faidx.linewidth + r.zstart % faidx.linebases
			stop_offset = faidx.offset + r.ostop // faidx.linebases * faidx.linewidth + r.ostop % faidx.linebases
			
			self.f.seek(start_offset)
			s = self.f.read(stop_offset - start_offset).decode()
			seq = "".join([c for c in s if c not in ["\r", "\n"]])
			if r.strand == "-":
				seq = _reverse_complement(seq)
			return FASTA(str(r), seq)
		elif isinstance(key, GenomicAnnotation):
			r = key.genomic_pos
			faidx = self.faidx_dict[r.name]
			if r.name not in self.faidx_dict:
				raise Exception("Seq name not found")
			faidx = self.faidx_dict[r.name]
			if r.zstart < 0 or r.ostop > faidx.length:
				raise Exception("Invalid range")
			start_offset = faidx.offset + r.zstart // faidx.linebases * faidx.linewidth + r.zstart % faidx.linebases
			stop_offset = faidx.offset + r.ostop // faidx.linebases * faidx.linewidth + r.ostop % faidx.linebases
			
			self.f.seek(start_offset)
			s = self.f.read(stop_offset - start_offset).decode()
			return FASTA(str(r), "".join([c for c in s if c not in ["\r", "\n"]]))			
		else:
			raise Exception("Unknown key format")
		
	def close(self):
		self.f.close()
		

class FASTAWriter(BaseWriter):
	def __init__(self, arg, maxcharac=-1):
		super(FASTAWriter, self).__init__(arg)
		self.maxcharac = maxcharac
	
	def write(self, fa):	
		self.f.write(">" + fa.name + "\n");
		if self.maxcharac > 0:
			s = "\n".join(textwrap.wrap(fa.seq, self.maxcharac, drop_whitespace=False, break_on_hyphens=False))
		else:
			s = fa.seq
		self.f.write(s + "\n");
		
class FASTQ(FASTA):
	def __init__(self, name, seq, quality):
		super(FASTQ, self).__init__(name, seq)
		if len(seq) != len(quality):
			raise Exception("Mismatch in length of sequence and quality")
		self.quality = quality

class FASTQReader(BaseReader):
	def _proceed_next_line(self):
		line = self.f.readline()
		if line == '':
			self._line = None
		else:
			line = line.rstrip("\r\n") # Auto-stripping
			self._line = line
			
	def _read(self):
		if self._line is None:
			return None
		words_array = self._line.split('@')
		name = words_array[1]
		self._proceed_next_line()
		seq = self._line		
		self._proceed_next_line()
		# "+"
		self._proceed_next_line()
		quality = self._line
		self._proceed_next_line()
		return FASTQ(name, seq, quality)
	

class FASTQWriter(BaseWriter):
	
	def initialize_header(self):
		pass
	
	def write(self, fa):	
		self.f.write("@" + fa.name + "\n");
		self.f.write(fa.seq + "\n");
		self.f.write("+\n");
		self.f.write(fa.quality + "\n");
	