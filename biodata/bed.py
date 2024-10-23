try:
	import pyBigWig
	_SUPPORT_BIGWIG = True
except:
	_SUPPORT_BIGWIG = False


# BED3, BED, BEDX, BEDGraph, BEDPE and all bed-derived file formats

from .baseio import BaseReader, BaseWriter, BaseIReader
from genomictools import GenomicAnnotation, StrandedGenomicAnnotation, GenomicPos, StrandedGenomicPos
from .tabix import TabixIReader
class BED3(GenomicAnnotation):
	'''
	The basic BED3 format, with only chrom, chromStart and chromEnd 
	'''
	__slots__ = 'chrom', 'chromStart', 'chromEnd'
	def __init__(self, chrom, chromStart, chromEnd):
		self.chrom = chrom
		self.chromStart = chromStart
		self.chromEnd = chromEnd
# 		
	@property
	def genomic_pos(self):
		return GenomicPos(self.chrom, self.chromStart + 1, self.chromEnd)
	
	@staticmethod
	def from_genomic_pos(r):
		return BED3(r.name, r.zstart, r.ostop)

class BED3Reader(BaseReader):
	'''
	'''
	def _is_empty_field(self, s):
		return s == "" or s == "."
	def _proceed_next_line(self):
		while True:
			line = self.f.readline()
			if line == '':
				self._line = None
				break
			line = line.rstrip("\r\n") # Auto-stripping
			if line != '' and not line.startswith("#"): # Auto comment removal. Blank lines are skipped by default
				self._line = line
				break
	def _metainfo_reader(self):
		self._proceed_next_line()
		while self._line is not None and (self._line.startswith("track") or self._line.startswith("browser")):
			self._proceed_next_line()
	
	def _read(self):
		line = self._line
		if line is None:
			return None
		self._proceed_next_line() 
		words_array = line.split('\t')
		chrom = words_array[0]
		chromStart = int(words_array[1]) 
		chromEnd = int(words_array[2])
		return BED3(chrom, chromStart, chromEnd)
	
class BED3Writer(BaseWriter):
	def write(self, bed3):
		'''
		Output the BED3
		'''
		if isinstance(bed3, BED3):
			super(BED3Writer, self).write("{}\t{}\t{}\n".format(bed3.chrom, bed3.chromStart, bed3.chromEnd))
		else: # assume GenomicPos
			r = bed3.genomic_pos
			super(BED3Writer, self).write("{}\t{}\t{}\n".format(r.name, r.zstart, r.ostop))


		
	
class BED(BED3, StrandedGenomicAnnotation):
	'''
	The standard 12-field BED format 
	'''
	__slots__ = "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"
	def __init__(self, chrom, chromStart, chromEnd, name=None, score=None, strand=None, thickStart=None, thickEnd=None, itemRgb=None, blockCount=None, blockSizes=None, blockStarts=None):
		super(BED, self).__init__(chrom, chromStart, chromEnd)
		self.name = name
		self.score = score
		self.strand = strand
		self.thickStart = thickStart
		self.thickEnd = thickEnd
		self.itemRgb = itemRgb
		self.blockCount = blockCount
		self.blockSizes = blockSizes
		self.blockStarts = blockStarts
	
	@property
	def stranded_genomic_pos(self):
		return StrandedGenomicPos(self.chrom, self.chromStart + 1, self.chromEnd, "." if self.strand is None else self.strand)
class BEDReader(BED3Reader):

	def _is_empty_field(self, s):
		return s is None or s == "" or s == "."
	
	def _read(self):
		line = self._line
		if line is None:
			return None
		self._proceed_next_line() 
		words_array = line.split('\t')
		chrom = words_array[0]
		chromStart = int(words_array[1]) 
		chromEnd = int(words_array[2])
		words_array.extend([None] * (12 - len(words_array)))
		
		name = words_array[3]
		score = None
		strand = None
		thickStart = None
		thickEnd = None
		itemRgb = None
		blockCount = None
		blockSizes = None
		blockStarts = None
		
		if not self._is_empty_field(words_array[4]):
			score = float(words_array[4])
		if words_array[5] != "":
			strand = words_array[5]
		if not self._is_empty_field(words_array[6]):
			thickStart = words_array[6]
		if not self._is_empty_field(words_array[7]):
			thickEnd = words_array[7]
		if not self._is_empty_field(words_array[8]):
			itemRgb = words_array[8]
		if not self._is_empty_field(words_array[9]):
			blockCount = int(words_array[9])
		if not self._is_empty_field(words_array[10]):
			blockSizes = list(map(int, filter(lambda a:a, words_array[10].split(","))))
		if not self._is_empty_field(words_array[11]):
			blockStarts = list(map(int, filter(lambda a:a, words_array[11].split(","))))
		
		if blockSizes is not None:
			if blockCount != len(blockSizes):
				raise Exception("Inconsistent blockCount and blockSizes")
			if blockCount != len(blockStarts):
				raise Exception("Inconsistent blockCount and blockStarts")
			if blockStarts[0] != 0:
				raise Exception("First block must start at 0")
			if chromStart + blockStarts[-1] + blockSizes[-1] != chromEnd:
				raise Exception("Last block must end at chromEnd")		 
		
		return BED(chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts)

	
	
class BEDWriter(BaseWriter):
	def write(self, bed):
		'''
		Output the BEDNode
		'''
		super(BEDWriter, self).write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(bed.chrom, bed.chromStart, bed.chromEnd, "" if bed.name is None else bed.name, "" if bed.score is None else bed.score, "" if bed.strand is None else bed.strand, "" if bed.thickStart is None else bed.thickStart, "" if bed.thickEnd is None else bed.thickEnd, "" if bed.itemRgb is None else bed.itemRgb, "" if bed.blockCount is None else bed.blockCount, "" if bed.blockSizes is None else ",".join(map(str, bed.blockSizes)), "" if bed.blockStarts is None else ",".join(map(str, bed.blockStarts))))


_bed_additional_fields = ["name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"]
def _block_func(s):
	return list(map(int, filter(lambda a:a, s.split(","))))
_bed_additional_field_funcs = [str, float, str, int, int, str, int, _block_func, _block_func] 		
class BEDXReader(BaseReader):
	'''
	A general BED-X class for customized BED file that extends the BEDx+ format.
	By default, input file is in bed3+ format. However, one could also apply this on bed6+, bed9+, etc.
	'''
	__slots__ = "fieldnames", "fieldfuncs", "BEDX"
	def __init__(self, arg, fieldnames=None, fieldfuncs=None, x=3, classname="BEDX", strandfield=None):
		'''
		fieldnames: List of names of the additional fields. 
		fieldfuncs: It can either be a list or a dict. 
		x: bedx+ 
		'''
		super(BEDXReader, self).__init__(arg)
		BEDXReaderSelf = self
		def _init(self, chrom, chromStart, chromEnd, *args):
			super(BEDXReaderSelf.BEDX, self).__init__(chrom, chromStart, chromEnd)
			if len(args) != len(self.__slots__):
				raise Exception("Inconsistent field names and entries")
			for slot, arg in zip(self.__slots__, args):
				setattr(self, slot, arg)
		if x < 3 or x > 12:
			raise Exception("Incorrect X")
		if fieldnames is None:
			fieldnames = [f"f{x + i}" for i in range(len(self._line.split("\t")) - x)]
		self.fieldnames = _bed_additional_fields[:x - 3] + fieldnames
		addition_field_funcs = {_bed_additional_fields[i]:_bed_additional_field_funcs[i] for i in range(x - 3)}
		if fieldfuncs is None:
			self.fieldfuncs = addition_field_funcs
		elif type(fieldfuncs) is dict:
			self.fieldfuncs = fieldfuncs
		else:
			self.fieldfuncs = _bed_additional_field_funcs[:x - 3] + fieldfuncs
		if strandfield is None:
			self.BEDX = type(classname, (BED3,), {"__slots__":self.fieldnames, "__init__":_init})
		else:
			@property
			def _stranded_genomic_pos(self, strandfield=strandfield):
				return StrandedGenomicPos(self.chrom, self.chromStart + 1, self.chromEnd, getattr(self, strandfield))
			self.BEDX = type(classname, (BED3, StrandedGenomicAnnotation), {"__slots__":self.fieldnames, "__init__":_init, "stranded_genomic_pos":_stranded_genomic_pos})
		
	def _read(self):
		line = self._line
		if line is None:
			return None
		self._proceed_next_line()
		words_array = line.split('\t')
		chrom = words_array[0]
		chromStart = int(words_array[1]) 
		chromEnd = int(words_array[2])
		if len(self.fieldnames) != len(words_array)-3:
			raise Exception("Number of fields in BED file does not match the number of field names")
		if self.fieldfuncs is None:
			return self.BEDX(chrom, chromStart, chromEnd, *words_array[3:])
		else:
			if type(self.fieldfuncs) is dict:
				args = [self.fieldfuncs[fieldname](obj) if fieldname in self.fieldfuncs else obj for obj, fieldname in zip(words_array[3:], self.fieldnames)]
			else:
				args = [fieldfunc(obj) for obj, fieldfunc in zip(words_array[3:], self.fieldfuncs)] 
			return self.BEDX(chrom, chromStart, chromEnd, *args)
class BEDXWriter(BaseWriter):
	def __init__(self, arg, fieldfuncs={}):
		super(BEDXWriter, self).__init__(arg)
		self.fieldfuncs = fieldfuncs
	def write(self, bed):
		'''
		Output the BEDX
		'''
		super(BEDXWriter, self).write(
			("\t".join(["{}"]*(3+len(bed.__slots__))) + "\n").format(
				bed.chrom, bed.chromStart, bed.chromEnd,
				*[(self.fieldfuncs[i] if i in self.fieldfuncs else str)(bed.__getattribute__(i)) for i in bed.__slots__]
			)
		)


class BEDGraph(BED3):
	__slots__ = "dataValue",
	def __init__(self, chrom, chromStart, chromEnd, dataValue):
		super(BEDGraph, self).__init__(chrom, chromStart, chromEnd)
		self.dataValue = dataValue
		
class BEDGraphReader(BaseReader):
	
	def __init__(self, arg, dataValueType=float):
		super(BEDGraphReader, self).__init__(arg)
		self.dataValueType = dataValueType
	def _metainfo_reader(self):
		self._proceed_next_line()
		while self._line is not None and (self._line.startswith("track") or self._line.startswith("browser")):
			self._proceed_next_line()
	
	def _read(self):
		line = self._line
		if line is None:
			return None
		self._proceed_next_line() 
		words_array = line.split('\t')
		chrom = words_array[0]
		chromStart = int(words_array[1]) 
		chromEnd = int(words_array[2])
		dataValue = self.dataValueType(words_array[3])
		return BEDGraph(chrom, chromStart, chromEnd, dataValue)
	
class BEDGraphWriter(BaseWriter):
	def __init__(self, arg, dataValueFunc=None):
		super(BEDGraphWriter, self).__init__(arg)
		if dataValueFunc is None:
			dataValueFunc = str
		self.dataValueFunc = dataValueFunc
		
	def write(self, bedgraph):
		'''
		Output the BEDGraph
		'''
		super(BEDGraphWriter, self).write("{}\t{}\t{}\t{}\n".format(bedgraph.chrom, bedgraph.chromStart, bedgraph.chromEnd, self.dataValueFunc(bedgraph.dataValue)))
	
	
class BEDGraphIReader(TabixIReader):
	def __init__(self, arg, tbi=None, dataValueType=float, missing_data_mode="zero"):
		super(BEDGraphIReader, self).__init__(arg, tbi)
		if isinstance(dataValueType, str): # Not very safe way to convert str into func
			dataValueType = eval(dataValueType, {})
		self.dataValueType = dataValueType
		if missing_data_mode == "zero":
			self.missing_data_value = 0
		elif missing_data_mode == "nan":
			self.missing_data_value = float("nan")
		
	def _parse_raw_entry(self, entry):
		return BEDGraph(entry[0], int(entry[1]), int(entry[2]), self.dataValueType(entry[3]))	
	def values_dict(self, r):
		r = GenomicPos(r)
		zstart = r.zstart
		ostop = r.ostop
		return {p+1 : bg.dataValue for bg in self.entries_iterator(r) for p in range(max(bg.genomic_pos.zstart, zstart), min(bg.genomic_pos.ostop, ostop))}		
	def values(self, r):
		r = GenomicPos(r)
		d = self.values_dict(r)
		return [d[i+1] if i+1 in d else self.missing_data_value for i in range(r.zstart, r.ostop)]

class BEDIReader(TabixIReader):
	def _is_empty_field(self, s):
		return s is None or s == "" or s == "."
	
	def _parse_raw_entry(self, entry):
		words_array = entry
		chrom = words_array[0]
		chromStart = int(words_array[1]) 
		chromEnd = int(words_array[2])
		words_array.extend([None] * (12 - len(words_array)))
		
		name = words_array[3]
		score = None
		strand = None
		thickStart = None
		thickEnd = None
		itemRgb = None
		blockCount = None
		blockSizes = None
		blockStarts = None
		
		if not self._is_empty_field(words_array[4]):
			score = float(words_array[4])
		if words_array[5] != "":
			strand = words_array[5]
		if not self._is_empty_field(words_array[6]):
			thickStart = words_array[6]
		if not self._is_empty_field(words_array[7]):
			thickEnd = words_array[7]
		if not self._is_empty_field(words_array[8]):
			itemRgb = words_array[8]
		if not self._is_empty_field(words_array[9]):
			blockCount = int(words_array[9])
		if not self._is_empty_field(words_array[10]):
			blockSizes = list(map(int, filter(lambda a:a, words_array[10].split(","))))
		if not self._is_empty_field(words_array[11]):
			blockStarts = list(map(int, filter(lambda a:a, words_array[11].split(","))))
		
		if blockSizes is not None:
			if blockCount != len(blockSizes):
				raise Exception("Inconsistent blockCount and blockSizes")
			if blockCount != len(blockStarts):
				raise Exception("Inconsistent blockCount and blockStarts")
			if blockStarts[0] != 0:
				raise Exception("First block must start at 0")
			if chromStart + blockStarts[-1] + blockSizes[-1] != chromEnd:
				raise Exception("Last block must end at chromEnd")		 
		
		return BED(chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts)
		
class BEDPE():
	'''
	'''
	def __init__(self, chrom1, start1, stop1, chrom2, start2, stop2, name=None, score=None, strand1=None, strand2=None):
		self.chrom1 = chrom1
		self.start1 = start1
		self.stop1 = stop1
		self.chrom2 = chrom2
		self.start2 = start2
		self.stop2 = stop2
		self.name = name
		self.score = score
		self.strand1 = strand1
		self.strand2 = strand2
	@property
	def genomic_pos1(self):
		return GenomicPos(self.chrom1, self.start1 + 1, self.stop1)
	@property
	def stranded_genomic_pos1(self):
		return StrandedGenomicPos(self.chrom1, self.start1 + 1, self.stop1, self.strand1)
	@property
	def genomic_pos2(self):
		return GenomicPos(self.chrom2, self.start2 + 1, self.stop2)
	@property
	def stranded_genomic_pos2(self):
		return StrandedGenomicPos(self.chrom2, self.start2 + 1, self.stop2, self.strand2)
	
class BEDPEReader(BaseReader):
	def __init__(self, arg):
		super(BEDPEReader, self).__init__(arg)
	def _proceed_next_line(self):
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
		if self._line is None:
			return None
		words_array = self._line.split('\t')
		self._proceed_next_line()
		name = None
		score = None
		strand1 = None
		strand2 = None
		chrom1 = words_array[0]
		start1 = int(words_array[1])
		stop1 = int(words_array[2])
		chrom2 = words_array[3]
		start2 = int(words_array[4])
		stop2 = int(words_array[5])
		if len(words_array) >= 7:
			name = words_array[6]
		if len(words_array) >= 8:
			score = words_array[7]
		if len(words_array) >= 9:
			strand1 = words_array[8]
		if len(words_array) >= 10:
			strand2 = words_array[9]
		return BEDPE(chrom1, start1, stop1, chrom2, start2, stop2, name, score, strand1, strand2)
		
	
	
class BEDPEWriter(BaseWriter):
	def __init__(self, arg):
		super(BEDPEWriter, self).__init__(arg)
			
	def write(self, bedpe):
		'''
		Output the BEDPENode
		Note that the start position is stored in 1-based but output in 0-based
		'''
		super(BEDPEWriter, self).write(
			"\t".join(list(map(str, [
				bedpe.chrom1, 
				bedpe.start1, 
				bedpe.stop1,
				bedpe.chrom2, 
				bedpe.start2, 
				bedpe.stop2, 
				"" if bedpe.name is None else bedpe.name, 
				"" if bedpe.score is None else bedpe.score, 
				"" if bedpe.strand1 is None else bedpe.strand1,
				"" if bedpe.strand2 is None else bedpe.strand2
				]))) + "\n")


class PINTSBidirectional(BED3):
	__slots__ = "confidence", "major_tss_pls", "major_tss_mns"
	def __init__(self, chrom, chromStart, chromEnd, confidence, major_tss_pls, major_tss_mns):
		super(PINTSBidirectional, self).__init__(chrom, chromStart, chromEnd)
		self.confidence = confidence
		self.major_tss_pls = major_tss_pls
		self.major_tss_mns = major_tss_mns
class PINTSBidirectionalReader(BaseReader):
	def __init__(self, arg):
		super(PINTSBidirectionalReader, self).__init__(arg)
	
	def _read(self):
		line = self._line
		if line is None:
			return None
		self._proceed_next_line() 
		words_array = line.split('\t')
		chrom = words_array[0]
		chromStart = int(words_array[1]) 
		chromEnd = int(words_array[2])
		confidence = set(words_array[3].split(","))
		major_tss_pls = list(map(int, words_array[4].split(",")))
		major_tss_mns = list(map(int, words_array[5].split(",")))
		return PINTSBidirectional(chrom, chromStart, chromEnd, confidence, major_tss_pls, major_tss_mns)

class BigBedIReader(BaseIReader):
	__slots__ = "bigbed", "fieldnames", "fieldfuncs", "BEDX"
	def __init__(self, arg, fieldnames=[], fieldfuncs=None, x=3, classname="BEDX"):
		'''
		fieldnames: List of names of the additional fields. 
		fieldfuncs: It can either be a list or a dict. 
		x: bedx+ 
		'''
		super(BigBedIReader, self).__init__(arg)
		if not _SUPPORT_BIGWIG:
			raise Exception("bigbed is not supported without pyBigWig.")
		BigBedIReaderSelf = self
		def _init(self, chrom, chromStart, chromEnd, *args):
			super(BigBedIReaderSelf.BEDX, self).__init__(chrom, chromStart, chromEnd)
			if len(args) != len(self.__slots__):
				raise Exception("Inconsistent field names and entries")
			for slot, arg in zip(self.__slots__, args):
				setattr(self, slot, arg)
		if x < 3 or x > 12:
			raise Exception("Incorrect X")
		self.fieldnames = _bed_additional_fields[:x - 3] + fieldnames

		addition_field_funcs = {_bed_additional_fields[i]:_bed_additional_field_funcs[i] for i in range(x - 3)}
		if fieldfuncs is None:
			self.fieldfuncs = addition_field_funcs
		elif type(fieldfuncs) is dict:
			self.fieldfuncs = fieldfuncs
		else:
			self.fieldfuncs = _bed_additional_field_funcs[:x - 3] + fieldfuncs
		self.BEDX = type(classname, (BED3,), {"__slots__":self.fieldnames, "__init__":_init})
		self.bigbed = pyBigWig.open(arg)

	def entries_iterator(self, r):
		r = GenomicPos(r)
		if r.name in self.bigbed.chroms():
			entries = self.bigbed.entries(r.name, r.zstart, r.ostop)
			if entries is not None:
				for zstart, ostop, s in entries:
					yield self.BEDX(r.name, zstart, ostop, *s.split())
					
	def entries(self, r):
		return list(self.entries_iterator(r))
		
	def __getitem__(self, key):
		return self.entries[key]
	
	def close(self):
		self.bigbed.close()

class ENCODENarrowPeak(BED3, StrandedGenomicAnnotation):
	def __init__(self, chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak):
		self.chrom = chrom
		self.chromStart = chromStart
		self.chromEnd = chromEnd
		self.name = name
		self.score = score
		self.strand = strand
		self.signalValue = signalValue
		self.pValue = pValue
		self.qValue = qValue
		self.peak = peak
		
	@property
	def stranded_genomic_pos(self):
		return StrandedGenomicPos(self.chrom, self.chromStart + 1, self.chromEnd, self.strand)

class ENCODENarrowPeakReader(BED3Reader):
	def _is_empty_field(self, s):
		return s is None or s == "" or s == "."	
	def _read(self):
		line = self._line
		if line is None:
			return None
		self._proceed_next_line() 
		words_array = line.split('\t')
		chrom = words_array[0]
		chromStart = int(words_array[1]) 
		chromEnd = int(words_array[2])
		name = words_array[3]
		score = float(words_array[4])
		strand = words_array[5]
		signalValue = float(words_array[6])
		pValue = float(words_array[7])
		qValue = float(words_array[8])
		peak = int(words_array[9])
		return ENCODENarrowPeak(chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak)
	
