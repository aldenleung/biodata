# biodata - A standard biological data processing package

The biodata package provides a standard API to access all different kinds of biological data using similar syntax. For each data type, data is processed by the corresponding reader (XX-Reader) and writer (XX-Writer) as a stream of entries. For example, `FASTAReader` is used to process FASTA file. A call to the method `read()` from `FASTAReader` yields a `FASTA` object. For indexed file, random access is supported through the XX-IReader. For example, an indexed FASTA file can be access by  `FASTAIReader`. 

## Basic usage

We will demonstrate the use of biodata package using FASTA file. 

```
>seq1
ACGT
>seq2
CCCGGGAAA
```

### Read the first entry

```python
from biodata.fasta import FASTAReader
with FASTAReader(filename) as fr:
	f = fr.read()
	print(f.name, f.seq) # seq1 ACGT
```

### Read entry by entry

```python
from biodata.fasta import FASTAReader
with FASTAReader(filename) as fr:
	for f in fr:
		print(f.name, f.seq)
# seq1 ACGT
# seq2 CCCGGGAAA
```

### Read all entries at once

```python
from biodata.fasta import FASTAReader
fasta_entries = FASTAReader.read_all(list, filename) # list of FASTA

seq_dict = FASTAReader.read_all(lambda fr: {f.name:f.seq for f in fr}, filename) 
# A dictionary with fasta name as key and fasta sequence as value
# {"seq1": "ACGT", "seq2": "CCCGGGAAA"}
```

### Peek an entry

```python
from biodata.fasta import FASTAReader
with FASTAReader(filename) as fr:
	f = fr.peek() # Only peek the entry without proceeding to the next entry
	print(f.name, f.seq) # seq1 ACGT
	f = fr.read() # Read the entry and proceed to the next entry
	print(f.name, f.seq) # seq1 ACGT
	f = fr.read()
	print(f.name, f.seq) # seq2 CCCGGGAAA
```

### Read an entry from StringIO

```python
# TextIOBase can be used as input
import io
from biodata.fasta import FASTAReader
FASTAReader.read_all(list, io.StringIO(">seq1\nACGT\n>seq2\nCCCGGGAAA\n"))
```

### Read an indexed file

```python
from biodata.fasta import FASTAIReader
from genomictools import GenomicPos, StrandedGenomicPos
from biodata.bed import BED

fir = FASTAIReader(filename, faifilename) # fai file can be created using 'samtools faidx filename'
f = fir[GenomicPos("seq2:1-4")] # Read from a region without strand
print(f.name, f.seq) # seq2:1-4 CCCG
f = fir[StrandedGenomicPos("seq2:1-4:-")] # Read from a region with strand
print(f.name, f.seq) # seq2:1-4:- CGGG
f = fir[BED("seq2", 0, 4, strand="-")] # Equivalent to StrandedGenomicPos but a BED entry is used
print(f.name, f.seq) # seq2:1-4:- CGGG
fir.close()
```

### Write entry by entry

```python
from biodata.fasta import FASTA, FASTAWriter
with FASTAWriter(output_file) as fw:
	fw.write(FASTA("seq1", "ACGT"))
	fw.write(FASTA("seq2", "CCCGGGAAA"))
```

### Write all entries at once

```python
from biodata.fasta import FASTA, FASTAWriter
fasta_entries = [FASTA("seq1", "ACGT"), fw.write(FASTA("seq2", "CCCGGGAAA"))]
FASTAWriter.write_all(fasta_entries, output_file)
```



# List of supported format

1. Delimited (`biodata.delimited`)
2. FASTA, FASTQ (`biodata.fasta`)
3. BED3, BED, BEDX, BEDGraph, BEDPE (`biodata.bed`)

Future supported formats. 

1. GFF (`biodata.gff`)
2. VCF (`biodata.vcf`)
3. BigBed (`biodata.bed`)
4. BigWig (`biodata.bigwig`)



## Extension of BaseReader

Users can extend the BaseReader and BaseWriter class easily.

```python
class ExampleNode(object):
	def __init__(self, value1, value2):
		self.value1 = value1
		self.value2 = value2

class ExampleNodeReader(BaseReader):
	def __init__(self, filename):
		super(ExampleNodeReader, self).__init__(filename)
	def _read(self):
		if self.line is None:
			return None
		words_array = self.line.split('\t')
		value1 = words_array[0]
		value2 = words_array[1]
		self.proceed_next_line()
		return ExampleNode(value1, value2)

filename = "SomeDocument.txt"
with ExampleNodeReader(filename) as er:
	for node in er:
		print(node.value1, node.value2)
```


