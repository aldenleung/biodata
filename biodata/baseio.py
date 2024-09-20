import sys
import bz2
import gzip
from io import TextIOBase, StringIO
try:
	import Bio.bgzf
	_SUPPORT_BGZ = True
except:
	_SUPPORT_BGZ = False
class BaseReader(object):
	'''
	The base class for all sequential Readers handling different text file types. By default it (1) reads line-by-line, (2) skips blank lines, (3) rstrips the line.
	
	Any subclass should override the _read function. Other functions such as _proceed_next_line and _metainfo_reader should be overrided if applicable. 
	
	The file name '-' is treated as stdin. If a file ends in gz / bz2, it will use gzip / bz2 to open the file with utf-8 encoding

	Example usage:

	.. code-block:: python
		
		# The following are equivalent in reading a files into lines.
		from biodata.baseio import BaseReader
		with BaseReader(filename) as br:
			lines = br.read_all()
		
		from biodata.baseio import BaseReader
		lines = []
		with BaseReader(filename) as br:
			for b in br:
				lines.append(b)
		
		from biodata.baseio import BaseReader
		lines = BaseReader.read_all(list, filename)

		# TextIOBase can be used as input
		import io
		from biodata.baseio import BaseReader
		with BaseReader(io.StringIO("Line1\nLine2\n")) as br:
			lines = br.read_all()
		
	.. code-block:: python
	
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
				words_array = self.line.split('\\t')
				value1 = words_array[0]
				value2 = words_array[1]
				self.proceed_next_line()
				return ExampleNode(value1, value2)
		
		filename = "SomeDocument.txt"
		with ExampleNodeReader(filename) as er:
			node = er.read()
			while node is not None:
				print(node.value1 + "\\t" + node.value2)
				node = er.read()
	
	
	'''
	
	__slots__ = "f", "_line", "_is_stored", "_stored_value"
	def __init__(self, arg):
		'''
		Either an arg or 
		'''
		if isinstance(arg, TextIOBase):
			self.f = arg 
		else:
			self.f = create_text_stream(arg, "r")
		self._is_stored = False
		self._metainfo_reader()
			
		
	def _proceed_next_line(self):
		'''
		Attempts to assign the next line to self._line. 
		It is recommended for any subclass to override and use this method on how to proceed to next line (e.g. Ignore comment line starting with #) 	
		'''
		while True:
			line = self.f.readline()
			if line == '':
				self._line = None
				break
			line = line.rstrip("\r\n") # Auto-stripping
			if line != '': # Auto comment removal. Blank lines are skipped by default
				self._line = line
				break
			
	
	def _metainfo_reader(self):
		'''
		Attempts to process meta info in the data file. 
		Subclasses can override this method to skip headers with other formats or process the comment / header lines (e.g. file format version)
		'''
		self._proceed_next_line()
	
	
	def _read(self):
		line = self._line
		self._proceed_next_line()
		return line
	
	def read(self):
		'''
		Returns the next object
		'''
		if self._is_stored:
			self._is_stored = False
			return self._stored_value
		else:
			return self._read()
	
	def peek(self):
		if not self._is_stored:
			self._is_stored = True 
			self._stored_value = self._read()
		return self._stored_value 
	
	def __iter__(self):
		'''
		Create a generator for reading data
		'''
		obj = self.read()
		while obj is not None:
			yield obj
			obj = self.read()
			
	def close(self):
		'''
		Close the file
		'''
		if self.f is not sys.stdin:
			self.f.close()
	
	def __enter__(self): 
		return self
	
	def __exit__(self, type, value, traceback): 
		self.close()
	
	@classmethod
	def read_all(cls, read_all_func, *args, **kwargs):
		with cls(*args, **kwargs) as br:
			return read_all_func(br)

class BaseWriter(object):
	'''
	The base text class for all Writers handling different text file types. By default it overwrites the target file.
	
	The file name '-' is treated as stdout. If a file ends in gz / bz2, it will use gzip / bz2 to open the file with utf-8 encoding
	
	Example usage:

	.. code-block:: python
	
		with BaseWriter(filename) as bw:
			bw.write("This is a sentence written to the specified file.")
		
	.. code-block:: python
		
		with BaseWriter(create_file_stream(filename, mode="a")) as bw:
			bw.write("This is a sentence appended to the specified file")
		
	.. code-block:: python
	
		class ExampleNode(object):
			def __init__(self, value1, value2):
				self.value1 = value1
				self.value2 = value2
		
		class ExampleNodeWriter(BaseWriter):
			def __init__(self, filename):
				super(ExampleNodeWriter, self).__init__(filename)
			def _initialize_header(self):
				super(ExampleNodeWriter, self).write("#Value1\\tValue2\\n")
			def write(self, node):
				super(ExampleNodeWriter, self).write(node.value1 + "\\t" + node.value2 + "\\n")
		
		node = ExampleNode("1", "This is the sentence to be written in the specified file.")
		filename = "SomeDocument.txt"
		with ExampleNodeWriter(filename) as ew:
			ew.write(node)
		
	'''
	
	__slots__ = "f"
	def __init__(self, arg, mode="w"):
		if isinstance(arg, TextIOBase):
			self.f = arg
		else:
			self.f = create_text_stream(arg, mode)
		self._initialize_header()
	
	def _initialize_header(self):
		'''
		This method is automatically called during initialization. Override this method to initialize and output the headers before writing the content.
		'''
		pass
	
	def write(self, s):
		'''
		Writes the string to the file
		'''
		self.f.write(s)
	
	def close(self):
		'''
		Closes the file
		'''
		if self.f is not sys.stdout and not isinstance(self.f, StringIO):
			self.f.close()
		
	def __enter__(self): 
		return self
	
	def __exit__(self, type, value, traceback): 
		self.close()

	@classmethod
	def write_all(cls, iterable, arg, *args, **kwargs):
		with cls(arg, *args, **kwargs) as bw:
			for i in iterable:
				bw.write(i)

class NullWriter(object):
	'''
	Creates a NULL writer, that does nothing when the user called the method write. This writer class is useful for optional output file. 
	

	Example usage:

	.. code-block:: python
	
		filename = None
		if filename is not None:
			bw = BaseWriter(filename)
		else:
			bw = NullWriter()
		bw.write("If filename is provided, this sentence will be written to the file. Otherwise, it will be written to nowhere")
		bw.close()
		
	'''
	def __init__(self):
		pass

	def write(self, s):
		pass
	
	def write_all(self, elements):
		pass
			
	def close(self):
		pass
		
	def __enter__(self): 
		return self
	
	def __exit__(self, type, value, traceback): 
		self.close()


def get_text_file_extension(filename, neglect_zip=True):
	'''
	Neglect the zip extension and return file extension 'txt', 'csv', etc.
	
	Only neglect one level of zip extension. "file.gz.bz2" will return 'gz'
	
	:param filename: The file name
	
	:returns: The file extension after decompression
	
	Example usage:
	
	.. code-block:: python
	
		get_text_file_extension("a.txt") # 'txt'
		
		get_text_file_extension("a.txt.gz") # 'txt'
		
		get_text_file_extension("a.csv.bz2") # 'csv'
		
		get_text_file_extension("a.txt.gz", neglect_zip = False) # 'gz'
		
		get_text_file_extension("example") # ''
		
		get_text_file_extension("a.csv.txt") # 'txt'
	''' 
	supported_zip_extensions = [".gz", ".bz2", ".bgz"]
	if (neglect_zip):
		for supported_zip_extension in supported_zip_extensions:
			if filename.endswith(supported_zip_extension):
				filename = filename[0:filename.rfind(supported_zip_extension)]
				break
	
	index = filename.rfind(".")
	if index == -1: # No file extension found
		return ""
	else:
		return filename[index+1:] 
	
def get_text_file_rootname(filename, neglect_zip=True):
	supported_zip_extensions = [".gz", ".bz2", ".bgz"]
	if (neglect_zip):
		for supported_zip_extension in supported_zip_extensions:
			if filename.endswith(supported_zip_extension):
				filename = filename[0:filename.rfind(supported_zip_extension)]
				break
	index = filename.rfind(".")
	if index == -1: # No file extension found
		return filename
	else:
		return filename[:index] 
	

	
def create_text_stream(filename, mode):
	'''
	Create a text file stream of the file.
	
	Other than a file name, filename can also be:
	 
	(1) The filename "-" is interpreted as stdin or stdout.
	(2) The filename can be urllib.request.Request. 
	
	A decompressed stream is returned for compressed file (.gz or .bz2).
	
	All encodings are assumed utf-8    
	
	:param filename: The file name
	
	:returns: a TextIOBase representing the file stream
	
	Example usage:
	
	.. code-block:: python
		
		create_file_stream("file.txt", mode="r")
		create_file_stream(urllib.request.Request(url), mode="r")
	
	'''
	import urllib.request
	if isinstance(filename, urllib.request.Request):
		import io
		request = filename
		result = urllib.request.urlopen(request)
		remote_file_name = result.info().get_filename()
		stream = io.BytesIO(result.read())
		if remote_file_name.endswith(".bgz"):
			f = gzip.open(stream, mode = mode + 't', encoding='utf-8')
		elif remote_file_name.endswith(".gz"):
			f = gzip.open(stream, mode = mode + 't', encoding='utf-8')		
		elif remote_file_name.endswith(".bz2"):
			f = bz2.open(stream, mode = mode + 't', encoding='utf-8')
		else:
			f = io.TextIOWrapper(stream, encoding='utf-8')
	else:
		if filename == "-":
			if mode == 'r':
				f = sys.stdin
			else:
				f = sys.stdout
		else:
			if filename.endswith(".bgz"):
				if mode == "r":
					f = gzip.open(filename, mode = mode + 't', encoding='utf-8')
				else:
					if not _SUPPORT_BGZ:
						raise Exception("bgzip is not supported without Bio.bgzf")
					f = Bio.bgzf.open(filename, mode=mode+'t')
						
			elif filename.endswith(".gz"):
				f = gzip.open(filename, mode = mode + 't', encoding='utf-8')			
			elif filename.endswith(".bz2"):
				f = bz2.open(filename, mode = mode + 't', encoding='utf-8')
			else:
				f = open(filename, mode)
	return f


class BaseIReader(object):
	'''
	The base class for all readers of indexed-files. 
	For a typical file types that allows quick query onto certain entries.
	There is no guaranteed on what are indexed. 
	Also, there may not be any unique key per entry.
	Some files may have an additional associated index file.
	
	fasta, region-related, etc.  
	'''
	def __init__(self, arg):
		'''
		'''
		pass

	def __getitem__(self, key):
		pass
	
	def __enter__(self): 
		return self
	
	def __exit__(self, type, value, traceback): 
		self.close()
