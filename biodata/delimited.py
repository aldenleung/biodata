'''
Created on Mar 31, 2021

@author: Alden
'''
import re

from .baseio import BaseReader, BaseWriter, get_text_file_extension
from collections import OrderedDict

def _guess_separator(f):
	# Auto separator detection based on file extension
	# This may be expanded in the future 
	file_extension = get_text_file_extension(f)
	if (file_extension == "tsv" or file_extension == "txt"): # Assume tab-delimited
		return "\t"
	elif (file_extension == "csv"): # Assume comma-delimited
		return ","
	else: 
		return "\t" # Default as tab-delimited

class DelimitedReader(BaseReader):
	'''
	A simple class for delimited file processing. Each line should have the same number of columns.
	
	This class targets for reading entry by entry (i.e. lines) instead of an entire table. 
	
	If the table has both column names and row names, and the header line contains exactly n_field_from_other_line - 1,
	 
	 
	One can define header by (1) setting header=True to automatically use the first non-comment line as header, or (2) using custom_header to define the header names. 
	The result object will be wrapped into a namedtuple  
	
	One can convert the element of each cell to specific type by defining custom_funcs.
	
	Example usage:
	
	.. code-block:: python
		
		import io
		from biodata.delimited import DelimitedReader
		
		print("===== Remove comment line =====")
		data = '# i am comment\n1\t2\n3\t4'
		print(data)
		print("-----------------------")
		with DelimitedReader(io.StringIO(data)) as dr:
			print(list(dr))
		with DelimitedReader(io.StringIO(data), comment_symbol="#", funcs=[int, int]) as dr:
			print(list(dr))
		print()
		print("===== Header =====")
		data = 'h1\th2\n1\t2\n3\t4'
		print(data)
		print("-----------------------")
		with DelimitedReader(io.StringIO(data)) as dr:
			print(list(dr))
		with DelimitedReader(io.StringIO(data), header=True) as dr:
			print(list(dr))
		with DelimitedReader(io.StringIO(data), header=True, funcs={"h1":int, "h2":int}) as dr:
			print(list(dr))
		print()
		print("===== header with comment symbol =====")
		data = '#h1\th2\n1\t2\n3\t4'
		print(data)
		print("-----------------------")
		with DelimitedReader(io.StringIO(data), header=True) as dr:
			print(list(dr))
		with DelimitedReader(io.StringIO(data), header=True, skip_header_comment_symbol="#") as dr:
			print(list(dr))
		print()
		print("===== Row name without header =====")
		data = 'r1\t1\t2\nr2\t3\t4\n'
		print(data)
		print("-----------------------")
		with DelimitedReader(io.StringIO(data)) as dr:
			print(list(dr))
		with DelimitedReader(io.StringIO(data), funcs=[str, int, int]) as dr:
			print(list(dr))
		with DelimitedReader(io.StringIO(data), funcs={1:int, 2:float}) as dr:
			print(list(dr))
		print()
		print("===== Row name with header =====")
		data = 'h1\th2\nr1\t1\t2\nr2\t3\t4'
		print(data)
		print("-----------------------")
		with DelimitedReader(io.StringIO(data)) as dr:
			print(list(dr))
		with DelimitedReader(io.StringIO(data), header=True, funcs={"h1":int, "h2":int}, additional_row_header="row_header") as dr:
			print(list(dr))
		print()
		
	.. code-block:: 
	
		===== Remove comment line =====
		# i am comment
		1	2
		3	4
		-----------------------
		[['# i am comment'], ['1', '2'], ['3', '4']]
		[[1, 2], [3, 4]]
		
		===== Header =====
		h1	h2
		1	2
		3	4
		-----------------------
		[['h1', 'h2'], ['1', '2'], ['3', '4']]
		[OrderedDict([('h1', '1'), ('h2', '2')]), OrderedDict([('h1', '3'), ('h2', '4')])]
		[OrderedDict([('h1', 1), ('h2', 2)]), OrderedDict([('h1', 3), ('h2', 4)])]
		
		===== header with comment symbol =====
		#h1	h2
		1	2
		3	4
		-----------------------
		[OrderedDict([('#h1', '1'), ('h2', '2')]), OrderedDict([('#h1', '3'), ('h2', '4')])]
		[OrderedDict([('h1', '1'), ('h2', '2')]), OrderedDict([('h1', '3'), ('h2', '4')])]
		
		===== Row name without header =====
		r1	1	2
		r2	3	4
		
		-----------------------
		[['r1', '1', '2'], ['r2', '3', '4']]
		[['r1', 1, 2], ['r2', 3, 4]]
		[['r1', 1, 2.0], ['r2', 3, 4.0]]
		
		===== Row name with header =====
		h1	h2
		r1	1	2
		r2	3	4
		-----------------------
		[['h1', 'h2'], ['r1', '1', '2'], ['r2', '3', '4']]
		[OrderedDict([('row_header', 'r1'), ('h1', 1), ('h2', 2)]), OrderedDict([('row_header', 'r2'), ('h1', 3), ('h2', 4)])]

	'''
	#__slots__ = "separator", "AutoDelimitedClass", "header_keys", "custom_funcs"
	
	def __init__(self, arg, *, separator=None, header=None, custom_header=None, funcs=None, comment_symbol=None, quote=None, skip_lines=None, skip_header_comment_symbol=None, additional_row_header=None):
		self.initialized = False
		super(DelimitedReader, self).__init__(arg)
		if isinstance(arg, str): 
			if separator is None:
				# Auto separator detection based on file extension
				self.separator = _guess_separator(arg)
			else:
				self.separator = separator
		else:
			if separator is None:
				self.separator = "\t" # By default use \t
			else:
				self.separator = separator
		
		# Raise exception for incorrect parameter settings	
		if header is not None and custom_header is not None:
			raise Exception()
		if header is None and skip_header_comment_symbol is not None:
			raise Exception()
		if header is None and additional_row_header is not None:
			raise Exception()
		
		self.header = header
		self.custom_header = custom_header
		self.funcs = funcs
		self.comment_symbol = comment_symbol
		self.quote = quote
		if quote is None:
			self.split_func = lambda line: line.split(self.separator)
		else: 
			self.split_func = lambda line: list(map("".join, re.findall(f'(?:(?<={self.separator})|(?<=^))(?:([^{self.quote}{self.separator}]*)|{self.quote}([^{self.quote}]*){self.quote})(?:(?={self.separator})|(?=$))', line)))
		#self.split_func = lambda line: list(csv.reader([line], delimiter=self.separator))[0]
		self.skip_lines = skip_lines
		self.skip_header_comment_symbol = skip_header_comment_symbol
		self.additional_row_header = additional_row_header
		self.initialized = True
		self._metainfo_reader()
		
			
	def _metainfo_reader(self):
		if not self.initialized:
			return
		if self.skip_lines is not None:
			for _ in range(self.skip_lines):
				self.f.readline()
		
		self._proceed_next_line()
		# Try to parse header
		self.header_keys = None
		if self.header is not None:
			line = self._line
			if self.skip_header_comment_symbol is not None:
				if not line.startswith(self.skip_header_comment_symbol):
					raise Exception()
				line = line[len(self.skip_header_comment_symbol):]
			words_array = self.split_func(line)
			if len(set(words_array)) != len(words_array):
				raise Exception()
			if self.additional_row_header is not None:
				self.header_keys = [self.additional_row_header] + words_array
			else:
				self.header_keys = words_array
			self._proceed_next_line()
			
		if self.custom_header is not None:
			self.header_keys = self.custom_header
		
			
		
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
			if line != '' and (self.comment_symbol is None or not line.startswith(self.comment_symbol)): # Auto comment removal. Blank lines are skipped by default
				self._line = line
				break
		
	def _read(self):
		if self._line is None:
			return None
		
		words_array = self.split_func(self._line)
		self._proceed_next_line()
		
		if self.header_keys is None:
			if self.funcs is not None:
				if type(self.funcs) is dict:
					words_array = [self.funcs[key](value) if key in self.funcs else value for key, value in enumerate(words_array)]
				else:
					if len(self.funcs) != len(words_array):
						raise Exception()
					words_array = [f(v) for v, f in zip(words_array, self.funcs)]
			return words_array
		else:
			if len(self.header_keys) != len(words_array):
				print(self.header_keys)
				print(words_array)
				raise Exception()
			if self.funcs is not None:
				words_array = [self.funcs[key](value) if key in self.funcs else value for key, value in zip(self.header_keys, words_array)]
			return OrderedDict((key, value) for key, value in zip(self.header_keys, words_array))


class DelimitedWriter(BaseWriter):
	'''
	A simple class for delimited file output.
	
	Example usage:
	
	.. code-block:: python
	
		with DelimitedWriter("Output.txt") as dw:
			dw.write(["A1", "A2", "A3"])
			dw.write(["B1", "B2", "B3"])
	'''
	__slots__ = "separator"
	def __init__(self, arg, separator=None):
		super(DelimitedWriter, self).__init__(arg)
		if separator is None:
			if isinstance(arg, str):
				self.separator = _guess_separator(arg)
			else:
				self.separator="\t" # Use \t
		else:
			self.separator = separator
			
	def write(self, arr):
		'''
		Output any iterable
		'''
		super(DelimitedWriter, self).write(self.separator.join(str(d) for d in arr) + "\n")
	
	
	
	