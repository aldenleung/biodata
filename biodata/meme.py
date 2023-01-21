from collections import OrderedDict
import re
import shlex
from .baseio import BaseReader

class MEMEAlphabetSymbol:
	def __init__(self, character, name, color):
		self.character = character
		self.name = name
		self.color = color
	
class MEMEAlphabet():
	def __init__(self, core_symbols, name=None, standard_like=None,
				ambig_symbols_dict={}, complements={}):
		self.core_symbols = core_symbols
		self.name = name
		self.standard_like = standard_like
		self.ambig_symbols_dict = ambig_symbols_dict
		self.complements = complements
	
	def __contains__(self, s):
		return s in self.core_symbols or s in self.ambig_symbols_dict
	
	def __len__(self):
		return len(self.core_symbols)
	
def _check_valid_alphabet_symbol(c):
	return re.match(r"^[A-Za-z0-9.\-*?]$", c) is not None
		
	
def _parse_alphabet(lines):
	'''
	Parse alphabet lines and return MEMEAlphabet
	'''
	def subparse(s):
		items = shlex.split(s)
		if len(items) > 1:
			character, name, color = items
			name = name[1:-1]
		else:
			character, name, color = items[0], None, None
		return MEMEAlphabetSymbol(character, name, color)
	
	name = None
	standard_like = None
	core_symbols = []
	complements = {}
	ambig_symbols_dict = {}
	for line in lines:
		line = line.strip()
		if line.startswith("#"):
			continue
		line = line.split("#")[0].strip()
		if len(line) == 0:
			continue
		if line.startswith("ALPHABET"):
			if "=" in line:
				for c in line.split("=")[-1].strip():
					core_symbols.append(MEMEAlphabetSymbol(c, None, None))
			else:
				name, standard_like = shlex.split(line)[1:]
				name = name[1:-1]
		elif "~" in line: # Complements ymbols
			tmp_entries = []
			for entry in line.split("~"):
				tmp_entries.append(subparse(entry.strip()))
			e1, e2 = tmp_entries
			complements[e1.character] = e2.character
			complements[e2.character] = e1.character
			core_symbols.extend(tmp_entries)
		elif "=" in line: # Ambiguous symbols
			k, v = line.split("=")
			ambig_symbols_dict[k.strip()] = v.strip()
		else:
			# Core symbols
			core_symbols.append(subparse(line.strip()))
	for symbol in core_symbols:
		if not _check_valid_alphabet_symbol(symbol.character):
			raise Exception(f"Invalid core symbol character - {symbol.character}")		
	for c in ambig_symbols_dict:
		if not _check_valid_alphabet_symbol(c):
			raise Exception("Invalid ambiguous symbol character")
	for characters in ambig_symbols_dict.values():
		for c in characters:
			if not any(c == symbol.character for symbol in core_symbols):
				raise Exception(f"Ambiguous symbol character matched to unknown core symbols - {c}")
	return MEMEAlphabet(core_symbols, name, standard_like,
				ambig_symbols_dict, complements)

	
	
	
class MEMEMotif():
	def __init__(self, identifier, alternate_name, matrix_dict, pwm, url):
		self.identifier = identifier
		self.alternate_name = alternate_name
		self.matrix_dict = matrix_dict
		self.pwm = pwm
		self.url = url
		
	@property	
	def motif_length(self):
		return len(self.pwm)
	@property
	def alphabet_length(self):
		return len(self.pwm[0])
	
class MEMEMotifReader(BaseReader):
	'''
	'''
	def _metainfo_reader(self):
		self._proceed_next_line()
		line = self._line
		meme_version = None
		alphabet = None
		strands = None
		background_letter_frequencies = None
		warnings = []
		while not line.startswith("MOTIF"):
			if line.startswith("MEME version"):
				meme_version = line.split("MEME version")[-1].strip()
				self._proceed_next_line()
			elif line.startswith("ALPHABET="):
				alphabet = _parse_alphabet([line])
				self._proceed_next_line()
			elif line.startswith("ALPHABET"):
				alphabet_lines = []
				while line != "END ALPHABET":
					alphabet_lines.append(line)
					self._proceed_next_line()
					line = self._line
				alphabet = _parse_alphabet(alphabet_lines)
				self._proceed_next_line()
			elif line.startswith("strands:"):
				strands = line.split("strands:")[-1].strip().split()
				self._proceed_next_line()
			elif line.startswith("Background letter frequencies"):
				pattern = re.compile(r"([^\s])\s(0\.[0-9]{3})")
				background_letter_frequencies = {}
				self._proceed_next_line()
				line = self._line
				while not line.startswith("MOTIF"):
					background_letter_frequencies.update({k:float(v) for k, v in pattern.findall(line)})
					self._proceed_next_line()
					line = self._line
				if sum(background_letter_frequencies.values()) != 1:
					warnings.append("Background letter frequencies not sum to 1")
				if alphabet is not None and len(background_letter_frequencies) != len(alphabet):
					warnings.append("Unmatched length of alphabet and background letter frequencies")
			else:
				self._proceed_next_line()
			line = self._line
		self.meme_version = meme_version
		self.alphabet = alphabet
		self.strands = strands
		self.background_letter_frequencies = background_letter_frequencies
		self.strands = strands
	
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
		while not line.startswith("MOTIF"):
			self._proceed_next_line() 
			line = self._line
			if line is None:
				return None
		s1 = line.split()
		if len(s1) == 2:
			identifier = s1[1]
			alternate_name = None
		elif len(s1) == 3:
			identifier, alternate_name = s1[1:3]
		else:
			raise Exception("Invalid motif line - " + line)
		
		
		while not line.startswith("letter-probability matrix:"):
			self._proceed_next_line()
			line = self._line
		pattern = re.compile("\\s(\\S+)=\\s(\\S+)")
		matrix_dict = OrderedDict([i.groups() for i in pattern.finditer(line[26:])])
		funcs = {"alength": int, "w": int, "nsites": int, "E": float}
		matrix_dict = {k:(funcs[k](v) if k in funcs else v) for k, v in matrix_dict.items()}
		self._proceed_next_line()
		line = self._line
		
		pwm = []
		while line is not None and not line.startswith("URL") and not line.startswith("MOTIF"):
			pwm.append(list(map(float, line.split())))
			self._proceed_next_line()
			line = self._line
		if line.startswith("URL"):
			url = line.split()[1]
			self._proceed_next_line()
		else:
			url = None 
		return MEMEMotif(identifier, alternate_name, matrix_dict, pwm, url)
