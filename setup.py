from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
	readme = readme_file.read()

requirements = ["genomictools>=0.0.8"]

setup(
	name="biodata",
	version="0.1.5",
	author="Alden Leung",
	author_email="alden.leung@gmail.com",
	description="A python package for common biological data I/O",
	long_description=readme,
	long_description_content_type="text/markdown",
	url="https://github.com/aldenleung/biodata/",
	packages=find_packages(),
	install_requires=requirements,
	classifiers=[
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"Programming Language :: Python :: 3.9",
		"Programming Language :: Python :: 3.10",
		"Programming Language :: Python :: 3.11",
		"License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
	],
)
