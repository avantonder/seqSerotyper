# seqSerotyper.R
seqSerotyper.R is a tool for assigning serotypes to pneumococcal whole genome sequences (WGS).

Author: Andries J van Tonder

e-mail: andries.vantonder@stcatz.ox.ac.uk

Copyright (C) 2016 University of Oxford

#General remarks
NOTE: This script has only been tested on Linux so probably wont work on Windows or OSX.

#Disclaimers

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Download
git clone --recursive https://github.com/avantonder/seqSerotyper/

#Dependencies
R

seqinr (R library)

transeq (BioLinux)

mview (BioLinux)

blat (Biolinux)

blastall (Biolinux)

#Running seqSerotyper.R
Usage: Rscript seqSerotyper.R -dataFile <input file> -srcDir <path>

Options:
  
-dataFile:	A file containing sample ids and file paths to fasta files

-srcDir:	The path of the directory that contains the reference files and related scripts


e.g. Rscript seqSerotyper.R -dataFile pneumo_fasta_ids.txt -srcDir /home/user/Software/seqSerotyper/bin

#Input file

The file containing the sample ids and file paths to the pneumococcal fasta files to be serotypes should 
be in the following tab-separated format:

id	filePath

ERR065287_A46670	/media/data/WGS_sequences/ERR065287_A46670.fa

#Output file

The script will create a file called Predicted_serotypes.txt in the seqSerotyper directory.  
The file lists the strain ids and associated serotypes.
