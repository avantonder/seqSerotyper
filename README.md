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
If you have git installed:

  git clone --recursive https://github.com/avantonder/seqSerotyper/

If not, click on the Clone or Download button above and select Download zip.  Move the file to your desired location and unzip the file and then cd into the seqSerotyper directory:

    unzip seqSerotyper-master.zip
  
    cd seqSerotyper-master

#Dependencies
##R and Perl

R can be downloaded from:

  https://www.r-project.org/

##seqinr (R library)

To install the seqinr library:

  From within R:

  If you have root access type sudo R in a terminal.  Once R loads type:

    install.packages("seqinr")

  and select a mirror from the dialogue window that pops up.  The library should then install and you can close R by typing quit()

  If install.packages doesn't work due to an issue with the version of R you have, you can install the library from source:

  Download the library source file from:

    https://cran.r-project.org/web/packages/seqinr/index.html

  In a terminal cd into the directory you downloaded the source file to and type the following:

    sudo R CMD INSTALL seqinr_3.3-1.tar.gz

seqSerotyper.R comes bundled with pre-compiled 64 bit executables for blat, blastall, mview and transeq.

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
The file lists the strain ids and associated serotypes:

  id  Serotype
  
  ERR065287_A46670  19A
