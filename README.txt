# fairly_FASTQ_overlap_remover

A Python script for removal of FASTQ reads (single or paired end) that appear to target the same template

python fairly_FASTQ_overlap_remover.py -h
usage: fairly_FASTQ_overlap_remover.py [-h] -1 INPUT_FILE_1 [-2 INPUT_FILE_2] [-s SHIFT] [-c CHUNK] [-g GZIP_COMPRESSION_LEVEL]

optional arguments:
  -h, --help            show this help message and exit
  -1  --input_file_1
                        Input file 1 in FASTQ format.
  -2 --input_file_2
                        Input file 2 in FASTQ format. Only required when mode is 'paired'
  -s --shift
                        Shift, or offset, from the 5' end of reads when choosing tags for sequence comparison.
  -c --chunk
                        The number of lines of the input file(s) read into memory when writing output files lacking duplicates. Very
                        roughly optimized for 16GB RAM and a read length of 150.
  -g --gzip_compression_level
                        Gzip compression level (1-9). Default is '1' (most reasonable for speed and size). Choose '0' if working with
                        uncompressed files.
                        
The purpose of this script is to remove likely duplicates before mapping FASTQ to a reference genome. 

It is implemented in Python 3 (version 3.9.4) with a Pandas dependency (Pandas 1.2.4 tested here).

For cleaning one input file (here representing a single read), only enter one input file ('-1') . For cleaning paired reads, please enter two input files ('-1' and '-2' arguments.
                       
For single-end read sets (single input file), the script uses 40bp tags obtained from the 5' end of the read to determine whether the entry may be a duplicate. For paired-end read sets (two input files), the script takes 20bp tags obtained from the 5' end of each read set, then concatenates them to determine whether the entry may be a duplicate. 

When likely duplicates are identified, that with the highest quality score is kept in the output file(s).

Please note that if 'N' is found in any analyzed tag, that sequence is never saved in the output file.

The 'shift parameter is the number of characters from the 5' end of input sequences from which tags should be extracted for duplicate determinations.

The 'chunk' parameter is the number of FASTQ text lines loaded into memory from input files while writing to output files. The current setting (1,000,000) is roughly optimized for 16 GB of RAM, but can be altered, if needed, by the user.

The 'gzip_compression_level' has large effects on the rate at which an output file is saved. A setting of '1' allows rapid saving, and higher settings offer little gain in file compression.

__

Copyright and Disclaimer

Copyright ©2021, Cory David Dunn. All Rights Reserved. Permission to use, copy, and distribute this software for educational, university-associated non-commercial research, and other not-for-profit purposes, without fee and without a signed licensing agreement, is hereby granted, provided that the above copyright notice, this paragraph and the subsequent two paragraphs appear within all copies and distributions. Contact the author at <cory.david.dunn@gmail.com> regarding the possibility of commercial licensing. 

IN NO EVENT SHALL CORY DAVID DUNN (THE AUTHOR) BE LIABLE TO ANY PARTY FOR ANY DAMAGES OR LOST PROFITS THAT MAY ARISE FROM THE USE OF THIS SOFTWARE AND/OR ITS DOCUMENTATION, EVEN IF THE AUTHOR HAS BEEN NOTIFIED OF THE POTENTIAL OF SUCH DAMAGE. 

THE AUTHOR SPECIFICALLY REFUSES TO ACKNOWLEDGE, AND ABSOLUTELY REJECTS, ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, IMPLIED WARRANTIES OF MERCHANTABILITY AND SUITABILITY FOR A PARTICULAR OBJECTIVE OR TASK. THE SOFTWARE AND ANY ASSOCIATED EXPLANATORY MATERIAL IS STRICTLY PROVIDED "AS IS". THE AUTHOR HAS NO COMMITMENT TO, OR ANY RESPONSIBILITY FOR, SERVICE, IMPROVEMENTS, CHANGES, OR COMMUNICATIONS RELATED TO THIS SOFTWARE OR ANY DERIVATIVE WORKS.
