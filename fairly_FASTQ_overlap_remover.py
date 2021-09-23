#!/usr/bin/env python
# coding: utf-8

# Copyright and Disclaimer

# Copyright ©2021, Cory David Dunn. All Rights Reserved. Permission to use, copy, and distribute this software for educational, \
# university-associated non-commercial research, and other not-for-profit purposes, without fee and without a signed licensing agreement, \
# is hereby granted, provided that the above copyright notice, this paragraph and the subsequent two paragraphs appear within all copies \
# and distributions. Contact the author at <cory.david.dunn@gmail.com> regarding the possibility of commercial licensing.

# IN NO EVENT SHALL CORY DAVID DUNN (THE AUTHOR) BE LIABLE TO ANY PARTY FOR ANY DAMAGES OR LOST PROFITS THAT MAY ARISE FROM THE USE \
# OF THIS SOFTWARE AND/OR ITS DOCUMENTATION, EVEN IF THE AUTHOR HAS BEEN NOTIFIED OF THE POTENTIAL OF SUCH DAMAGE.

# THE AUTHOR SPECIFICALLY REFUSES TO ACKNOWLEDGE, AND ABSOLUTELY REJECTS, ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, IMPLIED \
# WARRANTIES OF MERCHANTABILITY AND SUITABILITY FOR A PARTICULAR OBJECTIVE OR TASK. THE SOFTWARE AND ANY ASSOCIATED EXPLANATORY MATERIAL \
# IS STRICTLY PROVIDED "AS IS". THE AUTHOR HAS NO COMMITMENT TO, OR ANY RESPONSIBILITY FOR, SERVICE, IMPROVEMENTS, CHANGES, OR COMMUNICATIONS \
# RELATED TO THIS SOFTWARE OR ANY DERIVATIVE WORKS.

# Import dependencies

import pandas as pd
import logging
import time
import gzip
import argparse
from threading import Thread

version = '1.0'

def main():

    global start_time
    start_time = time.time() 

    ap = argparse.ArgumentParser()
    ap.add_argument('-1','--input_file_1',required=True,type=str,help="Input file 1 in FASTQ format.\n")
    ap.add_argument('-2','--input_file_2',required=False,type=str,default='NONE.FASTQ.GZ',help="Input file 2 in FASTQ format. Only required when mode is 'paired' \n")
    ap.add_argument('-s','--shift',required=False,type=int,default=0,help="Shift, or offset, from the 5' end of reads when choosing tags for sequence comparison.\n")
    ap.add_argument('-c','--chunk',required=False,type=int,default=1_000_000,help="The number of lines of the input file(s) read into memory when writing output files lacking duplicates. Very roughly optimized for 16GB RAM and a read length of 150.\n")
    ap.add_argument('-g','--gzip_compression_level',required=False,type=int,default=1,help="Gzip compression level (1-9). Default is '1' (most reasonable for speed and size). Choose '0' if working with uncompressed files.\n")

    args = vars(ap.parse_args())

    # Setting global variables

    global INPUT_FILENAME_1
    global INPUT_FILENAME_2
    global mode
    global offset
    global selected_chunk_size_in_sequences
    global compression_selected
    global mylogs

    INPUT_FILENAME_1 = args['input_file_1']
    INPUT_FILENAME_2 = args['input_file_2']
    
    # Set single or paired mode based upon input files

    if INPUT_FILENAME_2 == 'NONE.FASTQ.GZ':
        mode = 'single'
    elif INPUT_FILENAME_2 != 'NONE.FASTQ.GZ':
        mode = 'paired'
    
    offset = int(args['shift'])
    selected_chunk_size_in_sequences = args['chunk']
    compression_selected = args['gzip_compression_level']

    basename1 = INPUT_FILENAME_1.rsplit( ".")[0]
    basename2 = INPUT_FILENAME_2.rsplit( ".")[0]

    OUTPUT_FILENAME_1 = basename1 + '_fairly_output' + '.FASTQ'
    OUTPUT_FILENAME_2 = basename2+ '_fairly_output' + '.FASTQ'

    # Logging and streaming to console

    mylogs = logging.getLogger(__name__)
    mylogs.setLevel(logging.INFO)
    stream = logging.StreamHandler()
    stream.setLevel(logging.INFO)
    streamformat = logging.Formatter("%(message)s")
    stream.setFormatter(streamformat)
    mylogs.addHandler(stream)

    file = logging.FileHandler(basename1 + '_output.log')
    mylogs.addHandler(file)

    mylogs.info('___\n')
    mylogs.info('\nfairly_FASTQ_overlap_remover: Rapid removal of FASTQ entries apparently targeting the same coordinates on a template sequence\n')
    mylogs.info('Copyright ©2021, Cory David Dunn')
    mylogs.info('Contact: cory.david.dunn@gmail.com')
    mylogs.info('Version: ' + version)
    mylogs.info('___\n')
    
    # Check to make sure second input filename is selected in 'paired' mode

    if compression_selected != 0:

        OUTPUT_FILENAME_1 = OUTPUT_FILENAME_1 + '.gz'
        OUTPUT_FILENAME_2 = OUTPUT_FILENAME_2 + '.gz'

    # Print settings to console and log

    mylogs.info('Input filename 1: ' + INPUT_FILENAME_1)
    if mode == 'paired': 
        mylogs.info('Input filename 2: ' + INPUT_FILENAME_2)
    if offset != 0:
        mylogs.info("Character offset from 5' end(s): " + str(offset))
    if selected_chunk_size_in_sequences != 1_000_000:
        mylogs.info("'Chunk' size when writing clean output files: " + str(selected_chunk_size_in_sequences))
    if compression_selected != 1:
        mylogs.info("Compression level selected for gzip file (0 is no compression): " + str(compression_selected))
    mylogs.info('Output filename 1: ' + OUTPUT_FILENAME_1)
    if mode == 'paired': 
        mylogs.info('Output filename 2: ' + OUTPUT_FILENAME_2)

    mylogs.info('___\n')
    mylogs.info('Loading input sequence tags and quality scores.')

    if mode == 'single':
        
        # Here, directly take the 40bp tag and associated quality characters
        
        merge_FASTA20_to_FASTA40, mergeQ20_to_Q40 = import_FASTQ_tags(INPUT_FILENAME_1,mode)

        mylogs.info('40-mers retrieved from single provided file.')
        mylogs.info("Elapsed time: ~ " + str(int(time.time() - start_time)) + " seconds.")

    if mode == 'paired':

        list_of_FASTA1_lines, list_of_FASTQ1_QUAL_lines = import_FASTQ_tags(INPUT_FILENAME_1,mode)
        list_of_FASTA2_lines, list_of_FASTQ2_QUAL_lines = import_FASTQ_tags(INPUT_FILENAME_2,mode)
        
        # For paired set, merge the 20-bp tags into one string (FASTA and quality characters)
        
        merge_FASTA20_to_FASTA40  = [i + j for i, j in zip(list_of_FASTA1_lines, list_of_FASTA2_lines)]
        
        del list_of_FASTA1_lines
        del list_of_FASTA2_lines

        mylogs.info('Merged 20-mers from each of the input files.')
        mylogs.info("Elapsed time: ~ " + str(int(time.time() - start_time)) + " seconds.")

        mergeQ20_to_Q40 = [i + j for i, j in zip(list_of_FASTQ1_QUAL_lines, list_of_FASTQ2_QUAL_lines)]
        
        del list_of_FASTQ1_QUAL_lines
        del list_of_FASTQ2_QUAL_lines

    mylogs.info('Computing...')

    # Convert quality characters to numerical scores (sums across ASCII values) for all 40bp tags

    output_numer_Q_scores = FASTQ_quality_to_number(mergeQ20_to_Q40)

    # Merge the FASTA tags and quality scores into the same dataframe

    FASTA_DF = pd.DataFrame({ 'MERGE_FASTA': merge_FASTA20_to_FASTA40, 'MERGE_Q_number': output_numer_Q_scores })

    global length_of_FASTQ_in_sequences
    length_of_FASTQ_in_sequences = len(FASTA_DF)
    
    FASTA_DF.to_csv('ONE_OFFSET.csv')

    del merge_FASTA20_to_FASTA40
    del mergeQ20_to_Q40

    # Ensure that all tags are of length 40 and ensure that the tag harbors no 'N'(potentially indicating read quality and length)

    FASTA_DF['LEN'] = FASTA_DF['MERGE_FASTA'].apply(len).astype('int8')
    FASTA_DF['N_PRESENT'] = FASTA_DF['MERGE_FASTA'].str.contains('N')
    FASTA_DF = FASTA_DF[FASTA_DF['LEN'] == 40]
    FASTA_DF = FASTA_DF[FASTA_DF['N_PRESENT'] == False]

    # Sort the resulting dataframe 

    FASTA_DF.sort_values(by = 'MERGE_Q_number', ascending = False, inplace=True)
    FASTA_DF.drop_duplicates(subset='MERGE_FASTA', keep='first', inplace=True, ignore_index=False)

    global list_of_entries_to_save
    list_of_entries_to_save = FASTA_DF.index.to_list()
    total_left = len(list_of_entries_to_save)

    mylogs.info('Total FASTQ entries analyzed: ' + str(length_of_FASTQ_in_sequences))
    mylogs.info('Total FASTQ entries removed: ' + str(length_of_FASTQ_in_sequences - total_left))
    mylogs.info('Total FASTQ entries remaining: ' + str(total_left) + '.... writing to disk')

    del FASTA_DF

    mylogs.info('Writing output file(s).')

    # Write at least one file and...

    file_1_write = Thread(target=write_clean_FASTQ_to_GZIP, args=(INPUT_FILENAME_1,OUTPUT_FILENAME_1))
    file_1_write.start()

    # ...if the FASTQ run was in paired more, write the second file. Save simultaneously using threading.

    if mode == 'paired':
        
        file_2_write = Thread(target=write_clean_FASTQ_to_GZIP, args=(INPUT_FILENAME_2,OUTPUT_FILENAME_2))
        file_2_write.start()

    return()

def import_FASTQ_tags(inputfilename, mode):
    
    if mode == 'paired':
        chars_to_take = 20
        
    elif mode == 'single':
        chars_to_take = 40
    
    inputfile_F = gzip.open(inputfilename,'rt')

    list_of_FASTA_lines = []
    list_of_FASTQ_QUAL_lines = []

    # Move through FASTQ entries, save a 20-bp tag from each of pair or 40-bp tag from single read set
    
    count = 1
    
    for line in inputfile_F:
        
        # Take from the FASTA line
        
        if count % 4 == 2:
            twentychar5 = line[offset : (chars_to_take + offset)].upper()
            list_of_FASTA_lines.append(twentychar5)

        # Take from the quality character line

        elif count % 4 == 0:
            Qchar5 = line[offset : (chars_to_take + offset)]
            list_of_FASTQ_QUAL_lines.append(Qchar5)

        count += 1
    
    inputfile_F.close()
    
    return(list_of_FASTA_lines, list_of_FASTQ_QUAL_lines)

def FASTQ_quality_to_number(input_list_of_coded_Q_scores):

    list_of_Q_sums = []

    # Moving through each quality character string
    
    for stringQ in input_list_of_coded_Q_scores:
        
        string_to_sum = 0
        
        # For each character, convert to ASCII and sum across the 40bp tag
        
        for x in stringQ:
            asciicode = int(ord(x))
            string_to_sum += asciicode
        
        list_of_Q_sums.append(string_to_sum)
        
    return(list_of_Q_sums)

def write_clean_FASTQ_to_GZIP(inputfilename,outputfilename):
    
    # Open the input and output files. Gzip compression level is set at 1 [high compression levels slows down save considerably at little size benefit]
    
    inputfile_F = gzip.open(inputfilename,'rt')
    outputfile_F = gzip.open(outputfilename,mode='wt',compresslevel=compression_selected)

    # Loading the entire FASTQ file(s) into memory is implausible, so the input file must be handled in chunks.
    
    number_of_globals = ((length_of_FASTQ_in_sequences) // selected_chunk_size_in_sequences) + 1

    for i in range(number_of_globals):
        
        # Delete previous chunk dataframe, if it exists.
        
        try: 
            del FASTQ_DF
        except:
            pass
        
        # Build the next chunk dataframe, while also keeping track of the indexes compared to the initial tag (40bp) dataframe.
        
        section_number = 'section_' + str(i)
        globals()[section_number] = []
        new_border_left = i * selected_chunk_size_in_sequences
        new_border_right = ((i+1) * selected_chunk_size_in_sequences) - 1
        globals()[section_number][:] = [x for x in list_of_entries_to_save if (x >= new_border_left) and (x <= new_border_right)]

        
        # Given the window size, load all four lines for each FASTQ entry into a dataframe.
        
        LINEZ1 = []
        LINEZ2 = []
        LINEZ3 = []
        LINEZ4 = []

        count = 1

        for line in inputfile_F:
            
            # If the source file has ended, break the loop.

            if not line:
                break



            if count % 4 == 1:
                LINEZ1.append(line)
            elif count % 4 == 2:
                LINEZ2.append(line)
            elif count % 4 == 3:
                LINEZ3.append(line)
            elif count % 4 == 0:
                LINEZ4.append(line)

            count += 1

            if count == (selected_chunk_size_in_sequences * 4) + 1:
                break

        dict_to_add = {'LINEZ1': LINEZ1, 'LINEZ2': LINEZ2, 'LINEZ3': LINEZ3, 'LINEZ4': LINEZ4} 

        FASTQ_DF = pd.DataFrame(dict_to_add)

        del LINEZ1
        del LINEZ2
        del LINEZ3
        del LINEZ4

        # Given the index of non-duplicates originally identified (length of tag is 40bp and no 'N' characters in tag),
        # and the current window, generate a clean dataframe of FASTQ lines.
        
        current_index_list = FASTQ_DF.index.to_list()
        revised_index_list = [z+(i*selected_chunk_size_in_sequences) for z in current_index_list]
        maximum_sequence_index = length_of_FASTQ_in_sequences - 1
        revised_index_list_corrected = [x for x in revised_index_list if (x <= maximum_sequence_index)]
        FASTQ_DF = FASTQ_DF.rename(index=dict(zip(current_index_list,revised_index_list_corrected)))
        to_keep = list(set(globals()[section_number]))
        
        # Stack the clean FASTQ dataframe and write as a string to the output file
        
        FASTQ_DF_clean = FASTQ_DF[FASTQ_DF.index.isin(to_keep)]
        FASTQ_DF_clean_stacked = FASTQ_DF_clean.stack()
        FASTQ_DF_clean_stacked_list = FASTQ_DF_clean_stacked.to_list()
        FASTQ_DF_clean_stacked_string = ''.join(FASTQ_DF_clean_stacked_list)
        outputfile_F.write(FASTQ_DF_clean_stacked_string)   
    
    # Close input and output files, and report on elapsed time at end of each file save.
    
    inputfile_F.close()
    outputfile_F.close()
    
    mylogs.info("Output file: " + outputfilename + " closed. Elapsed time: ~ " + str(int(time.time() - start_time)) + " seconds.")
    
    return

# Run the main function

if __name__ == '__main__':
    main()

