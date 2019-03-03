"""
Sequencer extracts sequencing reads from result files. Applies a 
user-specified cutoff to the sequencing files for a signal threshold and 
discards everything below that by changing those basecalls to a gap ('-').

For help, run it with -h option. 
"""

import argparse
import sys
import logging

#Set up logging configuration
logging.basicConfig(level=logging.INFO,
                    filename = 'sequencer.log',
                    filemode = 'w',
                    format='%(name)s:%(levelname)s:%(message)s')

#Set up argument parser and extract the command-line arguments
parser = argparse.ArgumentParser(description='Process sequencing files')

parser.add_argument('-n',
                    metavar='nucleotides', 
                    help='Input file in FASTA format containing the nucleotide sequence')
parser.add_argument('-s',
                    metavar='scores',
                    help='Input file in FASTA format containing the nucleotide scores')
parser.add_argument('-f',
                    metavar='signal_cutoff',
                    type=int,
                    default=20,
                    help='Signal cutoff value under which a base call is rejected. Default value is 20')
parser.add_argument('-o',
                    metavar='output',
                    help='Output file for the results')
parser.add_argument('-r', 
                    help='Specifies that the sequencing file needs to be reverse complemented',
                    action='store_true')

arguments = parser.parse_args()

"""Fetch the sequence from a FASTA file containing a single entry.

Read a FASTA file containing a single sequence and return the sequence information
If the file contains more than one sequence, the extra sequences are not processed

Args:
  filename: Path to a file to be processed.

Returns:
  String containing the FASTA content with newlines stripped.
"""
def ReadFileContent(filename):
  with open(filename, 'r') as file:
    header_found = False;
    for line in file:
      #Desired content is always on the next line after the FASTA header, which 
      #itself always starts with a '>' character.
      if line.startswith('>'):
        header_found = True
        continue

      if header_found:
        return line.strip('\n')


"""Filter nucleotides based on the signal cutoff.

Filter two lists containing nucleotide sequence and corresponding scores based
on the specified quality cutoff. Any nucleotides with a score lower than the
cutoff is replalced by a gap ('-')

Args:
   nucleotide_list: A list of nucleotides
   scores_list: A list of nucleotide scores that correspond to the nucleotides in 
                 nucleotide_list
   quality_cutoff: Equal to -10*log10(P) value where P is the probability of a miscall.
                    For example, a signal cutoff value of 30 means there is a 
                    1/1000 chance of miscalled base. (log10(1/1000) = -3 * -10 = 30)

Returns:
   A list of filtered nucleotides with a score equal to or above the cutoff value.
   All base calls below the threshold are replaced by gaps ('-').
"""
def FilterSequence(nucleotide_list, quality_list, quality_cutoff):
  logging.debug(f'Input nucleotides list: {nucleotide_list}')
  logging.debug(f'Input quality scores list: {quality_list}')

  #Check whether either of the lists are empty (could indicate bad sequencing file)
  if not nucleotide_list:
    logging.error('Nucletide sequence list is empty, aborting')
    return

  if not quality_list:
    logging.error('Quality scores list is empty, aborting')
    return

  logging.debug(f'Nucleotide list length: {len(nucleotide_list)}, quality list length: {len(quality_list)}')
  if len(nucleotide_list) != len(quality_list):
    raise ValueError('Nucleotide and scores lists are not equal in length, indicating malformed data')


  logging.info("Signal threshold has been set to " + str(quality_cutoff))
  position = 0
  filtered_sequence = []

  for quality_score in quality_list:
    if int(quality_score) >= quality_cutoff:
      filtered_sequence.append(nucleotide_list[position])

    if int(quality_score) < quality_cutoff:
      filtered_sequence.append('-')

    position += 1

  return filtered_sequence

"""Write the results into a text file

Args:
  output_file: Name of the file to write results into
  data: String of data to be written
"""
def WriteOutput(output_file, data):
  with open(output_file, 'w') as output:
    output.write(data)
  print("Writing results into " + output_file)

"""Reverse complement a DNA sequence

Reverse complement a sequence. Non-DNA characters are unmodified and retained.

Args:
  sequence: DNA sequence

Returns:
  Reverse complemented DNA string
"""
def ReverseComplement(sequence):
  reverse_sequence = sequence[::-1]
  reverse_complement = ""
  for c in reverse_sequence:
    if c == 'A': reverse_complement += 'T'
    elif c == 'a': reverse_complement += 't'
    elif c == 'G': reverse_complement += 'C'
    elif c == 'g': reverse_complement += 'c'
    elif c == 'C': reverse_complement += 'G'
    elif c == 'c': reverse_complement += 'g'
    elif c == 'T': reverse_complement += 'A'
    elif c == 't': reverse_complement += 'a'
    else: reverse_complement += c

  return reverse_complement


def main():
  logging.info(f'Nucleotide sequence file: {arguments.n}')
  logging.info(f'Nucleotide scores file: {arguments.s}')

  #Nucleotide quality scores comes as a string of space-separated numbers so we need 
  #to split them by space to get a list of individual scores.
  nucleotide_scores = ReadFileContent(arguments.s).split(' ')
  nucleotide_sequence = ReadFileContent(arguments.n)

  if (arguments.r):
    logging.info("Reverse complement specified")

    nucleotide_scores = nucleotide_scores[::-1]

    nucleotide_sequence = ReverseComplement(nucleotide_sequence)

  nucleotide_sequence = list(nucleotide_sequence)

  signal_cutoff = arguments.f
  filtered_sequence = FilterSequence(nucleotide_sequence, nucleotide_scores, signal_cutoff)

  output_file = arguments.o
  nucleotide_data = ''.join(filtered_sequence)
  WriteOutput(output_file, nucleotide_data)

if __name__ == '__main__':
  main()