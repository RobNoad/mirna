# Author: Rodney Kincaid
# Date: May 3, 2011
# Last Modified: Octover 16, 2013
# Requirements
# -- Biopython : http://http://biopython.org
# -- RNAfold : http://http://www.tbi.univie.ac.at/RNA/ (RNAfold executable must be in Path)
# Usage Example
# -- python find_miRNA_stemloop.py genbank_formatted_sequences.gb > output_file.tab
# Output
# -- Prints to stdout tab separated Sequene Name, Sequence Description, Start Position, End Position, Sequence, Secondary Structure, Energy
#

from Bio.Seq import Seq, reverse_complement, transcribe, back_transcribe
from Bio import SeqIO
import re, argparse, subprocess

# Regular Expressions #
# Stem
stem_re = re.compile('([(][.(]{24,32})([(][.]{4,18}[)])([.)]{24,32})[)]')

# Maximum stem loop fold energy for reporting
maximum_fold_energy = -22.0

# Maximum number of nucleotides overlapping a CDS
max_CDS_overlap = 1

# Window size to fold
window = 200

# Window increment step
step = 50

# Process command line arguements
cmd_parser = argparse.ArgumentParser(description='Search for stemloop structures in non-CDS regions.')
cmd_parser.add_argument('file_name', metavar='file_name',help='Name of input Genbank formatted sequence file')
cmd_parser.add_argument('-CDS_overlap', dest='CDS_overlap', action='store', type=int, default=max_CDS_overlap, help='Max overlap with a CDS (default: 1)')
cmd_parser.add_argument('-stem_dG', dest='stem_dG', action='store', type=float, default=maximum_fold_energy, help='Max energy cutoff for a miRNA stem (default: -22)')
args = cmd_parser.parse_args()

# Set any relevant command line arguments
seqs_in = args.file_name
max_CDS_overlap = args.CDS_overlap
maximum_fold_energy = args.stem_dG

# Spawn RNAfold
rnafold = subprocess.Popen(["RNAfold", "-noPS"], bufsize=10000, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

# Read in and process seqs
seqs_input = open(seqs_in,"r")
for seq in SeqIO.parse(seqs_input, "genbank"):
	
	# Create a dictionary of the sequence and the reverse complement
	ss = {}
	ss['+'] = str(back_transcribe(str(seq.seq)))

	# Uncomment to search both strands
	#ss['-'] = str(reverse_complement(ss['+']))

	# Extract CDS regions from annotations
	ss_cds = [0 for x in range(len(ss['+']))]
	for f in seq.features :
		if f.type == 'CDS':
			if f.sub_features != []:
				# Subfeatures of the CDS -> Splicing
				for g in f.sub_features:
					for x in range(g.location.start.position,g.location.end.position):
						ss_cds[x] = 1
			else:
				# No subfeatures of the CDS -> Unspliced in the coding region
				for x in range(f.location.start.position,f.location.end.position):
					ss_cds[x] = 1
	ss_cds_rev = ss_cds[::-1]

	# Process the sequences
	for sn,s in ss.iteritems():
		length_of_s = len(s)
		this_start = 1
		this_end = this_start + window
		for x in range(this_start, length_of_s-window, step):
			outstring = []
			outstring.append(seq.name)
			outstring.append(seq.description)
			if sn == '+':
				cds = ss_cds
			else:
				cds = ss_cds_rev
			
			start_num = x
			end_num = x + window
			the_seq = s[start_num:end_num]

			
			if sum(cds[start_num:end_num]) < max_CDS_overlap:
				# Send the sequence to the RNAfold process
				rnafold.stdin.write(the_seq+"\n")
				rnafold.stdin.flush()
			
				# Read first line of RNAfold output (sequence)
				rna_seq = rnafold.stdout.readline().strip()

				# Read second line of RNAfold output (fold and energy)
				rnaout = rnafold.stdout.readline()

				f = rnaout.strip().split(' (')
				fold  = f[0].strip()
				energy = f[1].strip(')')
			
				# Find miRNA like stemloops within the larger fold
				found_miRNA= False
				hp = []
				for stem in stem_re.finditer(fold):
					# Send the sequence to the RNAfold process
					rnafold.stdin.write(rna_seq[stem.start(0):stem.end(0)]+"\n")
					rnafold.stdin.flush()
				
					# Read first line of RNAfold output (sequence)
					srna_seq = rnafold.stdout.readline().strip()

					# Read second line of RNAfold output (fold and energy)
					srnaout = rnafold.stdout.readline()

					sf = srnaout.strip().split(' (')					
					
					# Check for energy requirement and large internal bulges
					if float(sf[1].strip(')').strip()) <= maximum_fold_energy and  stem.group(1).find('....') == -1 and  stem.group(3).find('....') == -1:						
						hp.append(str(start_num+stem.start(0))+"\t"+str(start_num+stem.end(0))+"\t"+srna_seq+'\t'+sf[0].strip()+'\t'+sf[1].strip(')'))
						found_miRNA= True
			

				if found_miRNA:
					for mirnastem in hp:
						print "\t".join(outstring)+"\t"+mirnastem

					



