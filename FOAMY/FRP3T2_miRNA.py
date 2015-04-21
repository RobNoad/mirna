# Author: Rodney Kincaid
# Date: Nov 13, 2013

from Bio.Seq import Seq, reverse_complement, transcribe, back_transcribe
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from operator import add, mul
import numpy as np
import os, sys, re, argparse, subprocess

# Regexes
# Box B sequences: GTTCNANNC, GGTTSGNG, GKWCAAGTC
# Added GSTYS
# Box B sequences: GTTCNANNC, GGTTSGNG, GKWCAAGTC # Del GCTG and CCCT
box_b_re = re.compile('G[TG][TA]CAAGTC[ATCG]|[ATCG]GTT[CG][ATCG]{5}|[GA][GA]TCGA[GA][ATCG]{3}|GGAA[TC]GA[ATCG]{3}|GGGAAGG[ATCG]{3}|GGCTCTGG[ATCG]{2}')
# Box A sequences: TRGNNNNNNGR, TRGNNNNNGR, TRGCTC
box_a_re = re.compile('T[AG]G[ATCG]{14}')
# Terminator sequences: TTTT
box_t_re = re.compile('T[ATCG]TT|TT[ATCG]T')
box_u_re = re.compile('U[AUCG]UU|UU[AUCG]U')
# Stem
stem_re = re.compile('([(][.(]{19,26})([(][.]{4,18}[)])([.)]{19,26}[)])')
#Internal U's
internal_t_re = re.compile('TTTT|TT[ATCG]TTT|TTT[ATCG]TT')
#Internal bulge
internal_bulge_re = re.compile('[)(][.]{5,15}[)(]')

# Maximum stem loop fold energy for reporting
maximum_fold_energy = -15.0

# Maximum number of nucleotides overlapping a CDS
max_CDS_overlap = 50

# Maximum window lengths
window_a_to_b = 70
window_a_to_t = 150
window_b_to_t = 70

# Process command line arguements
cmd_parser = argparse.ArgumentParser(description='Search for ncRNAs with Pol III, type 2 promoters in a file of sequences.')
cmd_parser.add_argument('file_name', metavar='file_name',help='Name of input sequence file')
#cmd_parser.add_argument('-genbank', dest='genbank_format', action='store_const', const='genbank', default='fasta', help='Input file in GenBank format (default: fasta)')
#cmd_parser.add_argument('-CDS', dest='CDS_filter', action='store_true', default=False, help='Filter out results that overlap with CDS (*GenBank Format must also be selected*) (default: off)')
cmd_parser.add_argument('-wab', dest='wab', action='store', type=int, default=window_a_to_b, help='Max distance Box B to look for Box A (default: '+str(window_a_to_b)+')')
cmd_parser.add_argument('-wbt', dest='wbt', action='store', type=int, default=window_b_to_t, help='Max distance Box B to look for Terminator (default: '+str(window_b_to_t)+')')
cmd_parser.add_argument('-CDS_overlap', dest='CDS_overlap', action='store', type=int, default=max_CDS_overlap, help='Max overlap with a CDS (default: '+str(window_b_to_t)+')')
cmd_parser.add_argument('-stem_dG', dest='stem_dG', action='store', type=float, default=maximum_fold_energy, help='Max energy cutoff for a miRNA stem (default: '+str(maximum_fold_energy)+')')
#cmd_parser.add_argument('-stem_filter', dest='stem_filter', action='store', type=str, choices=['none','mirna','varna'], default='none', help='Apply a stem filter (default: none)')
args = cmd_parser.parse_args()

seqs_in = args.file_name
window_a_to_b = args.wab
window_b_to_t = args.wbt
max_CDS_overlap = args.CDS_overlap
maximum_fold_energy = args.stem_dG

# Spawn RNAfold
rnafold = subprocess.Popen(["RNAfold", "--noPS"], bufsize=10000, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

# Read in and process seqs
seqs_input = open(seqs_in,"r")
for seq in SeqIO.parse(seqs_input, "genbank"):
	
	# Create a dictionary of the sequence and the reverse complement
	ss = {}
	ss['+'] = str(back_transcribe(str(seq.seq)))
	ss['-'] = str(reverse_complement(ss['+']))

	# Extract CDS regions from annotations
	ss_cds = [0 for x in range(len(ss['+']))]
	for f in seq.features :
		if f.type == 'CDS':
			if f.sub_features != []:
				# Subfeatures of the CDS -> Splicing!
				for g in f.sub_features:
					for x in range(g.location.start.position,g.location.end.position):
						ss_cds[x] = 1
			else:
				# No subfeatures of the CDS -> Unspliced in the coding region
				for x in range(f.location.start.position,f.location.end.position):
					ss_cds[x] = 1
	ss_cds_rev = ss_cds[::-1]

	# Process the sequence and the reverse complement
	for sn,s in ss.iteritems():
		length_of_s = len(s)
		# Find and process each Box B like sequence
		for mb in box_b_re.finditer(s):
			mb_pos = mb.start()

			if mb_pos > window_a_to_b:
				window_start = mb_pos - window_a_to_b
			else:
				window_start = 0
	
			window_end = mb_pos+8

			# Find the Box A like sequence within the upstream window of this Box B like element
			for ma in box_a_re.finditer(s[window_start:mb_pos]):
				ma_pos = window_start+ma.start()
				
				# Find the terminator sequence \
				#mt = box_t_re.search(s[mb_pos-20:mb_pos+window_b_to_t]) 
				#if mt and not internal_t_re.search(s[ma_pos:ma_pos+40]) :
				for mt in box_t_re.finditer(s[mb_pos-20:mb_pos+window_b_to_t]):
					mt_pos =  mt.start()+mb_pos-20

					if not internal_t_re.search(s[ma_pos:mt_pos]):

						# Found a BLV pol III like sequence, print it and some relevant information
						#mt_pos =  mt.start()+mb_pos-20
						outstring = []
						outstring.append(seq.name)
						outstring.append(seq.description)
						outstring.append(sn)
						if sn == '+':
							loc = ":".join([str(ma_pos),str(mb_pos),str(mt_pos)])
							cds = ss_cds
						else:
							loc1 = ":".join([str(ma_pos),str(mb_pos),str(mt_pos)])
							loc2 = ":".join([str(length_of_s-mt_pos),str(length_of_s-mb_pos),str(length_of_s- ma_pos)])
							loc = loc1 + ' (' + loc2 + ')'
							cds = ss_cds_rev

						outstring.append(loc)

						if ma_pos-14 >= 0:
							the_seq = s[ma_pos-14:mt_pos+4]
							start_num = ma_pos-14
							end_num = window_end
						else:
							the_seq =  s[0:window_end]
							start_num = 0
							end_num = window_end
					
						outstring.append(str(len(the_seq)))
						outstring.append(ma.group(0))
						outstring.append(mb.group(0))
						outstring.append(mt.group(0))
					
						if sum(cds[start_num:end_num]) < max_CDS_overlap:
							# Send the sequence to the RNAfold process
							rnafold.stdin.write(the_seq+"\n")
							rnafold.stdin.flush()
					
							# Read first line of RNAfold output (sequence)
							rna_seq = rnafold.stdout.readline().strip()
							outstring.append(rna_seq)

							# Read second line of RNAfold output (fold and energy)
							rnaout = rnafold.stdout.readline()

							f = rnaout.strip().split(' (')
							fold  = f[0].strip()
							outstring.append(fold)
							energy = f[1].strip(')')
							outstring.append(energy)

							# Find miRNA like stemloops within the larger fold
							found_miRNA= False
							hp = []
							for stem in stem_re.finditer(fold):

								if not internal_bulge_re.search(stem.group(1)) and not internal_bulge_re.search(stem.group(3)):
									# Send the sequence to the RNAfold process
									rnafold.stdin.write(rna_seq[stem.start(0):stem.end(0)]+"\n")
									rnafold.stdin.flush()
						
									# Read first line of RNAfold output (sequence)
									srna_seq = rnafold.stdout.readline().strip()

									# Read second line of RNAfold output (fold and energy)
									srnaout = rnafold.stdout.readline()				

									sf = srnaout.strip().split(' (')								
									#print maximum_fold_energy, ":", float(sf[1].strip(')').strip())
									if float(sf[1].strip(')').strip()) <= maximum_fold_energy and box_u_re.search(rna_seq[+stem.end(0)-4:stem.end(0)+6]):
										hp.append(srna_seq+':'+sf[0].strip()+':'+sf[1].strip(')'))
										found_miRNA= True

							if found_miRNA:
								outstring.append(','.join(hp))
								print "\t".join(outstring)

							



