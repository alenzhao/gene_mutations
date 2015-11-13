
## codon-aminoacid table 
gencode = {
		'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
		'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
		'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
		'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
		'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
		'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
		'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
		'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
		'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
		'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

AAnames = {'A': 'Ala','I': 'Ile', 'L': 'Leu', 'V': 'Val', 'F': 'Phe', 'W': 'Trp', 'Y': 'Tyr', 'N': 'Asn', 'C': 'Cys', 'Q': 'Gln', 'M': 'Met', 'S': 'Ser', 'T': 'Thr', 'D': 'Asp', 'E': 'Glu', 'R': 'Arg' , 'H': 'His', 'K': 'Lys', 'G': 'Gly', 'P': 'Pro', '*': 'Ter' }; 

revcomp = {'A':'T', 'T':'A', 'C':'G', 'G':'C','a':'T', 't':'A', 'c':'G', 'g':'C', 'N':'N','n':'N' } 


def rev_seq(seq):
	revcomp = {'A':'T', 'T':'A', 'C':'G', 'G':'C','a':'T', 't':'A', 'c':'G', 'g':'C', 'N':'N','n':'N' }
	return "".join([revcomp[base] for base in reversed(seq)])