import numpy as np
import pandas

#This code will find partially and perefectly matched variants in two vcf files based on reference alleles and alternate alleles.

ref_file="/Users/rohamrazaghi/Documents/vcf/data.vcf"       # path to the individual VCF file
mut_file="/Users/rohamrazaghi/Documents/vcf/mutations.vcf"  # path to the list of mutation diseases/ second individual


# Open and read both files and split each line to strings
ref = open(ref_file)
mut = open(mut_file)
s_ref = ref.read()
s_mut = mut.read()
ref.close()
mut.close()

s_ref_split = s_ref.split('\n')
s_mut_split = s_mut.split('\n')
# define parameters
flag = 0
i_f = 0
i_g = 0
partially_matched = []
perfectly_matched = []
c = 0
# make a list of chromosomes for the ease of comparison later in the while loop
ch_set_g = set()
ch_name = []
for i in range(22):
    ch_name = np.append(ch_name,'%g'%(i+1))
ch_name = np.append(ch_name,'X')
ch_name = list(ch_name)


while flag == 0:                                # while loop
    
    if s_mut_split[i_g][0] == "#":                      # Neglecting the comment lines in both files
        i_g += 1
    elif s_ref_split[i_f][0] == "#":
        i_f += 1
                                                # splitting each line and compare every column
    else:
        mutation_splitted = s_mut_split[i_g].split()
        reference_splitted = s_ref_split[i_f].split()
        if mutation_splitted[0] == reference_splitted[0]:
            if mutation_splitted[1] == reference_splitted[1]:
                if mutation_splitted[3] == reference_splitted[3] and mutation_splitted[4] != reference_splitted[4]:

                    partially_matched.append([mutation_splitted[0], mutation_splitted[1] , mutation_splitted[3], mutation_splitted[4], reference_splitted[4]])

                elif mutation_splitted[3] == reference_splitted[3] and mutation_splitted[4] == reference_splitted[4]:
                    perfectly_matched.append([mutation_splitted[0], mutation_splitted[1], mutation_splitted[3], mutation_splitted[4]])
                i_g += 1
                i_f += 1
            else:
                if int(mutation_splitted[1]) < int(reference_splitted[1]):
                    i_g += 1
                else:
                    i_f += 1
        else:
            if ch_name.index(mutation_splitted[0]) < ch_name.index(reference_splitted[0]):
                i_g += 1
            else:
                i_f += 1
    if s_mut_split[i_g] == '' or s_ref_split[i_f] == '':
        flag = 1
        
#Output
print 'Below is the list of positions where the variants have identical Reference allele but different Alternate allele: '
print "In comparison, %i positions were found to be partially matched" % len(partially_matched)
name1 = ['chromosome #','position','Reference Allele','Alternate Allele in mutations file', 'Alternate Allele in Individual']

print pandas.DataFrame(name1,partially_matched, columns =[''])

print 'Below is the list of positions where the variants are perfectly matched: '
print "In comparison, %i positions were found to be perfectly matched" % len(partially_matched)
name2 = ['chromosome #','position','Reference Allele','Alternate Allele']

print pandas.DataFrame(name2,perfectly_matched, columns =[''])
