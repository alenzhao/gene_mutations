import os
import numpy as np
import pandas
from operator import itemgetter, attrgetter, methodcaller
from __builtin__ import list
import sys
from sys import exit


def projectrun(data_file, genes_list_file):


    print 'Starting [           ',
    print '\b'*12,
    sys.stdout.flush()


    os.system("grep '^#' " + data_file + " > ./data_comments.vcf")
    os.system("grep -v '^#' " + data_file + " > ./data_body.vcf")
    os.system("grep '^#' " + genes_list_file + " > ./genes_comments.vcf")
    os.system("grep -v '^#' " + genes_list_file + " > ./genes_body.vcf")

    data = open('./data_body.vcf')
    genes_list = open('./genes_body.vcf')
    s_data = data.read()
    s_genes = genes_list.read()
    data.close()
    genes_list.close()

    s_genes_splitted = s_genes.split('\n')
    s_data_splitted = s_data.split('\n')

    genes_array_dominant = []
    genes_values_dominant = []
    genes_array_recessive = []
    genes_values_recessive = []

    for i in range(len(s_genes_splitted)-1):
        line = s_genes_splitted[i].split('\t')
        if line[5] == 'AR':
            genes_array_recessive.append(line[0])
            genes_values_recessive.append([])
        elif line[5] == 'AD':
            genes_array_dominant.append(line[0])
            genes_values_dominant.append([])



    dominant_hash = {x:y for x, y in zip(genes_array_dominant, genes_values_dominant)}
    recessive_hash = {k:v for k, v in zip(genes_array_recessive, genes_values_recessive)}






    r = 0

    for i in range(len(s_data_splitted)-1):

        r += 1
        if r % 1000 == 0:
            print '\b.',
            sys.stdout.flush()
        tab_splitted = s_data_splitted[i].split()        # data splitted tabular
        if tab_splitted[9].split(':')[0] == '1/1':      # Genotype

                ann_field = tab_splitted[7].split('ANN=')[1].split(';')[0].split(',')   #split the ANN field in the annotated file

                for ii in range(len(ann_field)):                                 # Outputs the IDs in the ANN field according to the criteria

                    if (ann_field[ii].split('|')[2] == "HIGH") or ( ann_field[ii].split('|')[2] == "MODERATE"):

                        gene_string = ann_field[ii].split('|')[3]
                        try:
                            p = recessive_hash.keys().index(gene_string)
                            recessive_hash.values()[p].append(i)
                        except ValueError:
                            continue

        elif tab_splitted[9].split(':')[0] == '0/1':

            ann_field = tab_splitted[7].split('ANN=')[1].split(';')[0].split(',')   #split the ANN field in the annotated file
            for ii in range(len(ann_field)):                                 # Outputs the IDs in the ANN field according to the criteria

                    if (ann_field[ii].split('|')[2] == "HIGH") or ( ann_field[ii].split('|')[2] == "MODERATE"):

                        gene_string = ann_field[ii].split('|')[3]
                        try:
                            p = dominant_hash.keys().index(gene_string)
                            dominant_hash.values()[p].append(i)
                        except ValueError:
                            continue



    i_d = 0
    flag = 0
    while flag == 0:

        if dominant_hash.values()[i_d] == []:

            string = dominant_hash.keys()[i_d]
            dominant_hash.pop(string, 0)

        else:
            string2 = dominant_hash.keys()[i_d]
            v = dominant_hash.values()[i_d]
            v = set(v)
            v = list(v)
            dominant_hash[string2] = v
            i_d += 1

        try:
            u = dominant_hash.values()[i_d + 1]

        except IndexError:
            flag = 1


    i_r = 0
    flag = 0

    while flag == 0:
        unique = []
        if recessive_hash.values()[i_r] == []:

            string = recessive_hash.keys()[i_r]
            recessive_hash.pop(string, 0)

        else:
            string2 = recessive_hash.keys()[i_r]
            v = recessive_hash.values()[i_r]
            v = set(v)
            v = list(v)
            recessive_hash[string2] = v

            i_r += 1

        try:
            u = recessive_hash.values()[i_r + 1]

        except IndexError:
            flag = 1


    recessive_index = []

    for i in range(len(recessive_hash)):
        for ii in range(len(recessive_hash.values()[i])):
            recessive_index.append(recessive_hash.values()[i][ii])

    dominant_index = []

    for i in range(len(dominant_hash)):
        for ii in range(len(dominant_hash.values()[i])):
            dominant_index.append(dominant_hash.values()[i][ii])



    out_recessive = []

    unique_rec = set(recessive_index)
    unique_rec = list(unique_rec)
    unique_rec = sorted(unique_rec)
    for i in range(len(unique_rec)):
        out_recessive.append(s_data_splitted[unique_rec[i]])


    np.savetxt('./output_recessive1.vcf' , out_recessive , delimiter = " " , fmt= "%s")
    os.system("cat ./data_comments.vcf ./output_recessive1.vcf > ./output_recessive.vcf")

    out_dominant = []

    unique_dom = set(dominant_index)
    unique_dom = list(unique_dom)
    unique_dom = sorted(unique_dom)
    for i in range(len(unique_dom)):
        out_dominant.append(s_data_splitted[unique_dom[i]])

    np.savetxt('./output_dominant1.vcf' , out_dominant , delimiter = " " , fmt= "%s")
    os.system("cat ./data_comments.vcf ./output_dominant1.vcf > ./output_dominant.vcf")
    compound_het = []
    compound_het_genes = []

    for i in range(len(recessive_hash)):
        if len(recessive_hash.values()[i]) >= 2:
            compound_het_genes.append(recessive_hash.keys()[i])

            for ii in range(len(recessive_hash.values()[i])):

                compound_het.append(s_data_splitted[recessive_hash.values()[i][ii]])

    np.savetxt('./compound.het.data1.vcf' , compound_het , delimiter = " " , fmt= "%s")
    os.system("cat ./data_comments.vcf ./compound.het.data1.vcf > ./compound.het.data.vcf")
    for i in range(len(compound_het_genes)):
        if i == 0:
            os.system("grep -w '^" + compound_het_genes[i] + "' " + genes_list_file + " > ./compound.het.genes1.vcf")
        else:

            os.system("grep -w '^" + compound_het_genes[i] + "' " + genes_list_file + " >> ./compound.het.genes1.vcf")


    os.system("cat ./genes_comments.vcf ./compound.het.genes1.vcf > ./compound.het.genes.vcf")
    os.system("rm ./data_comments.vcf ./genes_comments.vcf ./data_body.vcf ./genes_body.vcf ./output_dominant1.vcf ./output_recessive1.vcf ./compound.het.data1.vcf ./compound.het.genes1.vcf")

    print '\b]  Done!',

if __name__ == '__main__':

    if len(sys.argv)!=3:
        print ("Invalid syntax, You should define two input files, the first file is the annotated vcf file and the second file is your genes list  ")
    else:
        data_file = sys.argv[1]
        genes_list_file = sys.argv[2]
        projectrun(data_file,genes_list_file)