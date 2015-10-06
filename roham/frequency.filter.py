import numpy as np
import sys

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print ("Invalid syntax, you are only allowed to specify one input file! ")
        exit()
    else:

        input_file = sys.argv[1]
        print "This script will filter your input VCF file (in ESP 6500 Exome data format)" \
              " based on minor allele frequency:"

        while True:
            try:
                maf = int(raw_input('Please choose your criteria for allele frequency:\n 1) Based on European American'
                                    ' allele count: Enter "0"\n 2) Based on African American allele count: Enter "1"\n '
                                    '3) Based on all of the allele counts: Enter "2"\n :  '))

            except ValueError:
                print "Sorry, I didn't understand that."
                continue

            if not (maf == 0 or maf == 1 or maf == 2):
                print "Sorry, Make sure your input is one of the following integers: 0, 1, or 2."
                continue
            else:
                break
        while True:
            try:
                threshold = float(raw_input('Now please enter your desired Threshold for the allele frequency'
                                            ' (in percentage): '))

            except ValueError:
                print "Sorry, I didn't understand that. Make sure your input is a float"
                continue

            else:
                break


print 'Starting [           ',
print '\b'*12,
sys.stdout.flush()


data_open = open(input_file)
data_read = data_open.read()

data_open.close()

data = data_read.split('\n')
output = []
flag = 0
i_d = 0
r = 0
while flag == 0:

    if r%10000 == 0:
        print '\b.',
        sys.stdout.flush()
    r += 1

    if data[i_d][0] == '#':
        output.append(data[i_d])
        i_d += 1
        continue

    line_splitted = data[i_d].split()
    info_column = line_splitted[7]
    MAF_column = info_column.split(';')[4]
    freq = MAF_column.split(',')[maf]
    # print freq
    if float(freq) < threshold:

        output.append(data[i_d])
    i_d += 1

    if data[i_d] == '':
        flag = 1


np.savetxt('./frequency_filter_output.vcf', output, delimiter = " ", fmt= "%s")

print '\b]  Done!',
