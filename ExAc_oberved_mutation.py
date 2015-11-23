import subprocess

def read_vcf(vcffile,chrom,exon,transcript_IDs):


    tabix_command = ["/Users/rohamrazaghi/htslib/tabix"]
    if vcffile == 'None': tabix_command.append("/Users/rohamrazaghi/Documents/data/ExAC.r0.3.sites.vep.vcf.gz")
    else: tabix_command.append(vcffile)

    tabix_command.append(chrom + ':' + exon[0] + '-' + exon[1])

    outfilename = 'EXAC.data.' + chrom + ':' + exon[0] + '-' + exon[1]
    f = open(outfilename,'wb')
    subprocess.call(tabix_command,stdout=f)
    f.close()
    silent_allele_count = 0
    missense_allele_count = 0
    stopgain_allele_count = 0
    splice_acceptor_allele_count = 0
    splice_donor_allele_count = 0

    File= open(outfilename,'r')
    for line in File:
        if len(line) == 0:
            break

        var = line.strip().split('\t')
        chrom = var[0]
        position = var[1]
        ref = var[3]
        alleles = var[4].split(',')
        info = var[7].split(';')

        if info[0].split('=')[0] == 'AC':
            AC = info[0].split('=')[1].split(',')
            for pr in range(len(AC)):
                AC[pr] = int(AC[pr])
        else:
            continue

        try:
            csq = var[7].split('CSQ=')[1].split(',')
        except IndexError:
            continue
        hash_t = {}
        for ind in range(len(alleles)):
            hash_t[alleles[ind]] = AC[ind]
        for i in range(len(csq)):
            vep = csq[i].split('|')
            allele = vep[0]
            trans_id = vep[2]
            prediction = vep[4]
            try:
                ac = hash_t[allele]
            except KeyError:
                continue
            for ii in range(len(transcript_IDs)):
                if trans_id == transcript_IDs[ii]:
                    if 'splice_acceptor_variant' in prediction:
                        splice_acceptor_allele_count += ac
                    if 'splice_donor_variant' in prediction:
                        splice_donor_allele_count += ac
                    if 'stop_gained' in prediction:
                        stopgain_allele_count += ac
                    if 'missense_variant' in prediction:
                        missense_allele_count += ac
                    if 'synonymous_variant' in prediction:
                        silent_allele_count += ac

    return [splice_acceptor_allele_count,splice_donor_allele_count,missense_allele_count,stopgain_allele_count,silent_allele_count]



# vcffile = 'None'
# transcriptids = ['ENST00000373443','ENST00000294517','ENST00000481886','ENST00000398167']
# exon = ['33546713','33586132']
# chrom = '1'
#
# [a,b,c,d,e] = read_vcf(vcffile,chrom,exon,transcriptids)
# print a
# print b
# print c
# print d
# print e
