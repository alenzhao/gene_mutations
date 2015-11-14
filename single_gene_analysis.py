from read_fasta_bychrom import read_chromosome, make_fasta_index
from possible_mutations import possible_mutations
from codontable import rev_seq


def single_gene_analysis(gene_name):

    fastafile = '/Users/rohamrazaghi/Documents/data/chromFa/fasta.hg19.all.fa'
    offset_index = make_fasta_index(fastafile)
    transcript_table = {}
    refseqfile = open('/Users/rohamrazaghi/gene_mutations/reftest')
    # gene_name = 'SGIP1'
    for line in refseqfile:
        if line[0] == '#':
            continue
        tx = line.split()
        if tx[12] == gene_name:

            exonstarts = tx[9].split(',')
            del exonstarts[-1]
            exonstarts = map(int, exonstarts)
            exonends = tx[10].split(',')
            del exonends[-1]
            exonends = map(int, exonends)
            # print exonends
            # print exonstarts
            exons = []
            for i in range(len(exonstarts)):
                exons.append([(exonstarts[i]-3), (exonends[i]+3)])
            transcript_table[(tx[12], tx[1])] = [tx[2], tx[3], exons]
    refseqfile.close()
    for i in range(len(transcript_table)):
        trans_id = transcript_table.keys()[i][1]
        chromosome = transcript_table.values()[i][0]
        strand = transcript_table.values()[i][1]
        exons_list = transcript_table.values()[i][2]
        # print trans_id,strand, exons_list
        # fastafile = '/Users/rohamrazaghi/Documents/data/chromFa/' + chromosome + '.fa'
        # offset_index = make_fasta_index(fastafile)
        [sequence, length] = read_chromosome(fastafile,offset_index,chromosome)
        exons_probs = []
        sum_acceptor = 0; sum_donor = 0; sum_missense = 0; sum_nonsense = 0; sum_silent = 0;
        if strand == '+':
            for ind in range(len(exons_list)):
                transcript = sequence[exons_list[ind][0]:exons_list[ind][1]]
                transcript = transcript.upper()
                # print transcript
                [acceptor, donor, missense, nonsense, silent] = possible_mutations(transcript)
                exons_probs.append([acceptor, donor, missense, nonsense, silent])
                sum_acceptor += acceptor; sum_donor += donor; sum_missense += missense; sum_nonsense += nonsense; sum_silent += silent;
        elif strand == '-':
            for ind in range(len(exons_list)):

                transcript = rev_seq(sequence[exons_list[ind][0]:exons_list[ind][1]])
                [acceptor, donor, missense, nonsense, silent] = possible_mutations(transcript)
                exons_probs.append([acceptor, donor, missense, nonsense, silent])
                sum_acceptor += acceptor; sum_donor += donor; sum_missense += missense; sum_nonsense += nonsense; sum_silent += silent;

        transcript_table.values()[i].append(exons_probs)
        transcript_table.values()[i].append(sum_acceptor)
        transcript_table.values()[i].append(sum_donor)
        transcript_table.values()[i].append(sum_missense)
        transcript_table.values()[i].append(sum_nonsense)
        transcript_table.values()[i].append(sum_silent)

    return transcript_table

# gene_name = 'SGIP1'
# test = get_transcript(gene_name)
#
# print test