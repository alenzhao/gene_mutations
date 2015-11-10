from codontable import gencode
from read_muttable import get_mutation_rates

# transcript = 'CAGTCATCGAGGATAGTAGTC'

def possible_mutations(transcript):

    mut_table = get_mutation_rates()



    missense_mut_prob = []
    nonsense_mut_prob = []
    silent_mut_prob = []
    donor_splice_mut_prob = []
    acceptor_splice_mut_prob = []


    bases = ['A','C','G','T']
    if transcript[1:3] == 'AG':
        trimer1 = transcript[0:3]
        for i in 'CGT':
            prob = mut_table[(trimer1, i)]
            acceptor_splice_mut_prob.append(prob)

        trimer2 = transcript[1:4]
        for i in 'ACT':
            prob = mut_table[(trimer2, i)]
            acceptor_splice_mut_prob.append(prob)
    else:
        acceptor_splice_mut_prob.append('.')

    #print acceptor_splice_mut_prob
    if transcript[-3:-1] == 'GT':
        trimer1 = transcript[-3:]
        for i in 'ACG':
            prob = mut_table[(trimer1, i)]
            donor_splice_mut_prob.append(prob)

        trimer2 = transcript[-4:-1]
        for i in 'ACT':
            prob = mut_table[(trimer2, i)]
            donor_splice_mut_prob.append(prob)
    else:
        donor_splice_mut_prob.append('.')

    #print donor_splice_mut_prob
    for i in range(3,len(transcript)-3,3):

        codon = transcript[i:i+3]
        prev = transcript[i-1]
        nxt = transcript[i+3]
        for ind in xrange(4):
            if codon[0] != bases[ind]:
                mutated = bases[ind] + codon[1:3]
                trimer = prev + codon[0:2]
                prob = mut_table[(trimer, bases[ind])]
                if gencode[mutated] == gencode[codon]:

                    silent_mut_prob.append(prob)

                elif gencode[mutated] == '*':

                    nonsense_mut_prob.append(prob)
                else:
                    missense_mut_prob.append(prob)

            if codon[1] != bases[ind]:
                mutated = codon[0] + bases[ind] + codon[2]
                trimer = codon
                prob = mut_table[(trimer, bases[ind])]
                if gencode[mutated] == gencode[codon]:

                    silent_mut_prob.append(prob)

                elif gencode[mutated] == '*':

                    nonsense_mut_prob.append(prob)
                else:
                    missense_mut_prob.append(prob)
            if codon[2] != bases[ind]:
                mutated = codon[0:2] + bases[ind]
                trimer = codon[1:3] + nxt
                prob = mut_table[(trimer, bases[ind])]
                if gencode[mutated] == gencode[codon]:

                    silent_mut_prob.append(prob)

                elif gencode[mutated] == '*':

                    nonsense_mut_prob.append(prob)
                else:
                    missense_mut_prob.append(prob)

    # print len(missense_mut_prob) + len(silent_mut_prob) + len(nonsense_mut_prob) + len(donor_splice_mut_prob) + len(acceptor_splice_mut_prob)
    sum_missense_probs = sum(missense_mut_prob)
    sum_nonsense_probs = sum(nonsense_mut_prob)
    sum_silent_probs = sum(silent_mut_prob)
    sum_accsplice_probs = sum(acceptor_splice_mut_prob)
    sum_donsplice_probs = sum(donor_splice_mut_prob)
    return [sum_accsplice_probs, sum_donsplice_probs, sum_missense_probs, sum_nonsense_probs, sum_silent_probs]

