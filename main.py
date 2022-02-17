from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    #parameters for NW
    gap_open=-10
    gap_extend=-1
    nw=NeedlemanWunsch("substitution_matrices/BLOSUM62.mat", gap_open, gap_extend)

    #use a dict to keep track of all the species alignment scores
    score_dict={}
    hs_gg_score, hs_gg_alignA, hs_gg_alignB = nw.align(hs_seq, gg_seq)
    score_dict['Gallus gallus']=hs_gg_score

    hs_mm_score, hs_mm_alignA, hs_mm_alignB = nw.align(hs_seq, mm_seq)
    score_dict['Mus musculus']=hs_mm_score

    hs_br_score, hs_br_alignA, hs_br_alignB = nw.align(hs_seq, br_seq)
    score_dict['Balaeniceps rex']=hs_br_score

    hs_tt_score, hs_tt_alignA, hs_tt_alignB = nw.align(hs_seq, tt_seq)
    score_dict['Tursiops truncatus']=hs_tt_score

    #sort the dictionary
    sorted_score_dict={k: v for k, v in sorted(score_dict.items(), key=lambda item: item[1])}
    
    #print the results
    print('SPECIES ORDERED BY SIMILARITY TO HUMAN BRD2')
    ordered_species=[x for x in sorted_score_dict.keys()][::-1] #need to reverse b/c larger score=more related (on the gene tree)
    for x,i in zip(ordered_species,range(len(ordered_species))):
        print('({}) {}'.format(i+1,x))
    
    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    print('ALIGNMENT SCORE AGAINST HUMAN BRD2')
    for x in ordered_species:
        print('{}: {}'.format(x, sorted_score_dict[x]))

if __name__ == "__main__":
    main()
