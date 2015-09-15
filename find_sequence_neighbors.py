#!/usr/bin/env python3

import itertools

def find_seq_neighbors(seq,max_num_mutations=2,alphabet=("A","T","G","C")):
    """
    Take a sequence and generate a list of all possible sequence neighbors 
    within max_num_mutations.  

        input:

        seq: sequence string (assumes all letters are within alphabet)
        max_num_mutations: integer indicating how many mutations away to walk
        alphabet: possible states at each site.  

        output:
    
        all_possible_neighbors: a list of strings containing all possible 
                                neighbors to seq.
    """
        
    wt_seq = list(seq)
    num_sites = len(wt_seq)
    all_possible_neighbors = []

    # For all possible numbers of mutations (0 through max_num...)
    for i in range(max_num_mutations+1):

        # Create a list of all possible combinations of i states given the
        # alphabet
        state_combos = list(itertools.combinations_with_replacement(alphabet,(i)))

        # Go through all possible "num_sites choose i" combinatons of site indexes
        for sites in itertools.combinations(range(num_sites),i):

            # Now go through possible state combinations for these sites.
            for j in range(len(state_combos)):

                # Do the actual mutations to the string
                mutated_seq = wt_seq[:]
                for k in range(i):
                    mutated_seq[sites[k]] = state_combos[j][k]

                all_possible_neighbors.append("".join(mutated_seq))

    return all_possible_neighbors
       

 
GENCODE = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
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

def translate(sequence):
    """
    Translate a nucleotide sequence into a protein sequence.  If there is a
    problem, write an "X" into the sequence.
    """

    try:
        return "".join([GENCODE[sequence[3*i:3*i+3]]
                        for i in range(len(sequence)//3)])
    except KeyError:
        out = []
        for i in range(len(sequence)//3):
            try:
                out.append(GENCODE[sequence[3*i:3*i+3]])
            except KeyError:
                out.append("X")
        return "".join(out)


seq = "GGTGGAGGTTCGGCCGAA" 
neighbors = find_seq_neighbors(seq,1)
translated = [translate(n) for n in neighbors]
unique_seq = dict([(t,()) for t in translated]).keys()

for u in unique_seq:
    print(u)
