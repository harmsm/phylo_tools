#!/usr/bin/env python3
__description__ = \
"""
createCloseAncestors.py

Take ancestral reconstruction codeml files and minimize the number of sequence 
differences between two reconstructed ancestors by: 
  1) choosing only one of "equivalent" residues (K for R, etc.) that changed
  2) setting ambiguous residues to the same value.  
This should minimize the number of sequence differences between the ancestors,
consistent with the constraints set by phylogenetic uncertainty and physical
chemical reasoning.  
"""
__author__ = "Michael J. Harms, harmsm@gmail.com"
__usage__ = "yo.py anc1_file anc2_file"
__date__ = "2015-06-21"

import sys

def readDatFile(dat_file,take_only_top=20,conv_cys_ser=False):
    """
    Read a codeml dat file.

    take_only_top is an integer that forces the program to only sample the top
    N states.  By default it is set to 20, so all states will be sampled. 
    conv_cys_ser will replace all cys residues with ser.
    """

    f = open(dat_file,'r')
    lines = f.readlines()
    f.close()

    # Remove blank lines and commented (#) lines
    lines = [l for l in lines if l[0] != "#" and l.strip() != ""]

    # Read the amino acids and posterior probabilities for each site into a
    # list.
    aa_list = []
    pp_list = []
    for l in lines:
        aa_list.append([])
        pp_list.append([])
        columns = l.split()[1:]
        for i in range(0,len(columns),2):

            if take_only_top and i > (2*take_only_top):
                break
                
            aa_list[-1].append(columns[i])
            pp_list[-1].append(float(columns[i+1]))
            
        # Normalize the values in pp to be between 0 and 1
        pp_list[-1] = [p/sum(pp_list[-1]) for p in pp_list[-1]]

        # Hacked code that converts any cys residues to ser
        if conv_cys_ser:
            if "C" in aa_list[-1]:
                if "S" in aa_list[-1]:
           
                    c_index = aa_list[-1].index("C")
                    c_pp = pp_list[-1].pop(c_index)
                    c_aa = aa_list[-1].pop(c_index)

                    s_index = aa_list[-1].index("S")
                    pp_list[-1][s_index] += c_pp
 
                else:
                    aa_list[-1][aa_list[-1].index("C")] = "S"


    # Combine amino acid and posterior probability data into a single list
    anc_data = [aa_list,pp_list]

    return anc_data


def fudgeAncestors(anc1_file,anc2_file,pp_cutoff=0.20):
    """
    Given two codeml files, find sequences for the ancestors with the idea of
    minimizing the sequence distance between the two by merging equivalent
    amino acids (e.g. R/K) and looking for ambiguities in the reconstructions.
    anc1 is favored, so it will be closest to its ML state, whereas anc2 will
    be tweaked to be as close to anc1 as possible.

        anc1_file and anc2_files: paml ancestor outputs (should have same 
                                  number of sites)
        pp_cutoff: posterior probability cutoff above which an alternate
                   state can be considered for the fudging

    """

    # read in the codeml fiels
    anc1 = readDatFile(anc1_file,conv_cys_ser=True)
    anc2 = readDatFile(anc2_file,conv_cys_ser=True)

    out_ml = []
    out_seq = []
    for i in range(len(anc1[0])):

        # Grab the ml aa at this site
        anc1_ml = anc1[0][i][0]
        anc2_ml = anc2[0][i][0]

        # record ml
        out_ml.append((anc1_ml,anc2_ml))

        # if the two ml aa are the same, record and go to next site
        if anc1[0][i][0] == anc2[0][i][0]:
            out_seq.append(anc1[0][i][0])
            continue

        # ---------------------------------------------------------------------
        # Look for equivalent amino acids (e.g. K/R, E/D, etc.)
        # ---------------------------------------------------------------------

        ml_set = [anc1_ml,anc2_ml]
        ml_set.sort()
        ml_set = tuple(ml_set)

        equivalent_sets = [("D","E"),
                           ("C","S"),
                           ("K","R"),
                           ("S","T"),
                           ("N","Q")]

        # if we're equivalent sets, record and go to next site
        if ml_set in equivalent_sets:
            out_seq.append(anc1_ml)
            continue 

        # ---------------------------------------------------------------------
        # Look for relatively high posterior probability alternate 
        # reconstructions that would make anc1 and and2 have the same state at
        # this position.
        # ---------------------------------------------------------------------

        anc1_poss = dict(zip(anc1[0][i],anc1[1][i]))
        anc2_poss = dict(zip(anc2[0][i],anc2[1][i]))

        possibilities = [] 
        try:
            if anc1_poss[anc2_ml] > pp_cutoff:
                possibilities.append(anc2_ml)
        except KeyError:
            pass
       
        try: 
            if anc2_poss[anc1_ml] > pp_cutoff:
                possibilities.append(anc1_ml)
        except KeyError:
            pass

        # If we found ambiguous residues, record and move to next amino acid.
        # Favor anc1
        if len(possibilities) > 0:
            if len(possibilities) == 2:
                out_seq.append(anc1_ml)
            else:
                out_seq.append(possibilities[0])
            continue

        # ---------------------------------------------------------------------
        # If we get here, this is real difference we can't fudge away.  Record.
        # ---------------------------------------------------------------------
    
        out_seq.append("{:s}|{:s}".format(anc1_ml,anc2_ml))

    for i in range(len(out_ml)):
        status = out_ml[i][1] == out_seq[i] 
        print(out_ml[i][0],out_ml[i][1],out_seq[i],status)

def main(argv=None):
    """
    """

    if argv == None:
        argv = sys.argv[1:]

    try:
        anc1_file = argv[0] 
        anc2_file = argv[1]
    except IndexError:
        err = "Incorrect arguments. Usage:\n\n{:s}\n\n".format(__usage__)
        raise IndexError(err)

    fudgeAncestors(anc1_file,anc2_file)
    
if __name__ == "__main__":
    main()
