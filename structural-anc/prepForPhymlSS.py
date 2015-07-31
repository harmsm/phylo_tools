#!/usr/bin/env python3
__description__ = \
"""
Take a fasta file with an alignment, grab the secondary structure and solvent
accessibility for each site, then create a stockholm file encoding secondary
structure call and accessiblity for input to phyml-ss.
"""
__author__ = "Michael J. Harms"
__date__ = "121101"
__usage__ = "prepForPhymlSS.py fasta_file calc_dir struct_alignment_dir"

import sys, os
from math import sqrt, ceil

aa_3to1 = {"ALA":"A",
           "CYS":"C",
           "ASP":"D",
           "GLU":"E",
           "PHE":"F",
           "GLY":"G",
           "HIS":"H",
           "ILE":"I",
           "LYS":"K",
           "LEU":"L",
           "MET":"M",
           "ASN":"N",
           "PRO":"P",
           "GLN":"Q",
           "ARG":"R",
           "SER":"S",
           "THR":"T",
           "VAL":"V",
           "TRP":"W",
           "TYR":"Y"}


class PhymlSSPrepError(Exception):
    """
    General error class for this module.
    """

    pass


def readFastaFile(fasta_file):
    """
    Read a fasta file into a dictionary keying sequence name to sequence.  All
    sequences must have the same length.
    """

    f = open(fasta_file)

    out_dict = {}
    fasta_buffer = []
    recording = False
    for l in f.readlines():

        # New record
        if l.startswith(">"):
            if recording:
                seq = "".join([b.strip() for b in fasta_buffer])
                out_dict[key] = seq
            else:
                recording = True

            # Strip out "struct_" annotation
            key = l[1:].strip()
            key = key.replace("struct_","")
               
            # Make sure the sequence name is unique 
            if out_dict.has_key(key):
                err = "Sequence \"{:s}\" was seen more than once!\n" % seq
                raise PhymlSSPrepError(err)

            fasta_buffer = []
    
        # should we record?
        else:
            if recording:
                fasta_buffer.append(l)

    if recording:
        seq = "".join([b.strip() for b in fasta_buffer])
        out_dict[key] = seq

    f.close()

    # make sure all sequences have the same length  
    current_length = -1
    for k in out_dict.keys():
        seq_length = len(out_dict[k])
        if current_length == -1:
            current_length = seq_length

        if seq_length != current_length:
            err = "Not all sequences have the same length! ({:s})\n" % k
            raise PhymlSSPrepError(err)

    return out_dict

def parseCalcOutput(calc_dir,seq_dict):
    """
    """

    current = os.getcwd()
    os.chdir(calc_dir) 

    out_dict = {}
    for k in seq_dict.keys():

        f = open(os.path.join(k,"{:s}.summary" % k),'r')
        lines = f.readlines()
        f.close()

        aa_out = [None for i in range(len(seq_dict[k]))]
        sa_out = [None for i in range(len(seq_dict[k]))]
        ss_out = [None for i in range(len(seq_dict[k]))]

        start = False
        current_aa_index = 0 
        for l in lines:

            col = l.split()
            aa = col[1]
            sa = float(col[2])
            ss = col[3].strip("\"")
            
            # We've reached the last amino acid in the alignment, but not
            # the last in the file. 
            if current_aa_index == len(seq_dict[k]):
                break

            # Skip over gaps 
            while current_aa_index < (len(seq_dict[k])-1) and \
                seq_dict[k][current_aa_index] == "-":
                current_aa_index += 1
        
            if aa != seq_dict[k][current_aa_index]:
                continue

            if seq_dict[k][current_aa_index] != aa:
 
                err = "Mismatch between sequence in pdb file and fasta file!\n"
                err += "Entry/pdb {:s}, column %i\n" % (k,current_aa_index)
                raise PhymlSSPrepError(err)

            aa_out[current_aa_index] = aa 
            sa_out[current_aa_index] = sa
            ss_out[current_aa_index] = ss
    
            current_aa_index += 1

        out_dict[k] = [aa_out,sa_out,ss_out]
            

    os.chdir(current)
 
    return out_dict

   

def yo(fasta_file,calc_dir,struct_align_dir,classification_dir):
    """
    """
    
    # Read fasta file
    seq_dict = readFastaFile(fasta_file)

    # Read calculated parameters for each site, aligned to fasta file
    calc_output = parseCalcOutput(calc_dir,seq_dict)

    num_align_columns = len(calc_output.values()[0][0])

    ### SOMEHOW CREATE LIST OF PDB FILES ### pdb_files

    stockholm_ss = []
    stockholm_sa = []
   
    # For every column in the alignment... 
    for i in range(num_align_columns):

        ss = []
        sa = []

        # Go through every pdb file...
        for p in pdb_files:

            # Read the calculated parameters for this site from that pdb
            if calc_output[p][3][i] != None:
                sa.append(calc_output[p][1][i])
                ss.append(calc_output[p][2][i])

        # If there is any structural information at all
        if len(sa) > 0:

            # calculate the mean and standard deviation of the solvent 
            # accessibility
            sa_mean = sum(sa)/(1.0*len(sa))
            sa_sd = [(s - sa_mean)**2 for s in sa]
            sa_sd = sqrt(sum(sa_sd)/(1.0*len(sa_sd)))

            # Record the most common secondary structure call.  If there is
            # no "most common" record it as ?
            final_ss = "?"
            observed_ss = dict((s,[]) for s in ss).keys()
            if len(observed_ss) == 1:
                final_ss = "{:s}"  % ss[0]
            else:

                # See how often each ss type is observed, recording it as
                # a list like [(ss1_counts,ss1),(ss2_counts,ss2)...]
                counts = [(len([t for t in ss if t == s]),s)
                          for s in observed_ss]

                # Sort the list 
                counts.sort(reverse=True)
    
                # If one ss element is the most common, record it
                if counts[0][0] > counts[1][0]:
                    final_ss = "{:s}" % counts[0][1]
                else:
                    final_ss = "?"

            # Call un-identified guys coil 
            if final_ss == "":
                final_ss = "C"

            stockholm_ss.append(final_ss)
            if sa_mean == 0.0:
                final_sa = 0
            elif sa_mean > 1.0:
                final_sa = 1.0
            else:
                final_sa = int(ceil(sa_mean*10))-1

            stockholm_sa.append("%i" % final_sa)

        # If there is no structural information, record "." (ambiguous)
        else:
            stockholm_sa.append(".")
            stockholm_ss.append(".")

    print("# {:s}".format(c))
    print("#=GR SS_cons\n{:s}".format("".join(stockholm_ss)))
    print("#=GR SA_cons\n{:s}".format("".join(stockholm_sa)))

    return [""]


def main(argv=None):
    """
    Main function.
    """
    
    if argv == None:
        argv = sys.argv[1:]

    try:
        fasta_file = argv[0]
        calc_dir = argv[1]
    except IndexError:
        err = "Insufficient number of arguments!\n\nUsage:\n\n{:s}\n\n" % __usage__
        raise PhymlSSPrepError(err)

    yo(fasta_file,calc_dir)

if __name__ == "__main__":
    main()
    
