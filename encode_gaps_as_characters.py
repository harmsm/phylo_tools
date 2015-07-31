#!/usr/bin/env python3
__description__ = \
"""
Take an alignment and create an R-readable table in which each site is mapped 
a gap (0) or non-gap (1) for each taxon.  This can then be fed to APE for a 
character reconstruction.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-07-28"
__usage__ = "encode_gaps_as_characters.py fasta_file"

import sys, string

def encode_gaps(fasta_file):
    """
    """
    
    character_map = dict([(char,1) for char in string.ascii_uppercase])
    character_map["-"] = 0   

    taxon = None
    sequence = []

    out_dict = {}
    with open(fasta_file,'r') as f:

        for line in f:

            if line.startswith(">"):
                if taxon:
                    character_array = [character_map[s] for s in sequence]
                    out_dict[taxon] = character_array[:]

                    sequence = []

                taxon = line[1:].strip()
        

            else:
                if taxon:
                    sequence.extend(list(line.strip()))

    if taxon:
        character_array = [character_map[s] for s in sequence]
        out_dict[taxon] = character_array[:]
       
    num_sites = len(list(out_dict.values())[0])

    header = ["{:12s} {:12s}".format("  ","taxon")]
    for i in range(num_sites):
        header.append("{:8s}".format("s{:d}".format(i)))
    
    out = ["".join(header)]
    for i, k in enumerate(list(out_dict.keys())):
        x = "".join(["{:8d}".format(j) for j in out_dict[k]])
        out.append("{:12d} {:12s}{:s}".format(i,k,x))

    return "\n".join(out)


def main(argv=None):
    
    if argv == None:
        argv = sys.argv[1:]

    try:
        fasta_file = argv[0]
    except IndexError:
        err = "incorrect arguments. Usage:\n\n{:s}\n\n".format(__usage__)
        raise IndexError(err)

    return encode_gaps(fasta_file)

if __name__ == "__main__":
    print(main())
