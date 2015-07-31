#!/usr/bin/env python
__description__ = \
"""
Take the contents of a directory containing multiple .dssp files, extract the 
residue-by-residue relative solvent accessibility and secondary structure 
assignment, calculate the average (across chains and models), and write out.
"""
__author__ = "Michael J. Harms"
__date__ = "121014"
__usage__ = "parseDSSP.py"

import os, string

SA_DICT = {"A":115.0,
           "R":225.0,
           "D":150.0,
           "N":160.0,
           "C":135.0,
           "E":190.0,
           "Q":180.0,
           "G": 75.0,
           "H":195.0,
           "I":175.0,
           "L":170.0,
           "K":200.0,
           "M":185.0,
           "F":210.0,
           "P":145.0,
           "S":115.0,
           "T":140.0,
           "W":255.0,
           "Y":230.0,
           "V":155.0,
           "c":190.0}   # Cysteine in disulfide


def readFile(input_file,chains_to_take=("A","B"),take_ss=("E","H")):
    """
    Read data from a DSSP file.
    """

    f = open(input_file,'r')
    lines = f.readlines()
    f.close()

    for i, l in enumerate(lines):
        if l.startswith("  #  RESIDUE"):
            start_index = i + 1
            break
    
    lines = lines[start_index:]
    out_dict = {}
    for l in lines:
        aacid = l[13]
        if aacid == "!":
            continue

        # A lowercase letter denotes a disulfide bond--standardeize to a 
        # lowercase c.
        if aacid in string.letters[:26]:
            aacid = "c"

        resid = int(l[5:10])

        # Only grab correct chains
        chain = l[11]
        if chain not in chains_to_take:
            continue

        # Key ss and asa to the residue/aacid pair
        ss    = l[16]
        if ss not in take_ss:
            ss = "C"

        asa   = float(l[34:38])
        key = (resid,aacid)
        if not out_dict.has_key(key):
            out_dict[key] = [(ss,asa)]
        else:
            out_dict[key].append((ss,asa))

    # Generate pretty output by going through residues in order
    out = []    
    out_keys = out_dict.keys()
    out_keys.sort()
    for k in out_keys:
        resid = k[0]
        aacid = k[1]

        ss = [v[0] for v in out_dict[k]]
    
        # If we have two chains that agree, take it; otherwise,
        # report that we don't know the ss
        final_ss = "?"
        if len(ss) == 2:
            if ss[0] == ss[1]:
                ss = ss[0]
            else:
                ss = "?"
        else:

            observed_ss = dict((s,[]) for s in ss).keys()

            if len(observed_ss) == 1:
                final_ss = observed_ss[0]
                continue

            # See how often each ss type is observed, recording it as
            # a list like [(ss1_counts,ss1),(ss2_counts,ss2)...]
            counts = [(len([t for t in ss if t == s]),s)
                      for s in observed_ss]

            # Sort the list 
            counts.sort(reverse=True)

            # If one ss element is the most common, record it
            if counts[0][0] > counts[1][0]:
                final_ss = ("%s" % counts[0][1])
            else:
                final_ss  = "?"

           
        # Calculate relative surface area 
        asa = [v[1] for v in out_dict[k]]   
        asa = sum(asa)/len(asa)
        asa = asa/SA_DICT[aacid]

        # Replace any 'c' with "C" for disulfide bonded cys.
        if aacid == "c": 
            aacid = "C" 

        out.append((int(resid),aacid,ss,asa))

    return out


def main(argv=None):
    """
    """

    file_list = [f for f in os.listdir(".") if f[-5:] == ".dssp"]
    out = []
    for f in file_list:
        out.append(readFile(f))

    for i in range(len(out[0])):
        for o in out[1:]:
            if o[i][:2] != out[0][i][:2]:
                err = "Not all .dssp files had the same amino acis in them!\n"
                raise IOError(err)

    ss = [[] for i in range(len(out[0]))]
    asa = [0. for i in range(len(out[0]))]
    for i in range(len(out[0])):
        for o in out:
            ss[i].append(o[i][2])
            asa[i] += o[i][3]

    for i in range(len(out[0])):
        asa[i] = asa[i]/len(out)

    # Find the most common ss call at each site by simple majority.  If there is
    # a tie append "?".  Otherwise, append "O"
    final_ss = []
    for i in range(len(ss)):
                
        observed_ss = dict([(s,[]) for s in ss[i]]).keys()

        if len(observed_ss) == 1:
            final_ss.append(observed_ss[0])
            continue

        # See how often each ss type is observed, recording it as
        # a list like [(ss1_counts,ss1),(ss2_counts,ss2)...]
        counts = [(len([t for t in ss if t == s]),s)
                  for s in observed_ss]

        # Sort the list 
        counts.sort(reverse=True)

        # If one ss element is the most common, record it
        if counts[0][0] > counts[1][0]:
            final_ss.append("%s" % counts[0][1])
        else:
            final_ss.append("?")

    for i in range(len(out[0])):
        print "%10s%10s%10.3f   \"%s\"" % (out[0][i][0],
                                           out[0][i][1],
                                           asa[i], final_ss[i])



if __name__ == "__main__":
    main() 
