#!/usr/bin/env python
__description__ = \
"""
Converts the .phy file used as input for a phyml-ss calculation (with its 
stockholm-style structural annotations) and converts it into a simple, 
paml-friendly fasta file with no structural information.
"""
__author__ = "Michael J. Harms"
__date__ = "130426"
__usage__ = "phy2fasta.py phy_file > fasta_output_file"

import sys

f = open(sys.argv[1],'r')
lines = f.readlines()
f.close()

lines = [l for l in lines[2:] if l.strip() != "" and not l.startswith("#")]

for l in lines:
    c = l.split()
    print ">%s\n%s" % tuple(c)

