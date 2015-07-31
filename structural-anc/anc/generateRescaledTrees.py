#!/usr/bin/env python
__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "130129"
__usage__ = "PLEASE FILL IN THE USAGE STRING"

import sys, string

from ancestorMixer import *
            
def main(argv=None):
    """
    """

    if argv == None:
        argv = sys.argv[1:]

    try:
        phyml_log = argv[0]
        tree_file = argv[1]
    except IndexError:
        err = "Insufficient number of arguments!  Usage:\n\n%s\n\n" % __usage__
        raise AncestralMixerError(err)

    p = PhymlOutput()
   
    # Read fractional population of each class, likelihood, etc. from phyml
    # log file
    p.readPhymlLogFile(phyml_log)
  
    # Write out individual trees rescaled by structural class 
    p.writeRescaledTrees(tree_file)


if __name__ == "__main__":
    main()
