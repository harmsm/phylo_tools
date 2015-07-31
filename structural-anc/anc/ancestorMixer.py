#!/usr/bin/env python
__description__ = \
"""
Takes the output from multiple paml ancestral reconstruction calculations and
mixes them accordind to the EX_EHO phyml OPT model.
"""
__author__ = "Michael J. Harms"
__date__ = "130129"
__usage__ = "PLEASE FILL IN THE USAGE STRING"

import sys, string, os, copy

class AncestralMixerError(Exception):
    """
    A general error class for this module.
    """
    
    pass

class PhymlOutput:
    """
    Class for reading phyml-ss output.
    """

    def __init__(self):
        """
        """
      
        # Organize internal matrices by amino acid... 
        self.matrix2aa = ["A","C","D","E","F","G","H","I","K","L","M","N",
                          "P","Q","R","S","T","V","W","Y"]
        self.aa2matrix = dict(zip(self.matrix2aa,range(20)))

        # Stuff read from log file 
        self.log_lk = None
        self.tree_scaler = {}
        self.class_fx = {}
        self.conf = None
        self.final_rate = None      

        # Stuff read from ancestors reconstructed on different class trees
        self.ancestors = {}
        self.node_list = None
        self.num_sites = 0

        # Stuff read from initial .phy file
        self.site_classes = None               
 

    def readAllAncestors(self,prefix="class_"):
        """
        Read reconstructed ancestors for each ex/eho model.
        """

        if self.class_fx == None:
            err = "No phyml log file has been read."
            raise AncestralMixerError(err)

        for c in self.class_fx.keys():

            current_dir = os.getcwd()
            os.chdir("%s%s" % (prefix,c))
            nodes = [int(f[4:-4]) for f in os.listdir(".") if f[-4:] == ".dat"]
            nodes.sort()  

            # Make sure the nodes for this tree are the same as the nodes in any
            # other loaded trees
            if self.node_list == None:
                self.node_list = nodes[:]
            else:
                if self.node_list != nodes:
                    err = "Different categories have different ancestral nodes."
                    raise AncestralMixerError(err)

            self.ancestors[c] = {}

            for n in nodes:

                # Read node .dat file
                f = open("node%i.dat" % n)
                node_lines = f.readlines()
                f.close()

                # Make sure this node has the same number of sites as other nodes
                if self.num_sites == 0:
                    self.num_sites = len(node_lines)
                else:
                    if self.num_sites != len(node_lines):
                        err = "Different nodes have different sequence lengths."
                        raise AncestralMixerError(err)

                m = [[0. for j in range(len(self.matrix2aa))]
                     for i in range(len(node_lines))]

                for i, l in enumerate(node_lines):
                    col = l.split()[1:]

                    for j in range(0,len(col),2):

                        aa = col[j]
                        pp = float(col[j+1])
                       
                        position = self.aa2matrix[aa]


                        m[i][position] = pp

                self.ancestors[c][n] = copy.deepcopy(m)

            os.chdir(current_dir)


    def printAncestors(self):
      
        for n in self.node_list:
            for i in range(self.num_sites):
        
                # Some python magic (which is sadly incomprehensible on first 
                # glance).  Basically makes a tuple of each classes pp for a given
                # aa at a given site.
                print zip(*(self.ancestors[k][n][i] for k in self.ancestors.keys())) 
 
        
    def readPhymlLogFile(self,phyml_file):
        """
        Grab rate scaling for each model from phyml file.
        """

        # Read lines from log file
        f = open(phyml_file,'r')
        lines = f.readlines()
        f.close()

        # Go through log file, parsing ... 
        for l in lines:
            if l.startswith(". Final log"):
                self.log_lk = float(l.split()[5])
                continue

            if l.startswith("E-bur"):
                self.tree_scaler["be"] = float(l.split()[1])
                self.class_fx["be"] = float(l.split()[2])
                continue
            if l.startswith("H-bur"):
                self.tree_scaler["bh"] = float(l.split()[1])
                self.class_fx["bh"] = float(l.split()[2])
                continue
            if l.startswith("O-bur"):
                self.tree_scaler["bo"] = float(l.split()[1])
                self.class_fx["bo"] = float(l.split()[2])
                continue
            if l.startswith("E-exp"):
                self.tree_scaler["ee"] = float(l.split()[1])
                self.class_fx["ee"] = float(l.split()[2])
                continue
            if l.startswith("H-exp"):
                self.tree_scaler["eh"] = float(l.split()[1])
                self.class_fx["eh"] = float(l.split()[2])
                continue
            if l.startswith("O-exp"):
                self.tree_scaler["eo"] = float(l.split()[1])
                self.class_fx["eo"] = float(l.split()[2])
                continue

            if l.startswith("Conf ="):
                self.conf = float(l.split()[2][:-1])
                continue

            if l.startswith("final rate =="):
                self.final_rate = float(l.split()[3])
                continue

        # Some quick error checking
        if None in [self.log_lk,self.conf,self.final_rate]:
            err = "Mangled phyml log file!\n"
            raise AncestralMixerError(err)

        class_check = sum([self.class_fx[k] for k in self.class_fx.keys()]) 
        if class_check < 0.99 or class_check > 1.01:
            err = "Class fractions don't add up to 1"
            raise AncestorMixerError(err)

        # Normalize tree scalers to final rate
        for k in self.tree_scaler.keys():
            self.tree_scaler[k] = self.tree_scaler[k]/self.final_rate
        

    def readClassesFromPhyFile(self,phy_file):
        """
        Read the EX_EHO classes from a .phy file.
        """

        f = open(phy_file,'r')
        lines = [l for l in f.readlines() if l.startswith("#=GR")]
        f.close()

        # Pull secondary structure and solvent accessibilities from the file, doing
        # some sanity checking along the way.

        mangled_file = False
        if len(lines) < 2:
            mangled_file = True

        if not mangled_file:
            ss = [l.split()[2] for l in lines if l.split()[1] == "SS_cons"]
            sa = [l.split()[2] for l in lines if l.split()[1] == "SA_cons"]

        if len(ss) != 1 or len(sa) != 1:
            mangled_file = True
        else:
            ss = ss[0]
            sa = sa[0]
            if len(ss) != len(sa):
                mangled_file = True

        if mangled_file:
            err = "Mangled phy file.  Must have #=GR entries for secondary "
            err += "structure (SS_cons)\n and solvent accessible (SA_cons)\n"

            raise AncestralMixerError(err)
   
        self.site_classes = []
        for i in range(len(ss)):

            # Find secondary structure
            if ss[i] == "H":
                sec_struct = "h"
            elif ss[i] == "E":
                sec_struct = "e"
            else:
                sec_struct = "o"

            # Find solvent accessibility
            if sa[i] in ["0","."]:
                exposed = "b"
            else:
                exposed = "e"

            # Record structural class for this site
            self.site_classes.append("%s%s" % (exposed,sec_struct))



    def createAmbiguousMixture(self):
   
        # make sure we've actually read a log file 
        if len(self.class_fx) == 0:
            err = "No classes fractions read from phyml log file!"
            raise AncestralMixerError(err)

        # Quick sanity check to make sure classes from log file and asr trees
        # are the same.
        class_keys = self.class_fx.keys()
        class_keys.sort()
        ancestor_keys = self.ancestors.keys()
        ancestor_keys.sort()

        if class_keys != ancestor_keys:
            err = "Structural classes read from phyml log file and from asr"
            err += " trees do no match."
            raise AncestralMixerError(err)
   

        # Go through each node...
        self.ancestors["mix"] = {}
        for n in self.node_list:

            m = [[0. for j in range(len(self.matrix2aa))]
                 for i in range(self.num_sites)]

            # Go through each site in that node...
            for i in range(self.num_sites):

                # Go through each possible state...
                for j in range(len(self.aa2matrix)):
        
                    pp = 0.
                    for c in class_keys:
                        pp += self.ancestors[c][n][i][j]*self.class_fx[c] 

                    m[i][j] = pp
                   
            self.ancestors["mix"][n] = copy.deepcopy(m) 

    def createFinalAncestors(self):
        """
        """

        # Do some error checking
        if self.site_classes == None:
            err = "No site classes read from original .phy file"
            raise AncestralMixerError(err)

        if len(self.site_classes) != self.num_sites:
            err = "Number of sites read from original .phy file and number of"
            err += "\nsites read from class asr trees differ."
            raise AncestralMixerError(err)
  
        self.ancestors["final"] = {}      
        for n in self.node_list:

            m = [[0. for j in range(len(self.matrix2aa))]
                 for i in range(self.num_sites)]
             
            for i in range(self.num_sites):

                site_class = self.site_classes[i]         
                for j in range(len(self.aa2matrix)):

                    pp = self.ancestors[site_class][n][i][j]*(1-self.final_rate)
                    pp += self.ancestors["mix"][n][i][j]*self.final_rate
                    m[i][j] = pp
                   
            self.ancestors["final"][n] = copy.deepcopy(m) 

    def writeFinalAncestors(self,output_dir="final_anc"):
        """
        """

        current_dir = os.getcwd()
   
        # Go into the output directory, making it if necessary 
        if not os.path.exists(output_dir):
            tree_depth = os.path.split(output_dir)
            for t in tree_depth:
                if t == "":
                    continue
    
                os.mkdir(t)
                os.chdir(t)
        else:
            os.chdir(output_dir)

        for n in self.node_list:

            f = open("node%i.dat" % n,"w")

            for i in range(self.num_sites):

                pp_list = zip(self.ancestors["final"][n][i],range(len(self.matrix2aa)))
                pp_list.sort(reverse=True)

                f.write("%i " % (i + 1))
                for pp in pp_list:
                    if pp[0] > 0:
                        f.write(" %s %.3f" % (self.matrix2aa[pp[1]],pp[0]))
                    else:
                        break
                f.write("\n")

            f.close()

        os.chdir(current_dir)


    def writeRescaledTrees(self,tree_file):
        """
        Take a tree file and write out tree files for each structual class,
        scaled according to scaling factors in phyml log file.
        """

        for c in self.class_fx.keys():

            scaling_ratio = self.tree_scaler[c]

            # Read the tree file
            f = open(tree_file,'r')
            data = f.read()
            f.close()

            # Do all of the rescaling business
            entries = data.split(":")
            for i in range(len(entries)):
                if entries[i][-1] == ",":
                    v = float(entries[i][-1])*scaling_ratio
                    entries[i] = "%f," % (v)
                else:
                    if entries[i][0] in string.digits:
                        col = entries[i].split(",")
                        if len(col) == 2:
                            v = float(col[0])*scaling_ratio
                            entries[i] = "%f,%s" % (v,col[1])
                        else:
                            col = entries[i].split(")")
                            v = float(col[0])*scaling_ratio
                            entries[i] = "%f)%s" % (v,col[1])

            # Write output
            f = open("%s_%s" % (c,tree_file),"w")
            f.write(":".join(entries))
            f.close()


def main(argv=None):
    """
    """

    if argv == None:
        argv = sys.argv[1:]

    try:
        phy_file = argv[0]
        phyml_log = argv[1]
    except IndexError:
        err = "Insufficient number of arguments!  Usage:\n\n%s\n\n" % __usage__
        raise AncestralMixerError(err)

    # Optional argument
    try:
        output_dir = argv[2]
    except IndexError:
        output_dir = "final_anc"
        
    
    p = PhymlOutput()
   
    # Read class of each site from initial phy file
    p.readClassesFromPhyFile(phy_file)

    # Read fractional population of each class, likelihood, etc. from phyml
    # log file
    p.readPhymlLogFile(phyml_log)

    # Read ancestors correspondin to each structural class
    p.readAllAncestors()

    # Create an ambiguous mixture -- the mixture of pp from each class ancestor
    # weighted by the frequency of that class
    p.createAmbiguousMixture()

    p.createFinalAncestors()

    p.writeFinalAncestors(output_dir)


if __name__ == "__main__":
    main()
