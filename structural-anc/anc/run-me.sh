#!/bin/bash

USAGE="run-me.sh file_root [excecute in directory containing all i/o from phyml-ss run]"
PAMLDAT="$HOME/local/lib/paml4.7/dat/"

# Parse command line
root=${1}
if [[ ! "${root}" ]]; then
    echo ${USAGE} 
    exit
fi

# Use nw utils to strip branch supports out of the phyml output tree (paml will
# choke on these)
nw_topology -b -I ${root}.phy_phyml_tree.txt > ${root}_paml.newick

# Generate trees for each structual class, scaled according to relative mutation
# rate.
./generateRescaledTrees.py ${root}.phy.log ${root}_paml.newick

# Convert the .phy file used by paml to a fasta file, tossing structrual info
# from the bottom
./phy2fasta.py ${root}.phy > ${root}.fasta

# Grab alpha for gamma distribution from log file
alpha=`grep Alpha ${root}.phy.log | tail -1 | awk '{print $8}' | sed 's/\]//'`
echo ${alpha}

# Use lazarus to generate ancestors for each structural class in phyml-ss
for x in `echo "be bo bh ee eo eh"`; do
    echo "Running lazarus on ${x}"
    lazarus_batch.py --alignment ${root}.fasta --tree ${x}_${root}_paml.newick --model $PAMLDAT/exeho_${x}.dat --branch_lengths fixed --asrv 8 --alpha ${alpha} --codeml &> ${x}.log
    mv tree1 class_${x}
done

# Now mix the ancestors spit out by each structural class lazarus run
./ancestorMixer.py ${root}.phy ${root}.phy.log 

cp class_be/tree1.txt final_anc
