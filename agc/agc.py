#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import textwrap
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default=400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default=10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default=100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default=8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file, "rt") as file:
        seq = ""
        for line in file:
            if line.startswith(">"):
                if len(seq) >= minseqlen:
                    yield seq
                seq = ""
                continue
            seq += line.strip()
    yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    seq_uniq = []
    occ = []
    sequences = read_fasta(amplicon_file, minseqlen)
    for sequence in sequences:
        if sequence in seq_uniq:
            index = seq_uniq.index(sequence)
            occ[index] = occ[index]+1
        else:
            seq_uniq.append(sequence)
            occ.append(1)

    zipped = sorted(zip(occ, seq_uniq), reverse=True)
    unique_sorted = []
    occ_sorted = []
    for occ, seq in zipped:
        unique_sorted.append(seq)
        occ_sorted.append(occ)

    for i in range(len(occ_sorted)):
        if occ_sorted[i] > mincount:
            yield [unique_sorted[i], occ_sorted[i]]


def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    compteur = 0
    lon = len(alignment_list[0])
    for i in range(lon):
        if alignment_list[0][i] == alignment_list[1][i]:
            compteur += 1
    return compteur/lon * 100


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount,chunk_size, kmer_size):
    ref = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    tpm = []
    for seq in ref:
        for seq2 in dereplication_fulllength(amplicon_file, minseqlen, mincount):
            align = nw.global_align(seq[0], seq2[0], gap_open=-1, gap_extend=-1,
                                    matrix=os.path.abspath(os.path.join(
                                        os.path.dirname(__file__), "MATCH")))
            if get_identity(align) <= 97:
                tpm.append(seq)
    return tpm


def write_OTU(otu_list, output_file):
    with open(output_file, 'w') as file:
        for i, seq in enumerate(otu_list):
            file.write(
                f">OTU_{i+1} occurrence:{seq[1]}\n{textwrap.fill(seq[0], width=80)}\n")


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    # write_OTU(abundance_greedy_clustering(args.amplicon_file,args.minseqlen, args.mincount, args.chunk_size, args.kmer_size),args.output_file)

# ==============================================================
# Chimera removal section
# ===============================================================


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))




if __name__ == '__main__':
    main()
