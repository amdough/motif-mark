import math
import re
import argparse
import cairo
import bioinfo 

def get_args():
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument("-f", help="Designates absolute file path to input fasta file", type=str, required =True)
    parser.add_argument("-m", help="Designates absolute file path to input motifs file", type=str, required=True)
    parser.add_argument("-o", help="Designates output file for motif-mark stats", type=str, required=True)

    return parser.parse_args()


def __main__():

    args = get_args()
    input_fasta = args.f 
    input_motifs = args.m 
    output = args.o 

    class Sequence:
        def __init__(self, seq):
            '''This class represents an input fasta sequence.'''
            self.seq = seq
    
    class Motif: 
        def __init__(self, motif):
            '''This class represents an input motif.'''
            self.motif = motif

        def find_motif(self, Sequence):
            '''This method finds the motif in the sequence and returns the number of times it is found.'''
            return len(re.findall(self.motif, Sequence.seq))

    class Feature:
        def __init__(self, name, start, end):
            '''This class represents a feature in the sequence.'''
            self.name = name
            self.start = start
            self.end = end






















if __name__ == "__main__":
    __main__()