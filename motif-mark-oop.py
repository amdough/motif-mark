import math
import re
import argparse
import pycairo
import bioinfo 


#----------------------------------------------------------------------------
# CLASS: Sequence
# This class represents an input fasta sequence.
# Stores header, sequence, and computes exon and intron sequences based on the feature class.
# Methods:
# 
# ----------------------------------------------------------------------------
class Sequence:
    def __init__(self, header:str, seq:str):
        '''This class represents an input fasta sequence.'''
        self.header = header
        self.seq = seq
        self.length = len(seq)
        self.exons = self._find_regions(uppercase=True)
        self.introns = self._find_regions(uppercase=False)

    def _find_regions(self, uppercase:bool) -> list[tuple[int, int]]:
        '''This method finds the regions of the sequence that are either uppercase (exons) or lowercase (introns) 
        Returns:
        a list of tuples with the start and end positions of each region.'''
        regions = []
        in_region = False
        start = 0
        for i, ch in enumerate(self.seq):
            is_upper = ch.isupper()
            if uppercase:
                in_target = is_upper
            else:
                in_target = not is_upper
            if in_target and not in_region:
                start = i
                in_region = True
            elif not in_target and in_region:
                regions.append((start, i))
                in_region = False
        if in_region:
            regions.append((start, len(self.seq)))
        return regions
    
    def __repr__(self):
        return f"Sequence(header={self.header!r}, length={self.length})"

#----------------------------------------------------------------------------
# CLASS: Motif
# This class represents an input motif.
#----------------------------------------------------------------------------

class Motif: 
    def __init__(self, name:str, start: int, end:int, gene:Sequence):
        '''This class represents an input motif.'''
        self.name = name
        self.start = start
        self.end = end
        self.gene = gene
    
    def __repr__(self):
        return f"Motif(name={self.name!r}, start={self.start}, end={self.end})"


#----------------------------------------------------------------------------
# Function: parse_fasta
# This function parses a FASTA file and returns a list of Seq objects. Multi lines are joined into one string.
#----------------------------------------------------------------------------

def parse_fasta(fasta_file: str) -> list[Sequence]:
    """Parses a FASTA file and returns a list of Seq objects. Multi lines are joined into one string."""
    genes = []
    header = None
    seq_lines = []
    with open(fasta_file, 'r') as fh:
        for line in fh:
            line = line.strip('\n')
            if line.startswith('>'):
                if header is not None:
                    genes.append(Sequence(header, ''.join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        genes.append(Sequence(header, ''.join(seq_lines)))
    return genes

def parse_motifs(motif_file: str) -> list[str]:
    """Return a list of motifs from the input motif file; one motif per line."""
    motifs = []
    with open(motif_file, 'r') as fh:
        for line in fh:
            motif = line.strip()
            if motif:
                motifs.append(motif)
    return motifs

def find_motif_occurrences(genes: list[Sequence], motif_names: list[str]) -> dict[str, list[Motif]]:
    """Finds the occurrences of each motif in the gene sequences and returns a dictionary mapping motif names to lists of Motif objects."""
    occurrences = {}
    for gene in genes:
        motif_list = []
        seq_upper = gene.seq.upper()
        for motif_name in motif_names:
            pattern = bioinfo.convert_motif(motif_name)






    return occurrences
    







#----------------------------------------------------------------------------
## Main function for motif-mark oop implementation.
#----------------------------------------------------------------------------
def get_args():
        parser = argparse.ArgumentParser(description="Visualizes the distribution of a motif in a fasta sequence. Outputs an image file with the visualization of the number of times the motif is found in the exon and intron sequences.")
        parser.add_argument("-f", help="Designates absolute file path to input fasta file", type=str, required =True)
        parser.add_argument("-m", help="Designates absolute file path to input motifs file", type=str, required=True)
        parser.add_argument("-o", help="Designates output file for motif-mark stats", type=str, required=True)

        return parser.parse_args()

def __main__():

    args = get_args()
    input_fasta = args.f 
    input_motifs = args.m 
    output = args.o 









if __name__ == "__main__":
    __main__()