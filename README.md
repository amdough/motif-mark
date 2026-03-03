# Motif-mark
## Author: Amanda Dougherty
## Python script that visualizes motifs on sequences using object-oriented programming and pycairo.
Given a FASTA file and a file list of motifs, ```motif-mark-oop.py``` will produce an output .png that depicts introns [horizontal lines], exons [rectangles], and motifs [colored rectangles at position on sequence] drawn to scale relative to the sequence length. 
### System Requirements 
Requires installation of python, pycairo
    - conda create -n my_pycairo pycairo
    - conda activate my_pycairo
### Usage

```python motif-mark-oop.py -f /Figure_1.fasta -m Fig_1_motifs.txt```

## Files: 

- [Main script](motif-mark-oop.py)
    - parses input, finds motif occurrences, generates the figure, and outputs script stats.

- [Bioinfo.py](bioinfo.py) 
    - Script needs a function `convert_motif(motif:str) -> str' to convert motif strings including IUPAC ambiguity codes (also defined within bioinfo.py)