# OoCA_Motif-mark
## Amanda Dougherty
## OOCA for python script using object-oriented code to visualize motifs on sequences (1 image per .fa file; .png file please)

### 2.9.2026
1. Created pycairo conda env
    - conda create -n my_pycairo
    - conda activate my_pycairo
    - conda install python=3.14.3
    - conda install -c conda-forge pycairo

2. Planned Classes & Interactions
- Sequence
    - = 1 fasta entry
    - string
    - **Potential Methods**
        - length()
        - add_feature(feat)

- Features
    - = a region or feature of interest within the sequence record
    - will have **subclasses** that designate type and style (case), and store location, label, and how it would appear on the canvas:
        - exon
        - intron
        - motif_hit

- Motif
    - = specific pattern of interest
    - will name and specify style for canvas
    - **Potential Methods**
        - find_motifs(sequence)
            - returns a list of all motifs

Other functions or potential classes: 
- a way to define and handle ambiguous nucleotides that can be referenced (dict)
- a fasta sequence parser - parse_line() - that can read in the sequence and recognize style/case (exon/intron) and store accordingly
- potentially some pycairo specific classes/functions to generalize drawing conventions (e.g. draw_sequence(); draw_intron/exon(); etc)