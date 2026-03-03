#!/usr/bin/env python3

"""
motif-mark-oop.py

Visualize motif occurrences within gene sequences using Pycairo.

1. Reads a FASTA file in which exons are UPPERCASE, and introns are lowercase letters
2. Reads a motif file and stores possible motifs - can include IUPAC ambiguity codes, converted to regex using bioinfo.convert_motif() (separate script)
3. Finds all motifs in each gene sequence from the FASTA file
4. Generates a .png image that visualizes all gene sequences, their introns/exons, and motifs (including overlaps)
5. Writes out a stats summary of counts in .tsv format

Usage: 

python motif-mark-oop.py -f /Figure_1.fasta -m Fig_1_motifs.txt

Outputs: 
Figure_1.png
Figure_1_motif_stats.tsv

Requires:
-pycairo
-bioinfo module with convert_motif() 

"""

import math
import os
import re
import argparse
import cairo
import bioinfo # script requires the bioinfo module, which contains the convert_motif function used to convert motif names to regular expression patterns.

#----------------------------------------------------------------------------
# Constants for pycairo drawing
#----------------------------------------------------------------------------

left_margin = 180
right_margin = 20
top_margin = 20
bottom_margin = 120

row_height = 140
baseline_offset = 70

exon_height = 18
motif_height = 18
motif_y_gap = 10

font_main = 14
font_small = 11

# Palette for up to ~5 motifs
PALETTE = [
    (0.20, 0.45, 0.80),
    (0.85, 0.35, 0.20),
    (0.20, 0.70, 0.35),
    (0.60, 0.30, 0.70),
    (0.70, 0.65, 0.20),
]
TARGET_SEQ_WIDTH = 2200  # pixels for the longest sequence region

COLOR_TEXT = (0.05, 0.05, 0.05)
COLOR_INTRON = (0.55, 0.55, 0.55)
COLOR_EXON = (0.15, 0.15, 0.15)
MOTIF_ALPHA = 0.55



#----------------------------------------------------------------------------
# Data Classes
#----------------------------------------------------------------------------

#----------------------------------------------------
# CLASS: Exon
# This class represents an exon region of a sequence.
#----------------------------------------------------

""" This class represents an exon region of a sequence. It stores the name of the exon, its start and end positions, and a reference to the parent Sequence object. """

class Exon: 
    def __init__(self, start:int, end: int, gene_number: int):
        self.start = start
        self.end = end
        self.gene_number = gene_number
    
    def draw(self, ctx:cairo.Context, x0: float, px_per_base: float, y_baseline: float):
       ctx.set_source_rgb(*COLOR_EXON)
       x = x0 + self.start * px_per_base
       w = max(1.0, (self.end - self.start) * px_per_base)
       y = y_baseline - exon_height / 2
       ctx.rectangle(x, y, w, exon_height)
       ctx.fill()



#----------------------------------------------------------------------------
# CLASS: Motif
# This class represents a single motif occurrence on a gene.
#----------------------------------------------------------------------------

class Motif:
    """
    This class represents a single motif occurrence on a gene. It stores start and end position, as well as a reference to the parent sequence. 
    """
    def __init__(self, name:str, start:int, end: int, gene_number: int, color: tuple[float, float, float]):
        self.name = name
        self.start = start
        self.end = end
        self.gene_number = gene_number
        self.color = color 

   
    def draw(self, ctx:cairo.Context, x0: float, px_per_base: float, y_baseline: float): # Draws the hit as a transparent rectangle with an outline for easier visibility onto the Cairo context
        r, g, b = self.color
        ctx.set_source_rgba(r, g, b, MOTIF_ALPHA)
        x = x0 + self.start * px_per_base
        w = max(1.0, (self.end - self.start) * px_per_base)
        y = (y_baseline - exon_height / 2) - motif_y_gap - motif_height
        ctx.rectangle(x, y, w, motif_height)
        ctx.fill()
        # Outline for visibility
        ctx.set_source_rgba(0, 0, 0, 0.35)
        ctx.set_line_width(1)
        ctx.rectangle(x, y, w, motif_height)
        ctx.stroke()

                   

#----------------------------------------------------------------------------
# CLASS: Gene
# This class represents an input fasta record.
# ----------------------------------------------------------------------------
class Gene:
    def __init__(self, header:str, seq:str, gene_number: int):
        '''This class represents an input fasta record. Stores header, sequence, and computes exon and intron sequences based on the feature class.'''
        self.header = header.strip()
        self.seq = seq.strip()
        self.gene_number = gene_number
        self.length = len(self.seq)

        self.exons: list[Exon] = self._find_exons()
        self.motifs: list[Motif] = []

    def _find_exons(self) -> list[Exon]:
        '''This method finds the exon regions of the sequence, which are represented by uppercase letters. 
        Returns:
        a list of Exon objects with the start and end positions of each exon.'''
        exons = []
        in_exon = False
        start = 0

        for i, ch in enumerate(self.seq):
            if ch.isupper() and not in_exon:
                start = i
                in_exon = True
            elif (not ch.isupper()) and in_exon:
                exons.append(Exon(start, i, self.gene_number))
                in_exon = False
        if in_exon:
            exons.append(Exon(start, len(self.seq), self.gene_number))
        return exons
        
    def title(self) -> str:
        """ Get the title from the header"""
        first = self.header.split()[0] if self.header else f"gene_{self.gene_number}"
        return f"{first} (len={self.length})"

    def y_baseline(self) -> float:
        """Calculate the y-coordinate for the baseline of this gene's row in the visualization."""
        row_top = top_margin + self.gene_number * row_height
        return row_top + baseline_offset
    
    def x_start(self) -> float:
        """Calculate the x-coordinate for the start of the sequence in the visualization."""
        return left_margin
    
    def draw(self, ctx: cairo.Context, px_per_base: float):
        """Draw the gene's exons and motifs on the given cairo context."""
        x0 = self.x_start()
        y_baseline = self.y_baseline()
        x1 = x0 + self.length * px_per_base

        # Labels
        ctx.set_source_rgb(*COLOR_TEXT)
        ctx.set_font_size(font_small)
        ctx.move_to(12, y_baseline + 4)
        ctx.show_text(self.title())

        # Intron Baseline (full gene length)
        ctx.set_source_rgb(*COLOR_INTRON)
        ctx.set_line_width(2)
        ctx.move_to(x0, y_baseline)
        ctx.line_to(x1, y_baseline)
        ctx.stroke()

        # Exons
        for exon in self.exons:
            exon.draw(ctx, x0, px_per_base, y_baseline)

        # Motifs - drawn after exons so they appear on top
        self.motifs.sort(key = lambda m: (m.start, -(m.end-m.start)))  # Sort by start position, then by length (longer first)
        for motif in self.motifs:
            motif.draw(ctx, x0, px_per_base, y_baseline)

        # Labels 
        ctx.set_source_rgb(*COLOR_TEXT)
        ctx.set_font_size(10)
        ctx.move_to(x0, y_baseline + 22)
        ctx.show_text("0")
        ctx.move_to(x1 - 30, y_baseline + 22)
        ctx.show_text(str(self.length))

#----------------------------------------------------------------------------
# Helper functions
#----------------------------------------------------------------------------

def parse_fasta(fasta_file: str) -> list[Gene]:
    """Parses a FASTA file and returns a list of Seq objects. Multi lines are joined into one string."""
    genes = []
    header = None
    seq_lines = []

    with open(fasta_file, 'r') as fh:
        for line in fh:
            line = line.strip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    genes.append(Gene(header, ''.join(seq_lines), gene_number=len(genes)))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())

    if header is not None:
        genes.append(Gene(header, ''.join(seq_lines), gene_number=len(genes)))
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

# def assign_motif_lanes(motifs: list[Motif]) -> int:
#     """Assigns lane numbers to motifs to avoid overlaps in the visualization. 
#     Modifies the Motif objects in place and returns the total number of lanes used."""
#     motifs.sort(key = lambda h: (h.start, -(h.end - h.start)))  # Sort by start position, then by length (longer first)
#     lanes: list[list[Motif]] = []

#     for motif in motifs:
#         placed = False
#         for lane_num, lane in enumerate(lanes):
#             if all(not motif.overlaps(other) for other in lane):
#                 motif.lane = lane_num
#                 lane.append(motif)
#                 placed = True
#                 break
#         if not placed:
#             motif.lane = len(lanes)
#             lanes.append([motif])
#     return len(lanes)


def add_motifs(genes: list[Gene], motif_names: list[str]):
   color_map = {m: PALETTE[i % len(PALETTE)] for i, m in enumerate(motif_names)}
   for gene in genes:
       seq_upper = gene.seq.upper()
       hits = []

       for motif in motif_names:
           pattern = bioinfo.convert_motif(motif)
           # detect overlaps
           regex = re.compile(rf"(?=({pattern}))", re.IGNORECASE)

           for match in regex.finditer(seq_upper):
               start = match.start()
               end = start + len(motif)
               hits.append(Motif(motif, start, end, gene.gene_number, color_map[motif]))
        
      
       gene.motifs = hits 
   return color_map

def is_in_exon(gene: Gene, start: int, end: int) -> bool:
    """
    Return True if the motif hit overlaps ANY exon region.
    Exons are uppercase runs stored as Exon objects in gene.exons.
    """
    for exon in gene.exons:
        # overlap test for half-open intervals [start,end) and [exon.start, exon.end)
        if not (end <= exon.start or exon.end <= start):
            return True
    return False


def write_motif_stats(genes: list[Gene], out_tsv: str) -> None:
    """
    Write a TSV of all motif hits + a summary section at the top.
    """
    # counts across all genes
    total_counts: dict[str, int] = {}
    exon_counts: dict[str, int] = {}
    intron_counts: dict[str, int] = {}

    # collect rows
    rows = []
    for gene in genes:
        for hit in gene.motifs:
            motif = hit.name
            total_counts[motif] = total_counts.get(motif, 0) + 1

            in_exon = is_in_exon(gene, hit.start, hit.end)
            if in_exon:
                exon_counts[motif] = exon_counts.get(motif, 0) + 1
            else:
                intron_counts[motif] = intron_counts.get(motif, 0) + 1

            # both coordinate styles (0-based and 1-based)
            rows.append((
                gene.header,
                gene.length,
                motif,
                hit.start,          # 0-based start
                hit.end,            # 0-based end (exclusive)
                hit.start + 1,      # 1-based start
                hit.end,            # 1-based end (inclusive if you interpret end-exclusive carefully)
                "exon" if in_exon else "intron"
            ))

    # write file
    with open(out_tsv, "w") as out:
        # summary
        out.write("# motif summary (all genes)\n")
        out.write("# motif\ttotal_hits\texon_hits\tintron_hits\n")
        motifs_sorted = sorted(total_counts.keys())
        for m in motifs_sorted:
            out.write(
                f"{m}\t{total_counts.get(m,0)}\t{exon_counts.get(m,0)}\t{intron_counts.get(m,0)}\n"
            )

        out.write("\n# hit list\n")
        out.write("gene_header\tgene_length\tmotif\tstart0\tend0_excl\tstart1\tend1\toverlaps_exon\n")
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")
            
            
#----------------------------------------------------------------------------
# Generating the Context drawing
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# CLASS: MotifFigure
# This class represents the overall figure for visualizing motifs in genes.
#----------------------------------------------------------------------------

class MotifFigure:
    def __init__(self, genes: list[Gene], motif_names: list[str], motif_colors: dict[str, tuple[float, float, float]], out_png: str):
        self.genes = genes
        self.motif_names = motif_names
        self.motif_colors = motif_colors
        self.out_png = out_png

    def render(self):
        max_len = max(g.length for g in self.genes)
        px_per_base = max(1.0, TARGET_SEQ_WIDTH / max_len)

        seq_width = math.ceil(max_len * px_per_base)
        width = int(left_margin + seq_width + right_margin)
        height = int(top_margin + len(self.genes) * row_height + bottom_margin)

        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        ctx = cairo.Context(surface)

        # background
        ctx.set_source_rgb(1, 1, 1)
        ctx.paint()

        # title
        ctx.set_source_rgb(*COLOR_TEXT)
        ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(font_main)
        ctx.move_to(left_margin, 28)
        ctx.show_text(f"Motif Mark: {os.path.basename(self.out_png)}")

        # draw genes
        for gene in self.genes:
            gene.draw(ctx, px_per_base)

        # legend
        self._draw_legend(ctx, height - 70)

        surface.write_to_png(self.out_png)

    def _draw_legend(self, ctx: cairo.Context, y: float):
        ctx.set_source_rgb(*COLOR_TEXT)
        ctx.set_font_size(font_small)
        ctx.move_to(12, y)
        ctx.show_text("Legend:")

        # exon sample
        ctx.set_source_rgb(*COLOR_EXON)
        ctx.rectangle(80, y - 12, 30, 14)
        ctx.fill()
        ctx.set_source_rgb(*COLOR_TEXT)
        ctx.move_to(120, y)
        ctx.show_text("Exon (UPPERCASE)")

        # intron sample
        ctx.set_source_rgb(*COLOR_INTRON)
        ctx.set_line_width(2)
        ctx.move_to(80, y + 18)
        ctx.line_to(110, y + 18)
        ctx.stroke()
        ctx.set_source_rgb(*COLOR_TEXT)
        ctx.move_to(120, y + 22)
        ctx.show_text("Intron (lowercase baseline)")

        # motif samples
        ctx.set_source_rgb(*COLOR_TEXT)
        ctx.move_to(12,y + 48)
        ctx.show_text("Motifs:")

        x = 80
        y2 = y + 48
        for m in self.motif_names:
            r,g,b = self.motif_colors[m]
            ctx.set_source_rgba(r,g,b, MOTIF_ALPHA) 
            ctx.rectangle(x, y2 - 10, 18, 10)
            ctx.fill()

            ctx.set_source_rgb(*COLOR_TEXT)
            ctx.move_to(x + 24, y2)
            ctx.show_text(m)
            x += 24 + 12 * len(m) + 40


#----------------------------------------------------------------------------
## Main 
#----------------------------------------------------------------------------
def get_args():
        parser = argparse.ArgumentParser(description="Visualizes the distribution of a motif in a fasta sequence. Outputs an image file with the visualization of the number of times the motif is found in the exon and intron sequences.")
        parser.add_argument("-f", help="Designates absolute file path to input fasta file", type=str, required =True)
        parser.add_argument("-m", help="Designates absolute file path to input motifs file", type=str, required=True)
        
        return parser.parse_args()

def __main__():

    args = get_args()

    genes = parse_fasta(args.f)
    motifs = parse_motifs(args.m)

    if not genes:
        raise SystemExit(f"Error: No valid gene sequences found in {args.f}")
    if not motifs:
        raise SystemExit(f"Error: No valid motifs found in {args.m}")
    
    motif_colors = add_motifs(genes, motifs)

    prefix, _ = os.path.splitext(args.f)
    out_png = prefix + ".png"

    outdir = "/Users/amandadougherty/bioinfo/Bi625/motif-mark/motif-mark"
  
    figure = MotifFigure(genes, motifs, motif_colors, out_png)
    figure.render()
    print(f"Wrote: {out_png}")

    # write stats TSV next to the PNG
    out_tsv = prefix + "_motif_stats.tsv"
    write_motif_stats(genes, out_tsv)
    print(f"Wrote: {out_tsv}")




if __name__ == "__main__":
    __main__()