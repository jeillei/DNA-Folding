from DNA import DNA, Strand
from utils import update_pos
from animation import animate_dna  

if __name__ == "__main__":
    strand = Strand('ACTG' * 30)
    dna = DNA(strand)
    animate_dna(dna, steps=500, save_path="outputs/dna_folding.mp4")
