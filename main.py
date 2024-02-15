# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    # hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    # gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    # mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    # br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    # tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    
    # TODO print all of the alignment score between each species BRD2 and human BRD2
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    species_seqs = {
        "Gallus gallus": read_fasta("./data/Gallus_gallus_BRD2.fa")[0],
        "Mus musculus": read_fasta("./data/Mus_musculus_BRD2.fa")[0],
        "Balaeniceps rex": read_fasta("./data/Balaeniceps_rex_BRD2.fa")[0],
        "Tursiops truncatus": read_fasta("./data/tursiops_truncatus_BRD2.fa")[0]
    }

    # Using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    
    alignment_scores = {}

    for species, seq in species_seqs.items():
        score, _, _ = nw.align(hs_seq, seq)
        alignment_scores[species] = score

    # Sort species by similarity (alignment score) in descending order
    sorted_species = sorted(alignment_scores.items(), key=lambda x: x[1], reverse=True)

    # Print species in order of most similar to human BRD
    print("Species in order of most similar to human BRD:")
    for species, score in sorted_species:
        print(f"{species}: {score}")

    # Print all of the alignment score between each species BRD2 and human BRD2
    print("\nAlignment scores between each species BRD2 and human BRD2:")
    for species, score in alignment_scores.items():
        print(f"Human BRD2 vs {species} BRD2: {score}")

    ## Run tests
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
    alignment_score, aligned_seq3, aligned_seq4 = nw.align(seq4, seq3)
    # assert alignment_score == 10
    # assert aligned_seq3 == "MAVHQLIRRP"
    # assert aligned_seq4 == "---MQLIRHP"
    print(alignment_score)
    print(aligned_seq3)
    print(aligned_seq4)

if __name__ == "__main__":
    main()
