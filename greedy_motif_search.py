#!/usr/bin/env python3
def build_profile(motifs, k):
    """
    Build a profile matrix from a list of motifs with Laplace's rule (pseudocounts).
    Returns a dictionary mapping each nucleotide to a list of probabilities for each column.
    """
    t = len(motifs)
    # Initialize counts with pseudocount of 1 for each nucleotide at each position
    counts = {nuc: [1] * k for nuc in "ACGT"}
    for motif in motifs:
        for i, nuc in enumerate(motif):
            counts[nuc][i] += 1
    # Total counts per column is t + 4 because of the 1 added for each nucleotide.
    profile = {}
    for nuc in "ACGT":
        profile[nuc] = [count / (t + 4) for count in counts[nuc]]
    return profile

def profile_probability(kmer, profile):
    """
    Given a k-mer and a profile matrix, return the probability of the k-mer.
    """
    prob = 1.0
    for i, nuc in enumerate(kmer):
        prob *= profile[nuc][i]
    return prob

def profile_most_probable_kmer(text, k, profile):
    """
    Find the profile-most probable k-mer in a string (ties break in favor of the first occurrence).
    """
    max_prob = -1
    most_probable = text[0:k]  # default: first k-mer
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        prob = profile_probability(kmer, profile)
        if prob > max_prob:
            max_prob = prob
            most_probable = kmer
    return most_probable

def score(motifs, k):
    """
    Compute the score of a list of motifs: the sum of mismatches between each column and the consensus.
    """
    score_val = 0
    t = len(motifs)
    for i in range(k):
        count = {"A": 0, "C": 0, "G": 0, "T": 0}
        for motif in motifs:
            count[motif[i]] += 1
        max_count = max(count.values())
        score_val += (t - max_count)
    return score_val

def greedy_motif_search(dna, k, t):
    """
    Greedy motif search using Laplace's rule for pseudocounts.
    """
    best_motifs = [seq[:k] for seq in dna]  # initialize with first k-mer of each sequence
    best_score = score(best_motifs, k)
    
    # iterate through all k-mers in the first sequence
    first_seq = dna[0]
    for i in range(len(first_seq) - k + 1):
        motifs = [first_seq[i:i+k]]
        # build motifs one by one for the remaining sequences
        for j in range(1, t):
            current_profile = build_profile(motifs, k)
            next_motif = profile_most_probable_kmer(dna[j], k, current_profile)
            motifs.append(next_motif)
        current_score = score(motifs, k)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs
    return best_motifs

def main():
    # Read input from file named "input"
    with open("input", "r") as infile:
        lines = [line.strip() for line in infile if line.strip()]
    
    # Parse first line to get k and t
    k, t = map(int, lines[0].split())
    dna = lines[1:1+t]
    
    # Run the greedy motif search
    result = greedy_motif_search(dna, k, t)
    
    # Write result to file named "output"
    with open("output", "w") as outfile:
        for motif in result:
            outfile.write(motif + "\n")

if __name__ == "__main__":
    main()
