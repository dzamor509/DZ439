import argparse
from collections import defaultdict, Counter

def read_sequences(input_file):
    """
    Reads a DNA sequence from a file (FASTA or plain text).
    Ignores header lines (starting with '>') and blank lines.

    Parameters:
        input_file (str): Path to the input file.

    Returns:
        str: A single concatenated uppercase DNA sequence.
    """
    sequence = []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue  # Skip headers and empty lines
            sequence.append(line.upper())
    return ''.join(sequence)


def get_kmers(sequence, k):
    """
    Extracts all valid k-mers from a DNA sequence and counts the frequency 
    of the nucleotide that immediately follows each k-mer.

    Parameters:
        sequence (str): The full DNA sequence.
        k (int): The length of each k-mer.

    Returns:
        dict: A defaultdict of Counters. Each k-mer maps to a Counter of A/C/G/T that follow it.
    """
    # Return empty dictionary if sequence is too short
    if len(sequence) < k:
        return {}

    kmer_dict = defaultdict(Counter)

    # Slide through the sequence to extract k-mers and their next base
    for i in range(len(sequence) - k):
        kmer = sequence[i:i+k]
        next_char = sequence[i+k]
        # Only count valid nucleotide combinations (A, C, G, T)
        if set(kmer + next_char).issubset({'A', 'C', 'G', 'T'}):
            kmer_dict[kmer][next_char] += 1

    # Include final k-mer if it's valid and has length k
    last_kmer = sequence[-k:]
    if len(last_kmer) == k and set(last_kmer).issubset({'A', 'C', 'G', 'T'}):
        if last_kmer not in kmer_dict:
            kmer_dict[last_kmer] = Counter()

    return kmer_dict


def write_output(kmer_dict, output_file):
    """
    Writes the k-mer statistics to a tab-separated output file.
    Each row includes the k-mer, its total count, and the frequency of A, C, G, T following it.

    Parameters:
        kmer_dict (dict): Dictionary of k-mer to following nucleotide counts.
        output_file (str): Path to output file.
    """
    with open(output_file, 'w') as f:
        header = ['k-mer', 'total_count', 'A', 'C', 'G', 'T']
        f.write('\t'.join(header) + '\n')

        for kmer, counter in sorted(kmer_dict.items()):
            total = sum(counter.values())
            # Write k-mer, total count, and counts of A, C, G, T
            row = [kmer, str(total)] + [str(counter.get(nuc, 0)) for nuc in 'ACGT']
            f.write('\t'.join(row) + '\n')


def parse_args():
    """
    Parses command-line arguments using argparse.

    Returns:
        argparse.Namespace: Contains the input file, output file, and k-mer size.
    """
    parser = argparse.ArgumentParser(
        description="Count k-mers and their next nucleotide from a DNA sequence file."
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input file with DNA sequences (FASTA or plain text)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output TSV file for k-mer statistics')
    parser.add_argument('-k', '--size', type=int, required=True,
                        help='Value of k for k-mer length')
    return parser.parse_args()


def main():
    """
    Main execution function:
    - Parses arguments
    - Reads input sequence
    - Extracts and counts k-mers
    - Writes output to file
    """
    args = parse_args()
    sequence = read_sequences(args.input)
    kmer_dict = get_kmers(sequence, args.size)
    write_output(kmer_dict, args.output)


if __name__ == '__main__':
    main()

