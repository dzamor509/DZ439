import argparse
from collections import defaultdict, Counter


def read_sequences(input_file):
    """
    Reads a sequence file (FASTA or plain text) and returns the combined sequence.

    Parameters:
        input_file (str): Path to input file.

    Returns:
        str: The combined DNA sequence in uppercase.
    """
    sequence = []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            sequence.append(line.upper())
    return ''.join(sequence)


def get_kmers(sequence, k):
    """
    Extracts k-mers and the character that follows each k-mer.

    Parameters:
        sequence (str): The DNA sequence.
        k (int): Length of each k-mer.

    Returns:
        dict: Dictionary where keys are k-mers and values are Counters of next characters.
    """
    kmer_dict = defaultdict(Counter)
    for i in range(len(sequence) - k):
        kmer = sequence[i:i+k]
        next_char = sequence[i+k]
        if set(kmer + next_char).issubset({'A', 'C', 'G', 'T'}):
            kmer_dict[kmer][next_char] += 1
    # Count final k-mer without a follower
    last_kmer = sequence[-k:]
    if set(last_kmer).issubset({'A', 'C', 'G', 'T'}):
        if last_kmer not in kmer_dict:
            kmer_dict[last_kmer] = Counter()
    return kmer_dict


def write_output(kmer_dict, output_file):
    """
    Writes k-mer frequencies and their next character frequencies to a file.

    Parameters:
        kmer_dict (dict): Dictionary with k-mers and following character counts.
        output_file (str): Path to output file.
    """
    with open(output_file, 'w') as f:
        header = ['k-mer', 'total_count', 'A', 'C', 'G', 'T']
        f.write('\t'.join(header) + '\n')
        for kmer, counter in sorted(kmer_dict.items()):
            total = sum(counter.values())
            row = [kmer, str(total)] + [str(counter.get(nuc, 0)) for nuc in 'ACGT']
            f.write('\t'.join(row) + '\n')


def parse_args():
    """
    Parses command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments with input, output, and size (k).
    """
    parser = argparse.ArgumentParser(description="Count k-mers and their followers from a DNA sequence file.")
    parser.add_argument('-i', '--input', required=True, help='Input file with DNA sequences (FASTA or plain text)')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file for k-mer statistics')
    parser.add_argument('-k', '--size', type=int, required=True, help='Value of k for k-mer length')
    return parser.parse_args()


def main():
    args = parse_args()
    sequence = read_sequences(args.input)
    kmer_dict = get_kmers(sequence, args.size)
    write_output(kmer_dict, args.output)


if __name__ == '__main__':
    main()
