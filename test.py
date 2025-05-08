# test_kmer_analysis.py

from collections import defaultdict, Counter
from pathlib import Path
import pytest
import sys

# If parse_args is in a separate file, update the import below accordingly
from kmer_analysis import parse_args  # <- adjust if module name is different

# Define the functions to test
def read_sequences(input_file):
    sequence = []
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            sequence.append(line.upper())
    return ''.join(sequence)

def get_kmers(sequence, k):
    if len(sequence) < k:
        return defaultdict(Counter)
    kmer_dict = defaultdict(Counter)
    for i in range(len(sequence) - k):
        kmer = sequence[i:i+k]
        next_char = sequence[i+k]
        if set(kmer + next_char).issubset({'A', 'C', 'G', 'T'}):
            kmer_dict[kmer][next_char] += 1
    last_kmer = sequence[-k:]
    if set(last_kmer).issubset({'A', 'C', 'G', 'T'}):
        if last_kmer not in kmer_dict:
            kmer_dict[last_kmer] = Counter()
    return kmer_dict

def write_output(kmer_dict, output_file):
    with open(output_file, 'w') as f:
        header = ['k-mer', 'total_count', 'A', 'C', 'G', 'T']
        f.write('\t'.join(header) + '\n')
        for kmer, counter in sorted(kmer_dict.items()):
            total = sum(counter.values())
            row = [kmer, str(total)] + [str(counter.get(nuc, 0)) for nuc in 'ACGT']
            f.write('\t'.join(row) + '\n')

# ------------------------------
#             TESTS
# ------------------------------

def test_read_sequences(tmp_path):
    test_file = tmp_path / "sample.fa"
    test_file.write_text(">header\nACGT\nTGCA\n")
    result = read_sequences(str(test_file))
    assert result == "ACGTTGCA"

def test_get_kmers_valid():
    sequence = "AATGC"
    k = 2
    result = get_kmers(sequence, k)
    expected = {
        "AA": Counter({"T": 1}),
        "AT": Counter({"G": 1}),
        "TG": Counter({"C": 1}),
    }
    for kmer in expected:
        assert result[kmer] == expected[kmer]
    assert "GC" in result
    assert result["GC"] == Counter()

def test_get_kmers_invalid_characters():
    sequence = "AACNXGT"
    k = 2
    result = get_kmers(sequence, k)
    assert "AA" in result
    assert "GT" in result
    assert "CN" not in result
    assert "NX" not in result
    assert "AC" not in result

def test_write_output(tmp_path):
    kmer_data = {
        "AA": Counter({"T": 2, "G": 1}),
        "AC": Counter({"C": 1}),
        "GT": Counter(),
    }
    output_file = tmp_path / "result.tsv"
    write_output(kmer_data, str(output_file))
    contents = output_file.read_text().strip().splitlines()
    assert contents[0] == "k-mer\ttotal_count\tA\tC\tG\tT"
    assert "AA\t3\t0\t0\t1\t2" in contents
    assert "AC\t1\t0\t1\t0\t0" in contents
    assert "GT\t0\t0\t0\t0\t0" in contents

def test_argument_parsing(monkeypatch):
    test_args = ["kmer_analysis.py", "-i", "input.fa", "-o", "output.tsv", "-k", "4"]
    monkeypatch.setattr(sys, "argv", test_args)
    args = parse_args()
    assert args.input == "input.fa"
    assert args.output == "output.tsv"
    assert args.size == 4

# ------------------------------
#      ADDITIONAL EDGE CASES
# ------------------------------

def test_k_equals_one():
    sequence = "ATCGAT"
    result = get_kmers(sequence, 1)
    assert result["A"]["T"] == 2
    assert result["T"]["C"] == 1
    assert result["C"]["G"] == 1
    assert result["G"]["A"] == 1
    assert "T" in result  # final k-mer

def test_k_equals_sequence_length():
    sequence = "ACGT"
    k = 4
    result = get_kmers(sequence, k)
    assert result == {"ACGT": Counter()}

def test_k_greater_than_sequence_length():
    sequence = "ACGT"
    k = 5
    result = get_kmers(sequence, k)
    assert result == {}

def test_empty_sequence():
    sequence = ""
    k = 3
    result = get_kmers(sequence, k)
    assert result == {}

def test_repeated_character_sequence():
    sequence = "AAAAAA"
    k = 2
    result = get_kmers(sequence, k)
    assert result["AA"]["A"] == 4
    assert sum(result["AA"].values()) == 4

def test_input_with_empty_lines(tmp_path):
    test_file = tmp_path / "seq.fa"
    test_file.write_text("ACGT\n\nTGCA\n\n")
    result = read_sequences(str(test_file))
    assert result == "ACGTTGCA"
