from collections import Counter
from pathlib import Path
import pytest
import sys

# Adjust import if your main script has a different name
from kmer_analysis import parse_args, read_sequences, get_kmers, write_output

# ------------------------------
#            TESTS
# ------------------------------

#  Test: Basic FASTA-style input parsing, skipping headers and joining lines
def test_read_sequences(tmp_path):
    test_file = tmp_path / "sample.fa"
    test_file.write_text(">header\nACGT\nTGCA\n")
    result = read_sequences(str(test_file))
    assert result == "ACGTTGCA"

#  Test: Valid sequence with k=2, ensures correct k-mers and following nucleotides
def test_get_kmers_valid():
    sequence = "AATGC"
    k = 2
    result = get_kmers(sequence, k)
    expected = {
        "AA": Counter({"T": 1}),
        "AT": Counter({"G": 1}),
        "TG": Counter({"C": 1}),
        "GC": Counter()  # Last k-mer with no following nucleotide
    }
    assert dict(result) == expected

#  Test: Handles invalid characters by skipping affected k-mers
def test_get_kmers_invalid_characters():
    sequence = "AACNXGT"
    k = 2
    result = get_kmers(sequence, k)
    assert "AA" in result
    assert "GT" in result
    assert "CN" not in result
    assert "NX" not in result
    assert "AC" not in result  # AC followed by 'N' is invalid

#  Test: Output formatting correctness for TSV file
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
    assert contents[1] == "AA\t3\t0\t0\t1\t2"
    assert contents[2] == "AC\t1\t0\t1\t0\t0"
    assert contents[3] == "GT\t0\t0\t0\t0\t0"

#  Test: Argument parsing logic via simulated command line input
def test_argument_parsing(monkeypatch):
    test_args = ["kmer_analysis.py", "-i", "input.fa", "-o", "output.tsv", "-k", "4"]
    monkeypatch.setattr(sys, "argv", test_args)
    args = parse_args()
    assert args.input == "input.fa"
    assert args.output == "output.tsv"
    assert args.size == 4

# ------------------------------
#      EDGE CASE TESTING
# ------------------------------

#  Test: k = 1; ensures nucleotide transitions are counted correctly
def test_k_equals_one():
    sequence = "ATCGAT"
    result = get_kmers(sequence, 1)
    assert result["A"]["T"] == 2
    assert result["T"]["C"] == 1
    assert result["C"]["G"] == 1
    assert result["G"]["A"] == 1
    assert "T" in result  # Final nucleotide still stored with empty context

#  Test: k equals length of sequence; final k-mer should exist, no next character
def test_k_equals_sequence_length():
    sequence = "ACGT"
    k = 4
    result = get_kmers(sequence, k)
    assert dict(result) == {"ACGT": Counter()}

#  Test: k > sequence length; function should return empty dict
def test_k_greater_than_sequence_length():
    sequence = "ACGT"
    k = 5
    result = get_kmers(sequence, k)
    assert dict(result) == {}

#  Test: Empty sequence; should return empty k-mer dictionary
def test_empty_sequence():
    sequence = ""
    k = 3
    result = get_kmers(sequence, k)
    assert dict(result) == {}

#  Test: Repeated characters; verifies frequency accumulation for same k-mer
def test_repeated_character_sequence():
    sequence = "AAAAAA"
    k = 2
    result = get_kmers(sequence, k)
    assert result["AA"]["A"] == 4
    assert sum(result["AA"].values()) == 4

#  Test: Sequence input with empty lines; ensures those lines are ignored
def test_input_with_empty_lines(tmp_path):
    test_file = tmp_path / "seq.fa"
    test_file.write_text("ACGT\n\nTGCA\n\n")
    result = read_sequences(str(test_file))
    assert result == "ACGTTGCA"
