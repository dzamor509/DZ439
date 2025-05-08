
# Kmer Analysis

## Introduction

This Python script processes a file containing DNA sequences and counts the frequency of all k-mers (substrings of length `k`) and the frequency of each nucleotide that follows them. It then outputs the results to a formatted tab-separated values (TSV) file.

## Features

* Reads DNA sequences from a FASTA or plain text file
* Allows user-specified `k` for k-mer length
* Outputs frequency of each k-mer and its immediate following nucleotides (A, C, G, T)
* Designed to skip invalid characters (e.g., N, X) in sequence data
* Fully tested with `pytest`, including edge cases

## Usage

### Command-line
python kmer_analysis.py -i <input_file> -o <output_file> -k <k_value>


## Output Format

The output TSV file includes:

* Each k-mer
* Total count of the k-mer
* Number of times each nucleotide (A, C, G, T) follows the k-mer


## Testing

This project includes a comprehensive `pytest` test suite in `tests.py`. To run tests:


### Covered Test Cases

* Argument parsing with `argparse`
* Valid input parsing from FASTA-style files
* Detection and handling of invalid characters
* Edge cases including:
  * `k = 1`
  * `k = length of sequence`
  * `k > length of sequence`
  * Empty input
  * Repeated characters
  * Input files with empty lines

## File Structure

kmer_analysis.py # Main script
reads.fa # Sample input file
test.py # Pytest-based unit tests
output.tsv # Example output (generated)


## Requirements

* Python 3.6+
* `pytest` library (install via `pip install pytest`)

## Author

Danielo Zamor


