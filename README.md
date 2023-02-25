# PyGDS-Final 

[![CI Build](https://github.com/bashgame/PyGDS-Final/actions/workflows/ci-build.yaml/badge.svg)](https://github.com/bashgame/PyGDS-Final/actions/workflows/ci-build.yaml)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python 3.9](https://img.shields.io/badge/Python-3.9-green.svg)](https://shields.io/)

Final Project for Python for Genomic Data Science

## Program requirements description
Write a program that taks as input a file containing DNA sequences in multi-
FASTA format, and computes the answers to the following questions:
* How many records are in the file?
* What are the lengths of the sequences in the file?
    * What is the longest sequence?
    * What is the shortest sequence?
    * Are there multiple longest/shortest sequences?
    * What are their identifiers?
* Which ORFs are present on forward Reading Strands 1, 2, and 3?
    * What is the length of the longest ORF in the FILE?
    * What is the identifier of the sequence containing the longest ORF?
    * For a given sequence identifier, what is the longest ORF contained in the
      sequence represented by that identifier?
    * What is the starting position of the longest ORF in the sequence that
      contains it? (character position in the str, not index of char)
* Given a length *n*, how many Repeats of length *n* occur in all sequences
  contained in the FASTA file?
    * How many times does each Repeat occur?
    * Which repeat of length *n* occurs most frequently?

### Definitions
* multi-FASTA format
    * Records are defined as a single-line header, starting with '>', followed
      by lines of sequence data.
    * The sequence identifier is the word following '>', the rest of the header
      is an optional description of the entry. There should be no whitespace
      between '>' and the identifier.
* A Reading Frame is a way of dividing the DNA sequence of nucleotides into a
  set of consecutive, non-overlapping triplets (or codons).
    * There are 6 possible Reading Frames based on starting position and
      direction.
    * The 3 forward (5' to 3') Reading Frames are:
        * Starting from the 0th character
        * Starting from the 1st character
        * Starting from the 2nd character
    * An ORF, or Open Reading Frame, is the part of a reading frame that has
      the potential to encode a protein.
        * An ORF always starts with ATG, the start codon
        * An ORF is always terminated by one of the stop codons (TAA, TAG, TGA)
* A Repeat is a substring of a sequence that occurs more than once in the
  sequence.
    * We will only be concerned with the forward strand of DNA
    * We will allow repeats to overlap themselves
        * eg ACACA contains 2 repeats of ACA

## Testing framework
Unit testing with unittest and coverage, linting with flake8  
Example file is located in tests/fixtures/dna.example.fasta
### TODO: Set up github workflow for automated testing on PRs

## Containerization
Includes devcontainer and Dockerfile for running in VSCode remote container
