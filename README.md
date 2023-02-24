# PyGDS-Final 
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
    * How many times does each Reapeat occur?
    * Which repeat of length *n* occurs most frequently?


## Testing framework
Unit testing with unittest and coverage, linting with pylint
### TODO: Set up github workflow for automated testing on PRs

## Containerization
Includes devcontainer and Dockerfile for running in VSCode remote container
