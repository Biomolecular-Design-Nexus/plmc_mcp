#!/usr/bin/env python3
"""
Remove gaps from query sequence in a2m alignment file.

This script removes all positions where the query sequence (first sequence)
contains gaps, effectively creating an alignment with a gap-free query.
"""

import sys
import argparse


def read_a2m(filename):
    """
    Read a2m file and return list of (header, sequence) tuples.
    """
    sequences = []
    current_header = None
    current_seq = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                # Start new sequence
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))

    return sequences


def remove_query_gaps(sequences):
    """
    Remove positions where the query (first) sequence has gaps.

    Args:
        sequences: List of (header, sequence) tuples

    Returns:
        List of (header, sequence) tuples with gaps removed
    """
    if not sequences:
        return sequences

    query_header, query_seq = sequences[0]

    # Find positions where query does NOT have gaps (both '.' and '-' are gaps in a2m)
    non_gap_positions = [i for i, char in enumerate(query_seq)
                         if char not in '.-']

    # Create new sequences keeping only non-gap positions
    cleaned_sequences = []
    for header, seq in sequences:
        cleaned_seq = ''.join(seq[i] for i in non_gap_positions)
        cleaned_sequences.append((header, cleaned_seq))

    return cleaned_sequences


def write_a2m(sequences, filename):
    """
    Write sequences to a2m file.
    """
    with open(filename, 'w') as f:
        for header, seq in sequences:
            f.write(f"{header}\n")
            # Write sequence in chunks of 80 characters for readability
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Remove gaps from query sequence in a2m alignment file'
    )
    parser.add_argument('input', help='Input a2m file')
    parser.add_argument('output', help='Output a2m file with query gaps removed')

    args = parser.parse_args()

    # Read input file
    print(f"Reading {args.input}...")
    sequences = read_a2m(args.input)
    print(f"Found {len(sequences)} sequences")

    if not sequences:
        print("Error: No sequences found in input file")
        sys.exit(1)

    # Get original lengths
    query_header, query_seq = sequences[0]
    original_length = len(query_seq)
    original_gaps = sum(1 for c in query_seq if c in '.-')

    print(f"Query sequence: {query_header}")
    print(f"Original alignment length: {original_length}")
    print(f"Query gaps to remove: {original_gaps}")

    # Remove query gaps
    cleaned_sequences = remove_query_gaps(sequences)

    new_length = len(cleaned_sequences[0][1])
    print(f"New alignment length: {new_length}")

    # Write output file
    print(f"Writing {args.output}...")
    write_a2m(cleaned_sequences, args.output)

    print("Done!")


if __name__ == '__main__':
    main()
