"""
Pseudolikelihood Maximization for Coevolution Analysis (plmc) for EV+Onehot modeling.

This MCP Server provides tools to generate PLMC model parameters and evolutionary couplings
needed for EV+Onehot fitness prediction modeling.

Tools:
1. plmc_generate_model: Generate PLMC model parameters and evolutionary couplings from A2M alignment
2. plmc_convert_a3m_to_a2m: Convert A3M alignment to A2M format and clean query gaps (required preprocessing)

Reference: https://github.com/debbiemarkslab/plmc
"""

# Standard imports
from typing import Annotated
from pathlib import Path
import os
from fastmcp import FastMCP
from datetime import datetime
import subprocess
import shutil

# Project structure
PROJECT_ROOT = Path(__file__).parent.parent.parent.resolve()
DEFAULT_INPUT_DIR = PROJECT_ROOT / "tmp" / "inputs"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "tmp" / "outputs"

INPUT_DIR = Path(os.environ.get("README_INPUT_DIR", DEFAULT_INPUT_DIR))
OUTPUT_DIR = Path(os.environ.get("README_OUTPUT_DIR", DEFAULT_OUTPUT_DIR))

# Ensure directories exist
INPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Timestamp for unique outputs
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# MCP server instance
readme_mcp = FastMCP(name="readme")

# PLMC directory - load from environment variable or use default
PLMC_DIR = os.path.expanduser(
    os.environ.get(
        "PLMC_DIR",
        str(PROJECT_ROOT / "repo" / "plmc")
    )
)

# Find plmc binary
# Try to find plmc in the environment or system PATH
PLMC_BIN = None
for possible_path in [
    os.path.join(PLMC_DIR, "bin", "plmc"),
    str(PROJECT_ROOT / "env" / "bin" / "plmc"),
    shutil.which("plmc")
]:
    if possible_path and os.path.exists(possible_path):
        PLMC_BIN = possible_path
        break

if PLMC_BIN is None:
    PLMC_BIN = "plmc"  # Fallback to PATH lookup

# Find reformat.pl script for a3m to a2m conversion
REFORMAT_PL = None
for possible_path in [
    str(PROJECT_ROOT / "env" / "bin" / "reformat.pl"),
    str(PROJECT_ROOT / "env" / "scripts" / "reformat.pl"),
    shutil.which("reformat.pl")
]:
    if possible_path and os.path.exists(possible_path):
        REFORMAT_PL = possible_path
        break

if REFORMAT_PL is None:
    REFORMAT_PL = "reformat.pl"  # Fallback to PATH lookup


def read_a2m(filename: str) -> list[tuple[str, str]]:
    """
    Read a2m file and return list of (header, sequence) tuples.

    Args:
        filename: Path to a2m file

    Returns:
        List of (header, sequence) tuples
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


def remove_query_gaps(sequences: list[tuple[str, str]]) -> list[tuple[str, str]]:
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


def write_a2m(sequences: list[tuple[str, str]], filename: str) -> None:
    """
    Write sequences to a2m file.

    Args:
        sequences: List of (header, sequence) tuples
        filename: Output file path
    """
    with open(filename, 'w') as f:
        for header, seq in sequences:
            f.write(f"{header}\n")
            # Write sequence in chunks of 80 characters for readability
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")


def convert_a3m_to_a2m(a3m_file: str, a2m_file: str) -> None:
    """
    Convert alignment from A3M format to A2M format using reformat.pl and clean query gaps.

    This function performs two steps:
    1. Converts A3M to A2M format using reformat.pl
    2. Removes positions where the query sequence has gaps (both '.' and '-')

    Cleaning query gaps is essential for downstream modeling as gaps in the query
    sequence can cause issues in fitness prediction workflows.

    Args:
        a3m_file: Path to input A3M file
        a2m_file: Path to output A2M file
    """
    print(f"Converting A3M to A2M: {a3m_file} -> {a2m_file}")
    print(f"  Using reformat.pl: {REFORMAT_PL}")

    cmd = [REFORMAT_PL, "a3m", "a2m", a3m_file, a2m_file]

    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr)

    print(f"Conversion complete: {a2m_file}")

    # Clean query gaps from the converted A2M file
    print(f"\nCleaning query gaps from A2M file...")
    sequences = read_a2m(a2m_file)
    print(f"  Found {len(sequences)} sequences")

    if sequences:
        query_header, query_seq = sequences[0]
        original_length = len(query_seq)
        original_gaps = sum(1 for c in query_seq if c in '.-')

        print(f"  Query sequence: {query_header}")
        print(f"  Original alignment length: {original_length}")
        print(f"  Query gaps to remove: {original_gaps}")

        # Remove query gaps
        cleaned_sequences = remove_query_gaps(sequences)
        new_length = len(cleaned_sequences[0][1])
        print(f"  New alignment length: {new_length}")

        # Write cleaned sequences back to file
        write_a2m(cleaned_sequences, a2m_file)
        print(f"\nQuery gap cleaning complete: {a2m_file}")
    else:
        print("  Warning: No sequences found in A2M file")


@readme_mcp.tool
def plmc_generate_model(
    alignment_path: Annotated[str | None, "Path to protein sequence alignment file in A2M format"] = None,
    focus_seq_id: Annotated[str, "Focus sequence identifier (protein name) for uppercase column modeling"] = None,
    output_dir: Annotated[str | None, "Output directory for model files (default: tmp/outputs)"] = None,
    out_prefix: Annotated[str | None, "Output file prefix (default: uniref100)"] = None,
    lambda_e: Annotated[float, "L2 regularization coefficient for couplings (default: 16.2 for ev+onehot)"] = 16.2,
    lambda_h: Annotated[float, "L2 regularization coefficient for single sites (default: 0.01 for ev+onehot)"] = 0.01,
    max_iterations: Annotated[int, "Maximum number of optimization iterations (default: 200 for ev+onehot)"] = 200,
    theta: Annotated[float, "Sequence reweighting threshold (default: 0.2 for ev+onehot)"] = 0.2,
) -> dict:
    """
    Generate PLMC model parameters and evolutionary couplings for EV+Onehot fitness prediction.

    This tool runs plmc with preset parameters optimized for EV+Onehot modeling and generates
    the required output files:
    - model_params: Binary parameter file containing model parameters
    - EC file: Text file containing evolutionary coupling scores

    Input is protein alignment file in A2M format (use plmc_convert_a3m_to_a2m if you have A3M).

    Preset parameters are based on ev_onehot/scripts/plmc.sh:
    - lambda_e: 16.2 (L2 regularization for couplings)
    - lambda_h: 0.01 (L2 regularization for single sites)
    - max_iterations: 200 (optimization iterations)
    - theta: 0.2 (sequence reweighting threshold)
    """
    # Input validation
    if alignment_path is None:
        raise ValueError("Path to protein alignment file (A2M format) must be provided")

    if focus_seq_id is None:
        raise ValueError("Focus sequence identifier must be provided")

    alignment_file = Path(alignment_path)
    if not alignment_file.exists():
        raise FileNotFoundError(f"Alignment file not found: {alignment_path}")

    # Set output directory
    if output_dir is None:
        out_dir = OUTPUT_DIR
    else:
        out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Set output prefix
    if out_prefix is None:
        out_prefix = "uniref100"

    # Output files (following ev_onehot naming convention)
    params_file = out_dir / f"{out_prefix}.model_params"
    couplings_file = out_dir / f"{out_prefix}.EC"

    # Run plmc with ev+onehot preset parameters
    # Command format from ev_onehot/scripts/plmc.sh:
    # plmc -o model_params -c EC -f PROTEIN -le 16.2 -lh 0.01 -m 200 -t 0.2 -g alignment.a2m
    cmd = [
        str(PLMC_BIN),
        "-o", str(params_file),
        "-c", str(couplings_file),
        "-f", focus_seq_id,
        "-le", str(lambda_e),
        "-lh", str(lambda_h),
        "-m", str(max_iterations),
        "-t", str(theta),
        "-g",  # Gap mode
        str(alignment_file)
    ]

    print(f"Running plmc to generate EV+Onehot model inputs...")
    print(f"  Using plmc binary: {PLMC_BIN}")
    print(f"  Focus sequence: {focus_seq_id}")
    print(f"  Alignment file: {alignment_file}")
    print(f"  Command: {' '.join(cmd)}")
    print()

    result = subprocess.run(cmd, check=True, capture_output=True, text=True)

    # Print output in real-time
    if result.stdout:
        print("PLMC Output:")
        print(result.stdout)
    if result.stderr:
        print("PLMC Stderr:")
        print(result.stderr)

    print(f"\nPLMC model generation complete!")

    return {
        "message": f"PLMC model generated successfully for {focus_seq_id}",
        "reference": "https://github.com/debbiemarkslab/plmc",
        "parameters": {
            "lambda_e": lambda_e,
            "lambda_h": lambda_h,
            "max_iterations": max_iterations,
            "theta": theta,
            "focus_seq": focus_seq_id
        },
        "artifacts": [
            {
                "description": "Model parameters file (for EV predictor)",
                "path": str(params_file.resolve())
            },
            {
                "description": "Evolutionary couplings file (for EV predictor)",
                "path": str(couplings_file.resolve())
            }
        ]
    }


@readme_mcp.tool
def plmc_convert_a3m_to_a2m(
    a3m_file_path: Annotated[str | None, "Path to input A3M alignment file"] = None,
    a2m_file_path: Annotated[str | None, "Path to output A2M alignment file (optional, will auto-generate if not provided)"] = None,
    out_prefix: Annotated[str | None, "Output file prefix (used if a2m_file_path not provided)"] = None,
) -> dict:
    """
    Convert protein sequence alignment from A3M format to A2M format and clean query gaps.

    This is a required preprocessing step before running plmc_generate_model, as plmc
    requires A2M format alignments for EV+Onehot modeling.

    This tool performs two operations:
    1. Converts A3M to A2M format using reformat.pl
    2. Removes all positions where the query (first) sequence contains gaps

    Query gap cleaning is essential because gaps in the query sequence can cause issues
    in downstream fitness prediction modeling workflows.

    Input is A3M alignment file and output is A2M alignment file with cleaned query sequence.
    """
    # Input validation
    if a3m_file_path is None:
        raise ValueError("Path to A3M alignment file must be provided")

    a3m_file = Path(a3m_file_path)
    if not a3m_file.exists():
        raise FileNotFoundError(f"A3M file not found: {a3m_file_path}")

    # Determine output file path
    if a2m_file_path is None:
        if out_prefix is None:
            # Use input filename without extension as prefix
            out_prefix = a3m_file.stem
        a2m_file = OUTPUT_DIR / f"{out_prefix}.a2m"
    else:
        a2m_file = Path(a2m_file_path)
        # Create parent directory if it doesn't exist
        a2m_file.parent.mkdir(parents=True, exist_ok=True)

    # Convert using the helper function
    convert_a3m_to_a2m(str(a3m_file), str(a2m_file))

    return {
        "message": f"Successfully converted A3M to A2M format and cleaned query gaps",
        "reference": "https://github.com/soedinglab/hh-suite",
        "artifacts": [
            {
                "description": "Input A3M alignment file",
                "path": str(a3m_file.resolve())
            },
            {
                "description": "Output A2M alignment file with query gaps removed (ready for plmc_generate_model)",
                "path": str(a2m_file.resolve())
            }
        ]
    }
