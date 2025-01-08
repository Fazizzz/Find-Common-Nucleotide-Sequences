"""
Find-Matching-gRNA-Seq.py

This script identifies common guide RNA (gRNA) sequences across multiple columns in a dataset.
It supports two modes:
1. `base` (default): Match guides in the 5'-3' direction.
2. `reverse`: Truncate sequences first, then convert them to reverse complements, and match in the 5'-3' direction.

In `reverse` mode, the output includes two columns:
- `Common Guides`: The original full-length sequence.
- `Reverse Complements`: The reverse complement of the truncated sequence used in matching.

Usage:
    python Find-Matching-gRNA-Seq.py -i <input_file> -o <output_file> [-g <guide_length>] [--mode base|reverse]

Requirements:
    - Python 3
    - Pandas library for data manipulation
    - Biopython for handling reverse complement sequences.

Author: M. Faizan Khalid
"""

import pandas as pd
import argparse
from Bio.Seq import Seq  # Biopython library for reverse complement functionality


def find_common_guides(input_file, output_file, guide_length=20, mode="base"):
    """
    Identifies common guides across multiple columns in a TSV file, using the specified mode.

    Parameters:
        input_file (str): Path to the input TSV file containing guide sequences across columns.
        output_file (str): Path to save the output TSV file containing common guide sequences.
        guide_length (int): Length to which each guide is truncated for matching (default: 20bp).
        mode (str): Mode for processing guides: 'base' for 5'-3', 'reverse' for reverse complement matching.
    """
    try:
        # Read input file as a DataFrame, allowing for rows with irregular numbers of columns
        df = pd.read_csv(input_file, sep='\t', dtype=str, low_memory=False, on_bad_lines='skip')
        print("Column names:", df.columns)  # Debugging step to verify columns in the dataset

        # Convert all values to lowercase for case-insensitive matching
        df = df.applymap(lambda x: x.lower() if isinstance(x, str) else '')

        # Truncate sequences to retain bases starting from the 5' end
        df_truncated = df.applymap(lambda x: x[:guide_length] if isinstance(x, str) and x != 'na' else '')

        if mode == "reverse":
            # Convert truncated sequences to reverse complements
            df_reversed = df_truncated.applymap(
                lambda x: str(Seq(x).reverse_complement()) if isinstance(x, str) and x != 'na' else ''
            )

            # Find common truncated reverse complements across columns
            common_reverse_complements = set(df_reversed.iloc[:, 0].dropna())
            for col in df_reversed.columns[1:]:
                common_reverse_complements.intersection_update(set(df_reversed[col].dropna()))

            # Map reverse complements back to their original sequences
            reverse_mapping = {}
            for col in df.columns:
                for original, truncated, reverse in zip(df[col], df_truncated[col], df_reversed[col]):
                    if reverse in common_reverse_complements:
                        reverse_mapping[reverse] = original

            # Create output with original sequences and their reverse complements
            output_data = {
                "Common Guides": [reverse_mapping[rev_comp] for rev_comp in common_reverse_complements],
                "Reverse Complements": list(common_reverse_complements),
            }
        else:  # mode == "base"
            # Find common truncated sequences across columns
            common_guides = set(df_truncated.iloc[:, 0].dropna())
            for col in df_truncated.columns[1:]:
                common_guides.intersection_update(set(df_truncated[col].dropna()))

            # Collect unique full-length guides matching the common truncated guides
            unique_common_guides = set(
                guide for col in df.columns for guide in df[col] if guide[:guide_length] in common_guides
            )

            # Create output with unique common guides
            output_data = {"Common Guides": list(unique_common_guides)}

        # Convert output data to a DataFrame
        df_output = pd.DataFrame(output_data)

        # Save the output to a TSV file
        df_output.to_csv(output_file, sep='\t', index=False)
        print(f"Output saved to {output_file}")

    except pd.errors.ParserError as e:
        print(f"ParserError: {e}")
        print("Error reading the file. Ensure that it is properly formatted as a TSV.")
    except Exception as e:
        print(f"Error: {e}")
        print("An error occurred while processing the guides. Please check the input file and parameters.")


def main():
    """
    Main function to parse command-line arguments and call the guide-finding function.
    """
    # Set up argument parsing for command-line usage
    parser = argparse.ArgumentParser(description="Find common guides across sites post-filtering.")
    parser.add_argument('-i', '--input_file', type=str, required=True, help="Path to the input TSV file")
    parser.add_argument('-o', '--output_file', type=str, required=True, help="Path to save the output TSV file")
    parser.add_argument('-g', '--guide_length', type=int, default=20, help="Guide length for matching (default: 20bp)")
    parser.add_argument('--mode', type=str, choices=['base', 'reverse'], default='base',
                        help="Mode for processing guides: 'base' (default) or 'reverse' for reverse complement matching.")

    args = parser.parse_args()

    # Execute the function with the specified arguments
    find_common_guides(
        args.input_file,
        args.output_file,
        args.guide_length,
        args.mode
    )


# Enable script execution from the command line
if __name__ == "__main__":
    main()

