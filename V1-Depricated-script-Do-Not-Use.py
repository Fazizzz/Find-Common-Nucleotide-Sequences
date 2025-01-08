"""
Find-Matching-gRNA-Seq.py

This script identifies common guide RNA (gRNA) sequences across multiple columns in a dataset.
It can handle sequence truncation for initial matching and has an optional second pass to search for reverse complement matches.

Usage:
    python Find-Matching-gRNA-Seq.py -i <input_file> -o <output_file> [-g <guide_length>] [--look_for_reverse_complement]

Requirements:
    - Python 3
    - Pandas library for data manipulation
    - Biopython for handling reverse complement sequences (only needed if using the reverse complement option)

Example:
    python Find-Matching-gRNA-Seq.py -i input_data.tsv -o output_common_guides.tsv -g 20 --look_for_reverse_complement

Author: M. Faizan Khalid
"""

import pandas as pd
import argparse
from Bio.Seq import Seq  # Biopython library for reverse complement functionality


def find_common_guides(input_file, output_file, guide_length=20, look_for_reverse_complement=False):
    """
    Identifies common guides across multiple columns in a TSV file, with optional sequence truncation and
    a second pass to search for reverse complement matches within each column.

    Parameters:
        input_file (str): Path to the input TSV file containing guide sequences across columns.
        output_file (str): Path to save the output TSV file containing common guide sequences.
        guide_length (int): Length to which each guide is truncated for matching (default: 20bp).
        look_for_reverse_complement (bool): Whether to perform a second pass for reverse complement matching.
    """
    try:
        # Read input file as a DataFrame, allowing for rows with irregular numbers of columns
        df = pd.read_csv(input_file, sep='\t', dtype=str, low_memory=False, on_bad_lines='skip')
        print("Column names:", df.columns)  # Debugging step to verify columns in the dataset

        # Convert all values to lowercase for case-insensitive matching
        df = df.applymap(lambda x: x.lower() if isinstance(x, str) else '')

        # First Pass: Truncate guides and find direct common matches across columns
        df_truncated = df.applymap(lambda x: x[:guide_length] if isinstance(x, str) and x != 'na' else '')
        common_guides = set(df_truncated.iloc[:, 0].dropna())
        
        # Intersect common guides across all columns
        for col in df_truncated.columns[1:]:
            common_guides.intersection_update(set(df_truncated[col].dropna()))

        # Collect unique full-length guides matching the common truncated guides
        unique_common_guides = set(
            guide for col in df.columns for guide in df[col] if guide[:guide_length] in common_guides
        )
        
        # Initialize output with unique common guides
        output_data = {'Common Guides': list(unique_common_guides)}

        # Second Pass: Reverse Complement Matching if flag is enabled
        reverse_complement_matches = set()
        if look_for_reverse_complement:
            for col in df.columns:
                # Generate reverse complements, then truncate in the 3'-5' direction
                truncated_reverse_complements = df[col].apply(
                    lambda x: str(Seq(x).reverse_complement())[:guide_length] if isinstance(x, str) and x != 'na' else ''
                )
                
                # Check each truncated reverse complement against truncated guides in all other columns
                for other_col in df_truncated.columns:
                    if other_col != col:
                        reverse_complement_matches.update(
                            guide for guide, rev_comp in zip(df[col], truncated_reverse_complements) 
                            if rev_comp in set(df_truncated[other_col])
                        )
            
            # Add reverse complements to the output data, ensuring uniqueness
            output_data['Reverse Complements'] = list(reverse_complement_matches)

        # Align list lengths for DataFrame creation
        max_length = max(len(output_data['Common Guides']), len(output_data.get('Reverse Complements', [])))
        output_data['Common Guides'] += [''] * (max_length - len(output_data['Common Guides']))
        if 'Reverse Complements' in output_data:
            output_data['Reverse Complements'] += [''] * (max_length - len(output_data['Reverse Complements']))

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
    parser.add_argument('--look_for_reverse_complement', action='store_true',
                        help="Enable second pass to search for reverse complements (default: False)")

    args = parser.parse_args()

    # Execute the function with the specified arguments
    find_common_guides(
        args.input_file,
        args.output_file,
        args.guide_length,
        args.look_for_reverse_complement
    )


# Enable script execution from the command line
if __name__ == "__main__":
    main()

