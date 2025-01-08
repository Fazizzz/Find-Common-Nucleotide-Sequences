# Find common gRNA/nucleotide sequences between sets 

This python script is a simple utility designed to allow for accurate identification of common gRNAs or any nucleotide sequences between sets of sequences. Often we have large sets of guides or probes and need to compare them for commonality based on the sequences alone. This python utility fits into that niche by allowing the user to provide large sets of sequences as individual columns in a single TSV file. Once run, this utility will output the sequences that are common between all the provided sets in the input file to a tsv file. The script additionally lets the user truncate the sequences in the 5`- 3` (left to right) direction to restrict the comparisons to a fraction of the sequence length. In the case gRNA’s this can be used to ignore the pam site at the end of the sequence. 

There is currently a reverse complement mode under development for this tool which allows the user to identify any reverse complement matches for sequences between the sets. The steps for that are currently to reverse complement the sequences and then truncate them which may cause issues. If using –-mode reverse you should avoid truncation and use the full length of your sequences while this feature is in development. It may also be helpful to remove pam sites before running this tool in  –-mode reverse to avoid errors in truncation while working with reverse complements. 

_______________________________________________________
## Table of Contents
- [Installation](#installation)
- [Requirements](#requirements)
- [Usage](#usage)
- [Input File Formats](#input-file-formats)
- [Output Format](#output-format)
- [Example](#example)
- [License](#license)
- [Acknowledgments](#acknowledgments)

_______________________________________________________
## Installation
Clone this repository to your local machine:

```bash
git clone https://github.com/yourusername/repo-name.git
cd repo-name
```
________________________________________________________
## Requirements

* `Python 3.8+`
* `Tab separated text file with sequences in rows provided under sets as columns. The script does expect your input to have a heading`

_______________________________________________________ 
## Usage

This script was intended to be used for gRNAs and can work with millions of sequences in a single file. It should also be possible to use this for probes or any set of nucleotide sequences. The truncation functionality is optional but defaults to 20bp because that’s the typical length of a gRNA sequence. Notably, this script is not case sensitive so if soft or hard masking are important factors in your analysis please keep in mind that the matching will convert everything to lowercase to ignore variation between sets provided. 

Run the script in base mode with:

```

python Find-Common-Sequences.py  -i {input.tsv} -o {output.tsv} -g {Truncate Length for comparison (default: 20)}

```

Run the script in reverse complement mode with:

```
python Find-Common-Sequences.py  -i {input.tsv} -o {output.tsv} -g {Full Length (default: 20)} --mode reverse

```

*Note: The reverse complement mode is still in development, so if you wish to use this mode please use the full length of your sequence or truncate your sequence manually before running the script. The current itteration of the script performs the reverse complement but then still truncates left to right which switches to 3`-5`. 


The script takes four main arguments:

*	`-h or --help:` Show this help message and exit
*	`-i or --input_file:` Path to TSV file with sequences.
*	`-o or --output_file:` Path to output TSV file with common sequences.
*	`-g or --guide_length` Desired length 5`-3` (left to right) for matching (default: 20bp). 
*	`--mode {base,reverse}:` Mode for processing: 'base' (default) or 'reverse' for reverse complement matching.
	Note: The input requires the columns to have headers for each set to function properly.

________________________________________________________

## Input File Format

# Example input file for comparison 

```
Set1	Set2
CTAATTTACAAACTGAACGAAGG	CTAATTTACAAACTGAACTTCGG
CTAATTTACAAACTGAACGACGG	CCTTCGTTCAGTTTGTAAATTAG
CTAATTTACAAACTGAACGTCGG	GGGGGGGGGGGGGGGGGGGGGGG
CTAATTTACAAACTGAACTTCGG	TTAATTTACAAACTGAACGAAGG
CTAATTTACAAACTGAACGAAGG	TTAATTTACAAACTGAACGAAGG
CCTTCGTTCAGTTTGTAAATTAA	TTAATTTACAAACTGAACGTAGG

```

_________________________________________________________

## Output Format

# Example output from base mode with default params

```
Common Guides
ccttcgttcagtttgtaaattag
ccttcgttcagtttgtaaattaa
ctaatttacaaactgaacttcgg

```


___________________________________________________________
## Example

Here's an example command for running the script run using the test data shown above:

```
python Find-Common-Sequences.py -i Test-file2.tsv -o Test-file2-output.tsv -g 20

```
___________________________________________________________

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](https://www.gnu.org/licenses/gpl-3.0.en.html#license-text) file for details.


___________________________________________________________

## Acknowledgments

* **M.Faizan Khalid** - *Author and current maintainer of the code*

This script was developed by Muhammad Faizan Khalid for use in CRISPR gRNA design and analysis. The script is intended to be a utility for use in bioinformatics analysis and analysis pipelines. Usage and implementation of this script is currently free, and this is not a commercial product. The author is not responsible for damages, upkeep and any issues with the script. It is provided as is. Please report any bugs or issues for future improvements which may be added into the repository. 
  
For citing this tool, please use Khalid M.Faizan or Khalid MF. You can follow my research using my [Google Scholar profile](https://scholar.google.com/citations?hl=en&user=qFZQ5wYAAAAJ&sortby=title&view_op=list_works&gmla=AL3_zigRWGX9g8Jc22idbBUMFuy7cVN_pEIyL6_DXSA-qWkJbcaONzhRNSmAwmQXKEm-3-WYGouZZC2pCE6zD9tZLxizbM7jQzzZMOgtkgsuL825u4lvSs9kwsccajhJbBg2Mrc37at_HCQ).

This project is made possible thanks to the open-source bioinformatics community for their resources and support.

