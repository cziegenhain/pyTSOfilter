# pyTSOfilter
python script to filter TSO artifact reads from aligned bam files
Current version: v.0.1.2

Important note: Only zUMIs-processed bam files are compatible. Assignment of UMI-reads to genes must be performed in stranded mode during zUMIs processing.

## Requirements
- pysam
- regex

## Usage
```
usage: pyTSOfilter.py [-h] --bam FILENAME --out FILENAME --fa FILENAME [--p P]
                      [--n_mismatch N_MISMATCH] [--ggg]

optional arguments:
  -h, --help            show this help message and exit
  --bam FILENAME        Path to input BAM file
  --out FILENAME        Path to output bam file
  --fa FILENAME         Path to genome reference fasta
  --p P                 Number of processes to use
  --n_mismatch N_MISMATCH
                        Number of mismatches allowed when matching TSO
  --ggg                 Require template-switching Gs as part of sequence match
```

