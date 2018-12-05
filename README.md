# SNP finder

SNP_finder.py reads in a FASTQ file containing aligned reads (.fq) and a reference sequence file (.fna), and outputs two results files: a file listing all reads, unmodified in FASTA format (.fna); and a file listing only the variable sites in VCF format (.vcf).

This program was developed to the specifications of the Program A assignment for Special Topic paper 247.782, 2018.

## Prerequisites

* Requires [Biopython](https://biopython.org/)   
* Clone this repository to your local machine using `https://github.com/annabehling/program_A`

## Usage

```
./SNP_finder.py [fastq_file] [outfile_base_name] [reffile]
```
Produces two files:
`outfile_base_name.fna` and `outfile_base_name.fq`

## Acknowledgements

* `fasta_file()` was adapted from the [Biopython wiki page](https://biopython.org/wiki/Converting_sequence_files)
* And thank you to [David Winter](https://github.com/dwinter/) for his help with developing this program

## License

* [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
