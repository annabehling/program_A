#!/usr/bin/env python

import sys #otherwise sys.exit gives sys is undefined
from Bio import SeqIO #SeqIO already has methods like convert and others
def fasta_file(fastq_file, outfile):
	"""convert fastq file to fasta file
	fastq_file : path to a .fq file (str)
	outfile : path to a .fna file that will be created (str)"""
	SeqIO.convert(fastq_file,'fastq',outfile,'fasta')

def variant_sites(fq_seqs, refseq): #to get sites that are variant, in case you want to do things other than make a vcf file
    """takes 2 sequences and returns the position of the variant site and the nucleotide variants
    fq_seqs : fastq sequences (iterable)
    refseq : path to the reference sequence file (str)
    returns index(1-based), variant, reference base (list of tuples)
    """
    var_sites = []
    var_indices = [] #keep track of sites that we've already found
    for rec in fq_seqs:
        for index,base in enumerate(rec): #enumerate goes over it keeping track of where we are
            if base != refseq[index]:
                if index not in var_indices:
                    var_sites.append( (index+1,base,refseq[index]) )
                    var_indices.append(index)
    return(var_sites)

def vars_to_vcf(var_sites, contig_name, outfile): #to write the set of variants to an outfile
    """takes a list of tuples and returns .vcf file with tab separated values of CHR, POS, REF, ALT, ID, FILTER, INFO
    var_sites : list of tuples
    contig_name : name of conttig on reference genome for vcf output
    outfile : path to newly created .vcf file (str)
    the columns ID, FILTER, INFO have been hardcoded"""
    outhandle = open(outfile, 'w')
    outhandle.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format("#CHROM", "POS", "REF", "ALT", "ID", "FILTER", "INFO"))
    for POS, REF, ALT in var_sites:
        outhandle.write('{}\t{}\t{}\t{}\t.\t.\t.\n'.format(contig_name, POS, REF, ALT))

if __name__ == '__main__': #only need this for command line executable. Not relevant for ipython notebook usage. only thing that gets called whehn you use the script
    try:
        fastq_file = sys.argv[1]
        base_name = sys.argv[2]
        reffile = sys.argv[3]
    except IndexError:
        print("Usage: ./SNP_finder.py [fastq_file] [outfile_base_name] [reffile]") #./ is there so if there's a program else where called SNP_finder.py, it will execute this one
        sys.exit(1)
    chromosome = SeqIO.read(reffile,format='fasta') #seqIO.parse to get any reference seq
    fq = SeqIO.parse(fastq_file,format='fastq') #seqIO.parse to get any fastq file
    result = variant_sites(fq, chromosome)
    vars_to_vcf(result, chromosome.id, base_name+'.vcf') #outputs .vcf file
    fasta_file(fastq_file, base_name+'.fna') #outputs .fna file

sys.exit(0)

#chmod u+x SNP_finder.py to make it executable. very last thing to do at the end of the script