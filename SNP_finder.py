import sys #otherwise sys.exit gives sys is undefined
from Bio import SeqIO
def fasta_file(fastq_file, outfile):
	"""convert fastq file to fasta file"""
	SeqIO.convert('infile','fastq','output','outfile')

def fastq_sequence(fastq_file, outfile):
    """takes a fastq file and writes all sequences to a new file
    infile : .fq file (str)
    outfile : new file containing output sequences"""
    with open (infile) as fastq_file:
        outhandle = open(outfile, "w")
        for i,line in enumerate(fastq_file):
            if i % 4 == 1:
                outhandle.write(line)

def variant_sites(fq_seqs, refseq): #to get sites that are variant
    """takes 2 sequences and returns the position of the variant site and the nucleotide variants
    fq_seqs : fastq sequences (iterable)
    refseq : reference sequence (str)
    returns index(1-based), variant, reference base
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

def vars_to_vcf(var_sites, chromosome, outfile): #to write the information to an outfile and get chromosome
    """takes a list of tuples and returns .vcf file with tab separated values of CHR, POS, REF, ALT, ID, FILTER, INFO
    var_sites : list of tuples
    chromosome : str
    the columns ID, FILTER, INFO have been hardcoded"""
    outhandle = open(outfile, 'w')
    outhandle.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format("#CHROM", "POS", "REF", "ALT", "ID", "FILTER", "INFO"))
    for POS, REF, ALT in var_sites:
        outhandle.write('{}\t{}\t{}\t{}\t.\t.\t.\n'.format(chromosome, POS, REF, ALT))

def vcf_file(fastq_file, reffile):
    """does overall stuff""" #call in def fasta_file , def fastq_sequence , def variant_sites , def_vars_to_vcf under this one function
    chromosome = SeqIO.read(reffile,format='fasta') #seqIO.parse to get any reference seq
    fq = SeqIO.parse(fastq_file,format='fastq') #seqIO.parse to get any fastq file
    result = variant_sites(fq, chromosome)
    return(result)
    #thing for outfile from end of january program
    #make executable using chmod u+x. very last thing to do at the end of the script
    #./ in readme

variable_sites = vcf_file("input.fq", "reference.fna") #include this bit
vars_to_vcf(variable_sites, "reference", "result.vcf") #and this bit

#only need this for command line executable. Not relevant for ipython notebook usage
if __name__ == '__main__':
    try:
        fastq_file = sys.argv[1]
        outfile = sys.argv[2]
        reffile = sys.argv[3]
        print('i got executed')
    except IndexError:
        print("Usage: ./program_A.py [fastq_file] [outfile] [reffile]")
        sys.exit(1)
    vcf_file(fastq_file, outfile, reffile)
sys.exit(0)