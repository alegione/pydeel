'''
Module      : Convert2fasta
Description : Tool to convert input files to fasta.
Copyright   : (c) Alistair Legione, 11 Jul 2019
License     : MIT
Maintainer  : legionea@unimelb.edu.au
Portability : POSIX

'''

from Bio import SeqIO

def gbkToFasta(gbk_filename, fasta_filename):

    input_handle  = open(gbk_filename, "r")

    output_handle = open(fasta_filename, "w")

    for seq_record in SeqIO.parse(input_handle, "genbank") :
        print "Converting GenBank record %s to FASTA" % seq_record.id
        output_handle.write(">%s %s\n%s\n" % (
               seq_record.id,
               seq_record.description,
               seq_record.seq))

    output_handle.close()
    input_handle.close()

    return None

def gbkTofaa(gbk_filename, faa_filename):

    input_handle  = open(gbk_filename, "r")

    output_handle = open(faa_filename, "w")

    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            output_handle.write(">%s from %s\n%s\n" % (
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_record.name,
                   seq_feature.qualifiers['translation'][0]))


    output_handle.close()
    input_handle.close()

    return None


def main():


if __name__ == '__main__':
    main()
