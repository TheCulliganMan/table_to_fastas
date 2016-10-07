#!/usr/bin/env python
from copy import deepcopy
from Bio import SeqIO
from Bio.Seq import Seq

def yield_positions(table_path):
    with open(table_path) as input_handle:
        for num, line in enumerate(input_handle):
            if num:
                split_line = line.split()
                yield split_line[:2]

def build_chrom_dict(table_path):
    positions = yield_positions(table_path)
    chrom_dict = {}
    for chrom, pos in positions:
        if chrom not in chrom_dict:
            chrom_dict[chrom] = []
        chrom_dict[chrom].append(int(pos))
    return chrom_dict

def yield_sequences(fasta_path, table_path, buff=50):
    chrom_dict = build_chrom_dict(table_path)
    with open(fasta_path) as input_handle:
        for record in SeqIO.parse(input_handle, 'fasta'):
            if record.id in chrom_dict:
                for pos in chrom_dict[record.id]:
                    if (pos - buff >= 0) and (pos + buff <= len(record.seq)):
                        new_record = deepcopy(record)
                        new_record.seq = new_record.seq[pos-buff:pos+buff]
                        new_record.id = "{}_{}".format(new_record.id, pos)
                        yield new_record

def write_sequences(fasta_path, table_path, fasta_output_path):
    sequence_generator = yield_sequences(fasta_path, table_path)
    with open(fasta_output_path, 'w+') as output_handle:
        SeqIO.write(sequence_generator, output_handle, 'fasta')

def main():
    table_path = "simus_variants_table_2.tsv"
    fasta_path = "MaSuRCA_contigs.fasta"
    output_path = "simus_output_array.fasta"
    write_sequences(fasta_path, table_path, output_path)






if __name__ == "__main__":
    main()
