#!/usr/bin/env python
from copy import deepcopy
from Bio import SeqIO
from Bio.Seq import Seq

#!/usr/bin/env python
import vcf

"""
vcf_file = "out.5simus.thinned.50k.maxmaf075hwe.vcf.filtered.vcf"
output_handle = open("table.tsv", "w+")
with open('vcf_filtered_out.vcf', 'w+') as output_handle:
    with open(vcf_file) as input_handle:
        vcf_reader = vcf.Reader(input_handle)
        vcf_writer = vcf.Writer(output_handle, vcf_reader)
        for num, record in enumerate(vcf_reader):
            if not num:
                output_handle.write(
                    "\t".join(["CHROM", "POS", "QUAL", "REF", "ALT"] + [sample.sample for sample in record.samples] + ["\n"])
                )
            samples_list =  [sample.gt_bases.replace("/","") for sample in record.samples]
            unique = set(samples_list)
            if len(unique) > 2:
                vcf_writer.write_record(record)
                output_handle.write(
                    "\t".join([str(i) for i in [record.CHROM, record.POS, record.QUAL, record.REF, ",".join([str(j) for j in record.ALT])]] + samples_list + ["\n"])
                )
"""

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
