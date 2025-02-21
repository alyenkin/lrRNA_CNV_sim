import cyvcf2
import pysam


def convert_record(variant: cyvcf2.Variant, ref: pysam.FastaFile) -> cyvcf2.Variant:
    var_class = variant.alt[0]
    if var_class not in ['<DEL>', '<DUP>']:
        return variant
    chrom = variant.CHROM
    pos = variant.POS

    end = variant.info['END']
    ref_seq = ref.fetch(chrom, pos - 1, pos)
    alt = ref.fetch(chrom, pos - 1, end)

    if variant.alt[0] == '<DEL>':
        ref_seq, alt = alt, ref_seq

    variant.REF = ref_seq
    variant.ALT = [alt]
    return variant

if __name__ == '__main__':
    in_vcf_file = ''
    out_vcf_file = ''

    in_vcf = cyvcf2.VCF(in_vcf_file)
    out_vcf = cyvcf2.Writer(out_vcf_file, in_vcf, mode = 'wz')

    ref_file = ''
    ref = pysam.FastaFile(ref_file)

    for record in in_vcf:
        out_vcf.write_record(convert_record(record, ref))

