from pathlib import Path
import utils
from utils import Transcript, Exon
import pysam
import subprocess
from typing import Dict, Tuple

def import_gtf(gtf_file) -> Transcript:
    exons = []
    inf = open(gtf_file, 'r')
    transcript_line = inf.readline()
    transcript_line_splits = transcript_line.split('\t')
    if transcript_line_splits[2] != 'transcript':
        raise Exception('First line is not a transcript')
    transcript_start = int(transcript_line_splits[3])
    for line in inf:
        s = line.split('\t')
        if s[2] != 'exon':
            continue
        exons.append(Exon(
            chrom_ref=s[0],
            start_ref=int(s[3]),
            end_ref=int(s[4]),
            start_gene=int(s[3]) - transcript_start,
            end_gene=int(s[4]) - transcript_start
        ))
    inf.close()
    return Transcript(exons)


def generate_fastq_file(filename: str, sequence: str, n_reads: int):
    with open(filename, 'w') as outf:
        for i in range(n_reads):
            outf.write(f'@read{i}\n{sequence}\n+\n{"I"*len(sequence)}\n') 


def import_bp_bed(bed_file: str) -> Dict[str, Tuple[int, int]]:
    bps = dict()
    with open(bed_file, 'r') as inf:
        for line in inf:
            s = line.strip().split('\t')
            bps[s[3]] = (int(s[1]), int(s[2]))
    return bps

def generate_dup_ref(ref: pysam.FastaFile, chrom: str, start: int, end: int, bp1: int, bp2: int) -> str:
    if (bp1 >= bp2):
        raise Exception('Breakpoint 1 must be less than breakpoint 2')
    if (start >= end):
        raise Exception('Start must be less than end')
    if (bp1 < start) or (bp2 > end):
        raise Exception('Breakpoint is outside of the duplication region')

    seq1 = ref.fetch(chrom, start, bp1)
    seq2 = ref.fetch(chrom, bp2, end)
    dupseq = ref.fetch(chrom, bp1, bp2)
    return seq1 + dupseq + dupseq + seq2


if __name__ == '__main__':
    ref_file = 'ref/hg38_chr7.fa'
    muc3a = import_gtf('ref/MUC3A.gtf')
    ref = pysam.FastaFile(ref_file)

    dup_bps = import_bp_bed('input/dup_bps.bed')

    for k, v in dup_bps.items():
        dup = utils.transcript_duplication(muc3a, v[0], v[1])
        generate_fastq_file(f'output/fastq/{k}.fastq', dup.getsequence(ref), 100)

    fastq_files = Path('output/fastq').glob('*.fastq')
    for fq in fastq_files:
        output_sam = f'output/alignments/{fq.stem}.sam'
        command = f"minimap2 -ax splice:hq -uf {ref_file} {str(fq)} > {output_sam}"
        subprocess.run(command, shell=True)

        command = f"samtools view -bS {output_sam} > {output_sam.replace('.sam', '.bam')}"
        subprocess.run(command)

    sqanti3_path = Path("~/programs/SQANTI/SQANTI3-5.2.2")

    bam_files = Path('output/alignments').glob('*.bam')
    for bam in bam_files:
        command = f"python {sqanti3_path}/sqanti3_qc.py -a {bam} -o {bam.stem}"
        subprocess.run(command, shell=True)
        





