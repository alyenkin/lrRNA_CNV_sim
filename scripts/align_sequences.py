from pathlib import Path
import subprocess

ref_file = 'ref/hg38_chr7.fa'

if __name__ == '__main__':
    fastq_files = Path('output/fastq').glob('*.fastq')
    for fq in fastq_files:
        output_sam = f'output/alignments/{fq.stem}.sam'
        command = f"minimap2 -ax splice:hq -uf {ref_file} {str(fq)} > {output_sam}"
        subprocess.run(command, shell=True)

        command = f"samtools view -bS {output_sam} > {output_sam.replace('.sam', '.bam')}"
        subprocess.run(command, shell=True)
