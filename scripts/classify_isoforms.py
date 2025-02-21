from pathlib import Path
import subprocess


sqanti3_path = "~/programs/SQANTI3-5.2.2"
stringtie_path = "~/programs/stringtie/stringtie"

full_ref_gtf = "ref/MUC3A_all.gtf"
ref_fa = "ref/hg38_chr7.fa"
if __name__ == '__main__':
    bam_files = Path('output/alignments').glob('*.bam')
    for bam in bam_files:
        assemble_file = Path(f"output/classify/assembled_transcriptome/{bam.stem}.gtf")
        command = f"""
            {stringtie_path} {bam} -L -o {assemble_file}
        """
        subprocess.run(command, shell=True)

        command = f"""
            python {sqanti3_path}/sqanti3_qc.py {assemble_file} \
                {full_ref_gtf} {ref_fa} \
                -o {bam.stem} -d output/classify/classified_isoforms
        """
        subprocess.run(command, shell=True)
