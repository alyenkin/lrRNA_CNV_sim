from pathlib import Path
import subprocess


sqanti3_path = "~/programs/SQANTI/SQANTI3-5.2.2"
if __name__ == '__main__':


    bam_files = Path('output/alignments').glob('*.bam')
    for bam in bam_files:
        command = f"python {sqanti3_path}/sqanti3_qc.py -a {bam} -o {bam.stem}"
        subprocess.run(command, shell=True)
