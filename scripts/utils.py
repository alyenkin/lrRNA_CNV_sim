import pysam
from typing import NamedTuple, List, Tuple


class Exon(NamedTuple):
    chrom_ref: str
    start_ref: int
    end_ref: int
    start_gene: int
    end_gene: int

    

class Transcript:
    exons: List[Exon]

    def __init__(self, exons: List[Exon]):
        self.exons = list(sorted(exons, key = lambda e: e.start_gene))

    def __str__(self):
        return str(self.exons)

    def __eq__(self, other):
        return self.exons == other.exons

    def split_exons(self, bp: int) -> Tuple[List[Exon], Exon | None, List[Exon]]:
        before = []
        interrupt = None
        after = []
        for e in self.exons:
            if e.end_ref < bp:
                before.append(e)
            elif e.start_ref > bp:
                after.append(e)
            else:
                interrupt = e
        return before, interrupt, after


    def split_exons_2(self, bp1: int, bp2: int) -> Tuple[List[Exon], Exon | None, List[Exon], Exon | None, List[Exon]]:
        bp1_before, bp1_interrupt, bp1_after = self.split_exons(bp1)
        bp2_before, bp2_interrupt, bp2_after = self.split_exons(bp2)

        dup_before = bp1_before
        dup_after = bp2_after

        dup_internal = [e for e in bp1_after if e in bp2_before]

        if bp1_interrupt == bp2_interrupt and bp1_interrupt is not None and len(dup_internal) > 0:
            raise Exception('Something is wrong')
        
        return dup_before, bp1_interrupt, dup_internal, bp2_interrupt, dup_after

    def getsequence(self, reference: pysam.FastaFile) -> str:
        return ''.join([reference.fetch(e.chrom_ref, e.start_ref, e.end_ref) for e in self.exons])


def duplication(t: Transcript, bp1: int, bp2: int) -> Transcript:
    if (bp1 > bp2):
        bp1, bp2 = bp2, bp1

    dist = bp2 - bp1

    dup_before, bp1_interrupt, dup_internal, bp2_interrupt, dup_after = t.split_exons_2(bp1, bp2)

    new_exons = []
    new_exons.extend(dup_before)
    for e in dup_internal:
        new_exons.append(e)

        new_exons.append(Exon(
            chrom_ref=e.chrom_ref, start_ref=e.start_ref, end_ref=e.end_ref, start_gene=e.start_gene + dist, end_gene = e.end_gene + dist
        ))

    if bp1_interrupt is not None and bp2_interrupt is not None and bp1_interrupt == bp2_interrupt:
        e = bp1_interrupt
        new_exons.append(Exon(
            chrom_ref=e.chrom_ref, start_ref=e.start_ref, end_ref=bp2, start_gene=e.start_gene, end_gene = bp2-e.start_ref + e.start_gene
        ))
        new_exons.append(Exon(
            chrom_ref=e.chrom_ref, start_ref=bp1, end_ref=e.end_ref, start_gene=bp2-e.start_ref + e.start_gene, end_gene = e.end_gene + dist
        ))
    else:
        if bp1_interrupt is not None:
            e = bp1_interrupt
            new_exons.append(e)
            new_exons.append(Exon(
                chrom_ref=e.chrom_ref, start_ref=bp1, end_ref=e.end_ref, start_gene=e.end_gene - (e.end_ref-bp1) + dist, end_gene = e.end_gene + dist
            ))

        if bp2_interrupt is not None:
            e = bp2_interrupt
            new_exons.append(Exon(
                chrom_ref = e.chrom_ref, start_ref=e.start_ref, end_ref=bp2, start_gene=e.start_gene, end_gene = bp2-e.start_ref + e.start_gene
            ))
            new_exons.append(Exon(
                chrom_ref = e.chrom_ref, start_ref=e.start_ref, end_ref=e.end_ref, start_gene=e.start_gene + dist, end_gene = e.end_gene + dist
            ))
    for e in dup_after:
        new_exons.append(Exon(
            chrom_ref=e.chrom_ref, start_ref=e.start_ref, end_ref=e.end_ref, start_gene=e.start_gene + dist, end_gene = e.end_gene + dist
        ))
    return Transcript(new_exons)







