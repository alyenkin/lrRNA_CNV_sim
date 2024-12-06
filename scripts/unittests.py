import pytest
import utils
from utils import Transcript, Exon


transcript_oneexon = Transcript([
    Exon('chr1', 0, 100, 0, 100)
])

transcript_simple = Transcript([
    Exon('chr1', 2000, 2100, 0, 100),
    Exon('chr1', 2300, 2400, 300, 400),
    Exon('chr1', 2500, 2600, 500, 600)
])

transcript_complex = Transcript([
    Exon('chr1', 2000, 2100, 0, 100),
    Exon('chr1', 2300, 2400, 300, 400),
    Exon('chr1', 2500, 2600, 500, 600),
    Exon('chr1', 2700, 2800, 700, 800),
    Exon('chr1', 2900, 3000, 900, 1000)
])

def test_bp_simple():
    splits = transcript_oneexon.split_exons(50)
    assert splits[0] == []
    assert splits[1] == Exon('chr1', 0, 100, 0, 100)
    assert splits[2] == []

def test_bp_before():
    splits = transcript_simple.split_exons(150)
    assert splits[0] == []
    assert splits[1] is None
    assert splits[2] == transcript_simple.exons

def test_bp2_simple():
    splits = transcript_oneexon.split_exons_2(50, 60)
    assert splits[0] == []
    assert splits[1] == Exon('chr1', 0, 100, 0, 100)
    assert splits[2] == []
    assert splits[3] == Exon('chr1', 0, 100, 0, 100)
    assert splits[4] == []

def test_bp2_before():
    splits = transcript_simple.split_exons_2(150, 250)
    assert splits[0] == []
    assert splits[1] is None
    assert splits[2] == []
    assert splits[3] is None
    assert splits[4] == transcript_simple.exons

def test_bp2_ied():
    splits = transcript_simple.split_exons_2(2150, 2450)
    assert splits[0] == [Exon('chr1', 2000, 2100, 0, 100)]
    assert splits[1] is None
    assert splits[2] == [Exon('chr1', 2300, 2400, 300, 400)]
    assert splits[3] is None
    assert splits[4] == [Exon('chr1', 2500, 2600, 500, 600)]


def assert_transcript(t1, t2):
    assert len(t1.exons) == len(t2.exons), f'Lengths differ: {len(t1.exons)} != {len(t2.exons)}'
    for e1, e2 in zip(t1.exons, t2.exons):
        assert e1 == e2, f'{e1} != {e2}'

def test_dup_basic():
    dup = utils.duplication(transcript_oneexon, 50, 60)
    expected = Transcript([
        Exon('chr1', 0, 60, 0, 60),
        Exon('chr1', 50, 100, 60, 110)
    ])
    assert_transcript(dup, expected)

def test_dup_intronic():
    dup = utils.duplication(transcript_simple, 2150, 2250)
    expected = Transcript([
        Exon('chr1', 2000, 2100, 0, 100),
        Exon('chr1', 2300, 2400, 400, 500),
        Exon('chr1', 2500, 2600, 600, 700)
    ])
    assert_transcript(dup, expected)

def test_dup_interexon():
    dup = utils.duplication(transcript_complex, 2350, 2750)
    expected = Transcript([
        Exon('chr1', 2000, 2100, 0, 100),
        Exon('chr1', 2300, 2400, 300, 400),
        Exon('chr1', 2500, 2600, 500, 600),
        Exon('chr1', 2700, 2750, 700, 750),
        Exon('chr1', 2350, 2400, 750, 800),
        Exon('chr1', 2500, 2600, 900, 1000),
        Exon('chr1', 2700, 2800, 1100, 1200),
        Exon('chr1', 2900, 3000, 1300, 1400)
    ])


