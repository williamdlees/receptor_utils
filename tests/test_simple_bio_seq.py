import sys
import warnings

import pytest
from Bio.Align import AlignInfo

from receptor_utils import simple_bio_seq as sbs


def test_chunks():
    assert list(sbs.chunks("ABCDEFG", 3)) == ["ABC", "DEF", "G"]
    assert list(sbs.chunks("", 3)) == []


def test_fasta_roundtrip_uppercases(tmp_path):
    path = str(tmp_path / "x.fasta")
    sbs.write_fasta(path, {"a": "acgt", "b": "ttGG"})
    assert sbs.read_fasta(path) == {"a": "ACGT", "b": "TTGG"}


def test_write_fasta_deprecated_argument_order(tmp_path):
    path = str(tmp_path / "x.fasta")
    with pytest.warns(DeprecationWarning):
        sbs.write_fasta({"a": "ACGT"}, path)
    assert sbs.read_fasta(path) == {"a": "ACGT"}


def test_read_single_fasta(tmp_path):
    path = tmp_path / "x.fasta"
    path.write_text(">a\nacg\n>b\nttt\n")
    assert sbs.read_single_fasta(str(path)) == "ACG"


def test_read_single_fasta_empty(tmp_path):
    path = tmp_path / "x.fasta"
    path.write_text("")
    assert sbs.read_single_fasta(str(path)) is None


IMGT = (
    ">X1|IGHV1-1*01|Homo sapiens|F|V-REGION|1..10|10 nt|\n" + "A" * 300 + "\n"
    ">X2|IGHV1-2*01|Homo sapiens|P|V-REGION|1..10|10 nt|\n" + "C" * 300 + "\n"
    ">X3|IGHV1-3*01|Homo sapiens|F|V-REGION|1..10|10 nt|\n" + "G" * 100 + "\n"
    ">X4|IGHV1-4/OR15*01|Homo sapiens|F|V-REGION|1..10|10 nt|\n" + "T" * 300 + "\n"
    ">X5|IGHJ4*01|Homo sapiens|F|J-REGION|1..10|10 nt|\n" + "ACGTACGT" + "\n"
    ">X6|IGHJ4*01|Mus musculus|F|J-REGION|1..10|10 nt|\n" + "TTTT" + "\n"
    ">X7|IGHG1*01|Homo sapiens|F|CH1|1..10|10 nt|\n" + "GGGG" + "\n"
)


@pytest.fixture
def imgt_file(tmp_path):
    p = tmp_path / "imgt.fasta"
    p.write_text(IMGT)
    return str(p)


def test_read_imgt_fasta_basic(imgt_file):
    res = sbs.read_imgt_fasta(imgt_file, ["Homo_sapiens".replace("_", " ")], chains=("IGHV", "IGHJ"))
    hs = res["Homo sapiens"]
    assert set(hs["IGHV"]) == {"IGHV1-1*01", "IGHV1-2*01", "IGHV1-3*01"}   # orphon excluded
    assert hs["IGHJ"] == {"IGHJ4*01": "ACGTACGT"}


def test_read_imgt_fasta_functional_only(imgt_file):
    res = sbs.read_imgt_fasta(imgt_file, ["Homo sapiens"], chains=("IGHV",), functional_only=True)
    # pseudogene and short V are dropped
    assert set(res["Homo sapiens"]["IGHV"]) == {"IGHV1-1*01"}


def test_read_imgt_fasta_orphon_and_species(imgt_file):
    res = sbs.read_imgt_fasta(imgt_file, ["Homo sapiens", "Mus musculus"], chains=("IGHV", "IGHJ"), include_orphon=True)
    assert "IGHV1-4/OR15*01" in res["Homo sapiens"]["IGHV"]
    assert res["Mus musculus"]["IGHJ"] == {"IGHJ4*01": "TTTT"}
    assert res["Mus musculus"]["IGHV"] == {}


def test_read_imgt_fasta_constant_names_get_seq_type(imgt_file):
    res = sbs.read_imgt_fasta(imgt_file, ["Homo sapiens"], chains=("CH",))
    assert res["Homo sapiens"]["CH"] == {"IGHG1_CH1*01": "GGGG"}


def test_sample_fasta_samples_subset():
    seqs = {str(i): "ACGT" for i in range(10)}
    sample = sbs.sample_fasta(seqs, 3)
    assert len(sample) == 3
    assert all(k in seqs for k in sample)


def test_sample_fasta_returns_all_when_not_enough():
    seqs = {"a": "A", "b": "C"}
    assert sbs.sample_fasta(seqs, 5) is seqs
    assert sbs.sample_fasta(seqs, 2) is seqs


def test_translate():
    assert sbs.translate("ATGAAATTT") == "MKF"
    assert sbs.translate("ATGAAATTTGG") == "MKF"                      # truncated
    assert sbs.translate("ATGAAATTTGG", truncate=False) == "MKFG"      # padded with N (GGN -> G)
    assert sbs.translate("ATG---TTT") == "M-F"


def test_translate_partial_codon_gaps():
    assert sbs.translate("ATGA-ATTT") == "M-F"
    assert sbs.translate("ATGA.ATTT") == "M-F"
    with pytest.raises(Exception):
        sbs.translate("ATGA-ATTT", ignore_partial_codon=False)


def test_reverse_complement():
    assert sbs.reverse_complement("AACG") == "CGTT"
    assert sbs.reverse_complement("") == ""


def test_nt_diff():
    assert sbs.nt_diff("ACGT", "ACGT") == 0
    assert sbs.nt_diff("ACGT", "ACCT") == 1
    assert sbs.nt_diff("acgt", "ACGT") == 0                 # case-insensitive
    assert sbs.nt_diff("ACNT", "ACGT") == 0                 # N ignored
    assert sbs.nt_diff("ACXT", "ACGT") == 0                 # X ignored
    assert sbs.nt_diff("ACGT", "ACGTAA") == 2               # length difference counts


def test_to_seq_records():
    recs = sbs.toSeqRecords({"a": "ACGT", "b": "TT"})
    assert [r.id for r in recs] == ["a", "b"]
    assert [str(r.seq) for r in recs] == ["ACGT", "TT"]


def test_dumb_consensus():
    seqs = {"1": "ACGT", "2": "ACGT", "3": "ACCT"}
    assert sbs.dumb_consensus(seqs, threshold=0.6) == "ACGT"
    assert sbs.dumb_consensus(seqs, threshold=0.9) == "ACXT"


def test_scored_consensus():
    # data from the former tests/test_scored_gap_consensus.py
    seqs = {"0": "abcdef", "1": "abcded", "2": "aacdef"}
    cons, score = sbs.scored_consensus(seqs)
    assert str(cons) == "aXcdeX"
    assert score == pytest.approx(2 / 3)


def test_scored_consensus_ambiguous_position():
    cons, score = sbs.scored_consensus({"0": "AAA", "1": "ACA"}, threshold=0.7)
    assert str(cons) == "AXA"
    assert score == pytest.approx(0.5)


def test_scored_gap_consensus_require_multiple():
    from Bio.Align import MultipleSeqAlignment
    msa = MultipleSeqAlignment(sbs.toSeqRecords({"0": "AAA", "1": "AAA"}))
    cons, score = sbs.scored_gap_consensus(msa, require_multiple=True)
    assert str(cons) == "AAA"
    assert score == 1.0


def test_closest_ref():
    ref = {"r1": "AAAA", "r2": "AAAT", "r3": "TTTT"}
    assert sbs.closest_ref("AAAA", ref) == ["r1"]
    assert sbs.closest_ref("AAAC", ref) == ["r1", "r2"]


def test_csv_roundtrip(tmp_path):
    path = str(tmp_path / "x.csv")
    rows = [{"a": "1", "b": "2"}, {"a": "3", "b": "4"}]
    sbs.write_csv(path, rows)
    assert sbs.read_csv(path) == rows


def test_csv_delimiter(tmp_path):
    path = str(tmp_path / "x.tsv")
    rows = [{"a": "1", "b": "2"}]
    sbs.write_csv(path, rows, delimiter="\t")
    assert "a\tb" in open(path).read()
    assert sbs.read_csv(path, delimiter="\t") == rows


def test_write_csv_scan_all(tmp_path):
    path = str(tmp_path / "x.csv")
    sbs.write_csv(path, [{"a": "1"}, {"a": "2", "b": "3"}], scan_all=True)
    assert sbs.read_csv(path) == [{"a": "1", "b": ""}, {"a": "2", "b": "3"}]


def test_write_csv_empty_creates_nothing(tmp_path):
    path = tmp_path / "x.csv"
    sbs.write_csv(str(path), [])
    assert not path.exists()
