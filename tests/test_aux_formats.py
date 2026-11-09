import json

import pytest

from receptor_utils import aux_formats as af
from receptor_utils import simple_bio_seq as simple

J4 = "ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
V_SEQ = "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGG" + "A" * 60


def delineation(**kw):
    d = {
        "delineation_scheme": "IMGT",
        "fwr1_start": 1, "fwr1_end": 75, "cdr1_start": 76, "cdr1_end": 99,
        "fwr2_start": 100, "fwr2_end": 150, "cdr2_start": 151, "cdr2_end": 174,
        "fwr3_start": 175, "fwr3_end": 288, "cdr3_start": 289,
    }
    d.update(kw)
    return d


def germline_set(alleles, wrap="list"):
    gs = {"locus": "IGH", "allele_descriptions": alleles}
    if wrap == "list":
        return {"GermlineSet": [gs]}
    if wrap == "items":
        return {"GermlineSet": {"items": [gs]}}
    return {"GermlineSet": gs}


def v_allele(label="IGHV1-1*01", **kw):
    a = {
        "label": label, "sequence_type": "V", "sequence": V_SEQ, "functional": True,
        "gene_start": 1, "gene_end": len(V_SEQ),
        "v_gene_delineations": [delineation()],
        "leader_1_start": 1, "leader_1_end": 6, "leader_2_start": 7, "leader_2_end": 15,
    }
    a.update(kw)
    return a


def j_allele(label="IGHJ4*01", **kw):
    a = {
        "label": label, "sequence_type": "J", "sequence": J4, "functional": True,
        "gene_start": 1, "gene_end": len(J4), "j_cdr3_end": 19, "j_codon_frame": 1,
    }
    a.update(kw)
    return a


def d_allele(label="IGHD1-1*01"):
    return {"label": label, "sequence_type": "D", "sequence": "GGTACAACTGGAACGAC", "gene_start": 1, "gene_end": 17}


# ---- ndm ----

def test_record_from_json():
    rec = af.record_from_json("IGHV1-1*01", delineation(cdr1_end=99), "VH")
    assert rec == "IGHV1-1*01\t1\t75\t76\t99\t100\t150\t151\t174\t175\t288\tVH\t0\n"


@pytest.mark.parametrize("wrap", ["list", "items", "dict"])
def test_ndm_from_json_wrappers(tmp_path, wrap):
    out = tmp_path / "x.ndm"
    assert af.ndm_from_json(germline_set([v_allele()], wrap), str(out), "VH") is True
    lines = out.read_text().splitlines()
    assert lines[0].startswith("#name\tFWR1 start")
    assert lines[1] == "IGHV1-1*01\t1\t75\t76\t99\t100\t150\t151\t174\t175\t288\tVH\t0"


def test_ndm_from_json_no_v_genes(tmp_path, capsys):
    out = tmp_path / "x.ndm"
    assert af.ndm_from_json(germline_set([j_allele()]), str(out), "VH") is False
    assert not out.exists()
    assert "No V gene alleles" in capsys.readouterr().out


def test_ndm_from_json_omits_alleles_without_delineations(tmp_path, capsys):
    out = tmp_path / "x.ndm"
    a = v_allele("IGHV1-2*01")
    del a["v_gene_delineations"]
    af.ndm_from_json(germline_set([v_allele(), a]), str(out), "VH")
    assert "IGHV1-2*01" not in out.read_text()
    assert "no delineations" in capsys.readouterr().out


def test_ndm_from_json_adds_fake_allele(tmp_path, capsys):
    out = tmp_path / "x.ndm"
    af.ndm_from_json(germline_set([v_allele("IGHV1-3")]), str(out), "VH")
    assert out.read_text().splitlines()[1].startswith("IGHV1-3*00\t")
    assert "fake allele" in capsys.readouterr().out


def test_ndm_from_json_ignores_other_delineation_scheme(tmp_path):
    out = tmp_path / "x.ndm"
    a = v_allele(v_gene_delineations=[delineation(delineation_scheme="Kabat")])
    af.ndm_from_json(germline_set([a]), str(out), "VH")
    assert len(out.read_text().splitlines()) == 1       # header only


def test_ndm_from_json_matches_sample_file(tmp_path, data_dir):
    ref = json.load(open(data_dir / "Human_IGH_VDJ_rev_2.json"))
    out = tmp_path / "x.ndm"
    assert af.ndm_from_json(ref, str(out), "VH")
    got = [ln for ln in out.read_text().splitlines() if not ln.startswith("#")]
    expected = [ln for ln in (data_dir / "Human_IGH_VDJ_rev_2.ndm").read_text().splitlines() if ln and not ln.startswith("#")]
    assert len(got) == len(expected) > 0
    assert got[0] == expected[0]
    assert set(got) == set(expected)


def test_record_from_gapped_seq():
    seq = "A" * 10 + "...." + "C" * 10          # gapped: 14 positions + 10
    rec = af.record_from_gapped_seq("S*01", seq, "VH", (5, 8, 11, 14, 20))
    fields = rec.rstrip("\n").split("\t")
    assert fields[0] == "S*01"
    assert fields[-2:] == ["VH", "0"]
    assert fields[1] == "1"
    assert len(fields) == 13
    # coordinates are ungapped positions within the sequence
    nums = [int(x) for x in fields[1:11]]
    assert max(nums) <= len(seq.replace(".", ""))


def test_ndm_from_fasta_omits_truncated(tmp_path, write_fasta_file, capsys):
    full = "A" * 30
    fasta = write_fasta_file("g.fasta", {"IGHV1-1*01": full, "IGHV1-2*01": "." * 3 + "A" * 27})
    out = tmp_path / "x.ndm"
    af.ndm_from_fasta(fasta, str(out), "VH", (5, 10, 15, 20, 28))
    text = out.read_text()
    assert "IGHV1-1*01" in text and "IGHV1-2*01" not in text
    assert "truncated at the 5' end" in capsys.readouterr().out


# ---- aux ----

@pytest.mark.parametrize("wrap", ["list", "items", "dict"])
def test_aux_from_json(tmp_path, wrap):
    out = tmp_path / "x.aux"
    assert af.aux_from_json(germline_set([j_allele(), v_allele()], wrap), str(out), False) is True
    rows = simple.read_csv(str(out), delimiter="\t")
    assert rows == [{"#name": "IGHJ4*01", "j_codon_frame": "0", "chain_type": "JH", "j_cdr3_end": "17", "extra_bps": "0"}]


def test_aux_from_json_skips_bad_alleles(tmp_path, capsys):
    out = tmp_path / "x.aux"
    off_codon = j_allele("IGHJ1*01", j_cdr3_end=20)
    missing = j_allele("IGHJ2*01", j_codon_frame=None)
    ok = j_allele("IGHJ3*01")
    assert af.aux_from_json(germline_set([off_codon, missing, ok]), str(out), False)
    names = [r["#name"] for r in simple.read_csv(str(out), delimiter="\t")]
    assert names == ["IGHJ3*01"]
    captured = capsys.readouterr().out
    assert "not on a codon boundary" in captured and "Unrecognisable parameter" in captured


def test_aux_from_json_no_j_returns_false(tmp_path):
    out = tmp_path / "x.aux"
    assert af.aux_from_json(germline_set([v_allele()]), str(out), False) is False
    assert not out.exists()


def test_aux_from_json_fake_allele_and_sorting(tmp_path, capsys):
    out = tmp_path / "x.aux"
    af.aux_from_json(germline_set([j_allele("IGHJ6*01"), j_allele("IGHJ4")]), str(out), False)
    names = [r["#name"] for r in simple.read_csv(str(out), delimiter="\t")]
    assert names == ["IGHJ4*00", "IGHJ6*01"]
    assert "fake allele" in capsys.readouterr().out


def test_results_from_seqs_canonical_motif():
    (res,) = af.results_from_seqs({"IGHJ4*01": J4})
    assert res.motif == "WGQG"
    assert res.canonical_junc is True
    assert res.frame == 3
    assert res.notes == ""
    assert res.stop_codon is False or res.stop_codon is True


def test_results_from_seqs_non_canonical_motifs():
    # CGxG is accepted only as a non-canonical fallback
    seq = "AAATGTGGTAAAGGT" + "ACC"
    (res,) = af.results_from_seqs({"TRJP1*01": seq})
    assert res.canonical_junc is False
    assert "non-canonical motif" in res.notes


def test_results_from_seqs_no_motif():
    (res,) = af.results_from_seqs({"IGHJ9*01": "AAAAAAAAAAAA"})
    assert res.notes == "No j motif found"


def test_aux_from_seqs_writes_annotations(tmp_path, capsys):
    out = tmp_path / "x.aux"
    af.aux_from_seqs({"IGHJ4*01": J4, "IGHJ9*01": "AAAAAAAAAAAA"}, str(out), True)
    rows = simple.read_csv(str(out), delimiter="\t")
    assert [r["#name"] for r in rows] == ["IGHJ4*01"]
    assert rows[0]["chain_type"] == "JH"
    assert rows[0]["j_codon_frame"] == "2"
    printed = capsys.readouterr().out
    assert "MiAIRR" in printed and "no solutions" in printed


def test_aux_from_fasta(tmp_path, write_fasta_file):
    fasta = write_fasta_file("j.fasta", {"IGHJ4*01": J4})
    out = tmp_path / "x.aux"
    af.aux_from_fasta(fasta, str(out), False)
    assert simple.read_csv(str(out), delimiter="\t")[0]["#name"] == "IGHJ4*01"


def test_aux_from_fasta_sample_file(tmp_path, data_dir):
    out = tmp_path / "x.aux"
    seqs = simple.read_fasta(str(data_dir / "jh_seqs.fasta"))
    af.aux_from_fasta(str(data_dir / "jh_seqs.fasta"), str(out), False)
    got = simple.read_csv(str(out), delimiter="\t")
    assert 0 < len(got) <= len(seqs)
    assert all(r["chain_type"] == "JH" and r["j_cdr3_end"] for r in got)


# ---- mixcr ----

def test_make_mixcr_ref_v_d_j():
    ref = germline_set([v_allele(), d_allele(), j_allele()])
    (rec,) = af.make_mixcr_ref(ref, "uri", 9606, "Homo sapiens")
    assert rec["taxonId"] == 9606 and rec["speciesNames"] == ["Homo sapiens"]
    genes = {g["name"]: g for g in rec["genes"]}
    assert genes["IGHV1-1*01"]["geneType"] == "V"
    assert genes["IGHV1-1*01"]["anchorPoints"]["CDR1Begin"] == 75
    assert genes["IGHV1-1*01"]["anchorPoints"]["CDR3Begin"] == 285
    assert genes["IGHD1-1*01"]["anchorPoints"] == {"DBegin": 0, "DEnd": 17}
    assert genes["IGHJ4*01"]["anchorPoints"] == {"JBegin": 0, "FR4Begin": 21, "FR4End": len(J4)}
    assert genes["IGHJ4*01"]["chains"] == ["IGH"]
    assert genes["IGHJ4*01"]["baseSequence"] == "uri#IGHJ4*01"
    assert len(rec["sequenceFragments"]) == 3
    assert all(g["isFunctional"] for g in rec["genes"])


@pytest.mark.parametrize("flag", [False, "F", "false"])
def test_make_mixcr_ref_non_functional(flag):
    (rec,) = af.make_mixcr_ref(germline_set([j_allele(functional=flag)]), "uri", 1, "x")
    assert rec["genes"][0]["isFunctional"] is False


def test_make_mixcr_ref_v_without_delineations_skipped(capsys):
    a = v_allele()
    del a["v_gene_delineations"]
    (rec,) = af.make_mixcr_ref(germline_set([a]), "uri", 1, "x")
    assert rec["genes"] == []
    assert "no delineations" in capsys.readouterr().out


def test_make_mixcr_ref_bad_allele_skipped(capsys):
    (rec,) = af.make_mixcr_ref(germline_set([j_allele(gene_start=None)]), "uri", 1, "x")
    assert rec["genes"] == []
    assert "Unrecognisable parameter" in capsys.readouterr().out


# ---- leaders ----

def test_leader_from_json_nt_and_aa():
    ref = germline_set([v_allele(), j_allele()])
    nt = af.leader_from_json(ref, "ignored", False)
    assert nt == {"IGHV1-1*01": V_SEQ[0:6] + V_SEQ[6:15]}
    aa = af.leader_from_json(ref, "ignored", True)
    assert aa == {"IGHV1-1*01": simple.translate(V_SEQ[:15])}


def test_leader_from_json_missing_leader(capsys):
    a = v_allele()
    del a["leader_1_start"]
    assert af.leader_from_json(germline_set([a]), "ignored", False) == []
    assert "leader not available" in capsys.readouterr().out


def test_write_leaders_to_fasta(tmp_path):
    out = tmp_path / "l.fasta"
    assert af.write_leaders_to_fasta(germline_set([v_allele()]), str(out), False) is True
    assert simple.read_fasta(str(out)) == {"IGHV1-1*01": V_SEQ[:15]}


def test_write_leaders_to_fasta_none_available(tmp_path):
    out = tmp_path / "l.fasta"
    assert af.write_leaders_to_fasta(germline_set([j_allele()]), str(out), False) is False
    assert not out.exists()
