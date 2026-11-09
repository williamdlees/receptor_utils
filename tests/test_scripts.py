import json
import runpy
import sys

import pytest

from receptor_utils import simple_bio_seq as simple

from test_aux_formats import J4, germline_set, j_allele, v_allele, d_allele


def out_of(capsys):
    return capsys.readouterr().out


# ---- name_allele ----

def test_name_allele_exact_match(run_script, capsys, data_dir):
    ref = simple.read_fasta(str(data_dir / "V_gapped.fasta"))
    seq = ref["IGHV1-18*01"].replace(".", "")
    assert run_script("name_allele", [data_dir / "V_gapped.fasta", seq]) == 0
    assert out_of(capsys).splitlines()[0] == "IGHV1-18*01"


def test_name_allele_novel_snp_and_gene_filter(run_script, capsys, data_dir):
    ref = simple.read_fasta(str(data_dir / "V_gapped.fasta"))
    seq = ref["IGHV1-18*01"].replace(".", "")
    mut = seq[:150] + ("A" if seq[150] != "A" else "C") + seq[151:]
    assert run_script("name_allele", [data_dir / "V_gapped.fasta", mut, "-g", "IGHV1-18"]) == 0
    out = out_of(capsys)
    assert out.splitlines()[0].startswith("IGHV1-18*01_")
    assert "Gapped sequence:" in out


def test_name_allele_reverse_complement(run_script, capsys, data_dir):
    ref = simple.read_fasta(str(data_dir / "V_gapped.fasta"))
    seq = ref["IGHV1-18*01"].replace(".", "")
    assert run_script("name_allele", [data_dir / "V_gapped.fasta", simple.reverse_complement(seq), "-r"]) == 0
    assert out_of(capsys).splitlines()[0] == "IGHV1-18*01"


def test_name_allele_d_gene(run_script, capsys, data_dir):
    d = simple.read_fasta(str(data_dir / "D.fasta"))
    assert run_script("name_allele", [data_dir / "D.fasta", d["IGHD1-1*01"]]) == 0
    assert out_of(capsys).splitlines()[0] == "IGHD1-1*01"


def test_name_allele_unknown_gene_exits(run_script, capsys, data_dir):
    assert run_script("name_allele", [data_dir / "D.fasta", "ACGT", "-g", "NOPE"]) == 0
    assert "No reference genes" in out_of(capsys)


def test_name_allele_mixed_loci_exits(run_script, capsys, write_fasta_file):
    ref = write_fasta_file("mixed.fasta", {"IGHJ1*01": "ACGT", "IGKJ1*01": "ACGT"})
    assert run_script("name_allele", [ref, "ACGT"]) == 0
    assert "more than one locus" in out_of(capsys)


def test_name_allele_empty_reference_exits(run_script, capsys, tmp_path):
    empty = tmp_path / "empty.fasta"
    empty.write_text("")
    assert run_script("name_allele", [empty, "ACGT"]) == 0
    assert "No valid sequences" in out_of(capsys)


def test_name_allele_missing_args(run_script):
    assert run_script("name_allele", []) == 2


# ---- gap_sequences (gap_inferred) ----

def test_gap_sequences(run_script, tmp_path, data_dir, capsys):
    ref = simple.read_fasta(str(data_dir / "V_gapped.fasta"))
    ungapped = ref["IGHV1-18*01"].replace(".", "")
    inferred = tmp_path / "inferred.fasta"
    simple.write_fasta(str(inferred), {"new1": ungapped})
    out = tmp_path / "out.fasta"
    assert run_script("gap_inferred", [inferred, data_dir / "V_gapped.fasta", out]) == 0
    assert simple.read_fasta(str(out)) == {"new1": ref["IGHV1-18*01"]}


def test_gap_sequences_removes_bad_reference(run_script, tmp_path, data_dir, capsys):
    ref = simple.read_fasta(str(data_dir / "V_gapped.fasta"))
    good = ref["IGHV1-18*01"]
    refs = tmp_path / "refs.fasta"
    simple.write_fasta(str(refs), {"IGHV1-18*01": good, "BAD*01": "A" * len(good)})
    inferred = tmp_path / "inferred.fasta"
    simple.write_fasta(str(inferred), {"s": good.replace(".", "")})
    out = tmp_path / "out.fasta"
    assert run_script("gap_inferred", [inferred, refs, out]) == 0
    assert "Removing BAD*01 from reference" in out_of(capsys)
    assert simple.read_fasta(str(out))["s"] == good


def test_gap_sequences_missing_args(run_script):
    assert run_script("gap_inferred", ["only_one"]) == 2


# ---- make_igblast_ndm ----

def test_make_igblast_ndm_from_fasta(run_script, tmp_path, data_dir):
    out = tmp_path / "x.ndm"
    assert run_script("make_igblast_ndm", [data_dir / "V_gapped.fasta", "VH", out]) == 0
    lines = out.read_text().splitlines()
    assert lines[0].startswith("#name")
    assert len(lines) == 1 + len(simple.read_fasta(str(data_dir / "V_gapped.fasta")))
    assert lines[1].split("\t")[-2] == "VH"


def test_make_igblast_ndm_from_fasta_custom_coords(run_script, tmp_path, data_dir):
    out = tmp_path / "x.ndm"
    assert run_script("make_igblast_ndm", [data_dir / "V_gapped.fasta", "VH", out, "-c", "79,114,166,195,313"]) == 0
    assert out.exists()


def test_make_igblast_ndm_from_json(run_script, tmp_path, data_dir):
    out = tmp_path / "x.ndm"
    assert run_script("make_igblast_ndm", [data_dir / "Human_IGH_VDJ_rev_2.json", "VH", out]) == 0
    assert out.read_text().count("\n") > 100


@pytest.mark.parametrize("extra", [
    ["-c", "1,2,3"],                 # wrong count
    ["-c", "1,2,3,4,x"],             # non-integer
    ["-c", "5,4,6,7,8"],             # not ascending
])
def test_make_igblast_ndm_bad_coords(run_script, tmp_path, data_dir, extra, capsys):
    assert run_script("make_igblast_ndm", [data_dir / "V_gapped.fasta", "VH", tmp_path / "x.ndm"] + extra) == 1
    assert "Error" in out_of(capsys)


def test_make_igblast_ndm_coords_not_allowed_for_json(run_script, tmp_path, data_dir):
    assert run_script("make_igblast_ndm", [data_dir / "Human_IGH_VDJ_rev_2.json", "VH", tmp_path / "x.ndm", "-c", "1,2,3,4,5"]) == 1


def test_make_igblast_ndm_bad_chain(run_script, tmp_path, data_dir, capsys):
    assert run_script("make_igblast_ndm", [data_dir / "V_gapped.fasta", "XX", tmp_path / "x.ndm"]) == 1
    assert "chain must be one of" in out_of(capsys)


def test_make_igblast_ndm_bad_extension(run_script, tmp_path):
    f = tmp_path / "ref.txt"
    f.write_text("x")
    assert run_script("make_igblast_ndm", [f, "VH", tmp_path / "x.ndm"]) == 1


# ---- annotate_j ----

def test_annotate_j_fasta(run_script, tmp_path, write_fasta_file):
    fasta = write_fasta_file("j.fasta", {"IGHJ4*01": J4})
    out = tmp_path / "x.aux"
    assert run_script("annotate_j", [fasta, out]) == 0
    rows = simple.read_csv(str(out), delimiter="\t")
    assert rows[0]["#name"] == "IGHJ4*01"


def test_annotate_j_verbose(run_script, tmp_path, write_fasta_file, capsys):
    fasta = write_fasta_file("j.fasta", {"IGHJ4*01": J4})
    assert run_script("annotate_j", [fasta, tmp_path / "x.aux", "-v"]) == 0
    assert "MiAIRR" in out_of(capsys)


def test_annotate_j_json(run_script, tmp_path):
    src = tmp_path / "set.json"
    src.write_text(json.dumps(germline_set([j_allele()])))
    out = tmp_path / "x.aux"
    assert run_script("annotate_j", [src, out]) == 0
    assert simple.read_csv(str(out), delimiter="\t")[0]["#name"] == "IGHJ4*01"


def test_annotate_j_bad_extension(run_script, tmp_path, capsys):
    f = tmp_path / "j.txt"
    f.write_text("x")
    assert run_script("annotate_j", [f, tmp_path / "x.aux"]) == 1
    assert "must be a .json or .fasta" in out_of(capsys)


# ---- extract_refs ----

IMGT_TEXT = (
    ">X1|IGHV1-1*01|Homo sapiens|F|V-REGION|1..10|10 nt|\n" + "A.A" * 100 + "\n"
    ">X2|IGHD1-1*01|Homo sapiens|F|D-REGION|1..10|10 nt|\nGG.TT\n"
    ">X3|IGHJ1*01|Homo sapiens|F|J-REGION|1..10|10 nt|\nTT.GG\n"
    ">X4|IGHG1*01|Homo sapiens|F|CH1|1..10|10 nt|\nCCCC\n"
)


def test_extract_refs_from_local_file(run_script, tmp_path, monkeypatch):
    ref = tmp_path / "imgt.fasta"
    ref.write_text(IMGT_TEXT)
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)
    assert run_script("extract_refs", ["Homo sapiens", "-r", ref]) == 0
    names = sorted(p.name for p in work.iterdir())
    assert "Homo_sapiens_IGHV_gapped.fasta" in names
    assert "Homo_sapiens_IGHV.fasta" in names
    assert "Homo_sapiens_IGHD.fasta" in names
    assert "Homo_sapiens_IGHJ.fasta" in names
    assert "Homo_sapiens_CH.fasta" in names
    assert simple.read_fasta(str(work / "Homo_sapiens_IGHD.fasta")) == {"IGHD1-1*01": "GGTT"}
    assert "." in simple.read_fasta(str(work / "Homo_sapiens_IGHV_gapped.fasta"))["IGHV1-1*01"]


def test_extract_refs_functional_only_and_locus(run_script, tmp_path, monkeypatch):
    ref = tmp_path / "imgt.fasta"
    ref.write_text(
        ">X1|IGKV1-1*01|Homo sapiens|F|V-REGION|1..10|10 nt|\n" + "ACG" * 100 + "\n"
        ">X2|IGKJ1*01|Homo sapiens|F|J-REGION|1..10|10 nt|\nTTGG\n"
    )
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)
    assert run_script("extract_refs", ["Homo sapiens", "-r", ref, "-L", "IGK", "-F"]) == 0
    assert (work / "Homo_sapiens_IGKV_gapped.fasta").exists()
    assert (work / "Homo_sapiens_IGKJ.fasta").exists()


def test_extract_refs_downloads_when_no_ref_given(run_script, tmp_path, monkeypatch):
    import io
    import scripts.extract_refs as er

    class FakeResp(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    calls = []

    def fake_urlopen(url):
        calls.append(url)
        return FakeResp(IMGT_TEXT.encode())

    monkeypatch.setattr(er.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.chdir(tmp_path)
    assert run_script("extract_refs", ["Homo sapiens"]) == 0
    assert len(calls) == 1 and calls[0].startswith("https://www.imgt.org/")
    assert (tmp_path / "Homo_sapiens_IGHV_gapped.fasta").exists()


# ---- fix_macaque_gaps ----

def test_fix_macaque_gaps_removes_blank_codons(run_script, tmp_path, write_fasta_file):
    seq = "".join(f"{i:03d}".replace("0", "A").replace("1", "C").replace("2", "G").replace("3", "T") for i in range(1, 31))
    seq_list = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    seq_list[15] = "..."          # codon 16
    seq_list[27] = "..."          # codon 28
    gapped = "".join(seq_list)
    bad = "".join(seq_list[:15] + ["AAA"] + seq_list[16:])      # codon 16 not blank
    infile = write_fasta_file("in.fasta", {"ok": gapped, "bad": bad})
    out = tmp_path / "out.fasta"
    assert run_script("fix_macaque_gaps", [infile, out, "IGH"]) == 0
    res = simple.read_fasta(str(out))
    assert list(res) == ["ok"]
    assert len(res["ok"]) == len(gapped) - 6
    assert "." not in res["ok"]


def test_fix_macaque_gaps_invalid_chain(run_script, tmp_path, write_fasta_file, capsys):
    infile = write_fasta_file("in.fasta", {"a": "ACGT"})
    assert run_script("fix_macaque_gaps", [infile, tmp_path / "o.fasta", "TRB"]) == 0
    assert "'chain' must be one of" in out_of(capsys)
    assert not (tmp_path / "o.fasta").exists()


# ---- identical_seqs ----

def test_identical_seqs(run_script, write_fasta_file, capsys):
    f = write_fasta_file("x.fasta", {"a": "ACGTACGT", "b": "ACGTACGT", "c": "ACGT", "d": "TTTT"})
    assert run_script("identical_seqs", [f]) == 0
    out = out_of(capsys)
    assert "Identical sequences:\na, b" in out
    assert "c (ACGT)is a sub-sequence of a (ACGTACGT)" in out
    assert "d" not in out.split("Sub-sequences:")[1].replace("(", "").split("sub-sequence")[0].split()[:1]


def test_identical_seqs_no_duplicates(run_script, write_fasta_file, capsys):
    f = write_fasta_file("x.fasta", {"a": "AAAA", "b": "CCCC"})
    assert run_script("identical_seqs", [f]) == 0
    assert out_of(capsys) == "Identical sequences:\nSub-sequences:\n"


# ---- rev_comp ----

def test_rev_comp(run_script, capsys):
    assert run_script("reverse_complement", ["AACGT"]) == 0
    assert out_of(capsys).strip() == "ACGTT"


def test_rev_comp_missing_arg(run_script):
    assert run_script("reverse_complement", []) == 2


# ---- at_coords ----

@pytest.fixture
def coords_fasta(write_fasta_file):
    return write_fasta_file("c.fasta", {"s1": "AACCGGTT", "s2": "TTTTAAAA"})


def test_at_coords_named(run_script, coords_fasta, capsys):
    assert run_script("at_coords", [coords_fasta, 3, 6, "-n", "s1"]) == 0
    assert out_of(capsys).strip() == "CCGG"


def test_at_coords_single_sequence_default(run_script, write_fasta_file, capsys):
    f = write_fasta_file("one.fasta", {"s": "AACCGGTT"})
    assert run_script("at_coords", [f, 1, 2]) == 0
    assert out_of(capsys).strip() == "AA"


def test_at_coords_reverse_complement_result(run_script, coords_fasta, capsys):
    assert run_script("at_coords", [coords_fasta, 1, 4, "-n", "s1", "-r"]) == 0
    assert out_of(capsys).strip() == "GGTT"


def test_at_coords_extract_from_reverse_complement(run_script, coords_fasta, capsys):
    assert run_script("at_coords", [coords_fasta, 1, 3, "-n", "s1", "-xr"]) == 0
    assert out_of(capsys).strip() == "AAC"       # revcomp of AACCGGTT is AACCGGTT


def test_at_coords_unknown_name(run_script, coords_fasta, capsys):
    assert run_script("at_coords", [coords_fasta, 1, 3, "-n", "nope"]) == 0
    assert "not found in file" in out_of(capsys)


@pytest.mark.parametrize("start,end", [(5, 2), (0, 3), (-1, 3)])
def test_at_coords_invalid_range_prints_nothing(run_script, coords_fasta, capsys, start, end):
    assert run_script("at_coords", [coords_fasta, start, end, "-n", "s1"]) == 0
    assert out_of(capsys) == ""


# ---- merge_fasta ----

def test_merge_fasta(run_script, tmp_path, write_fasta_file, capsys):
    f1 = write_fasta_file("1.fasta", {"a": "AAAA", "b": "CCCC"})
    f2 = write_fasta_file("2.fasta", {"c": "AAAA", "d": "GGGG"})
    out = tmp_path / "out.fasta"
    assert run_script("merge_fasta", [f1, f2, out]) == 0
    assert simple.read_fasta(str(out)) == {"a": "AAAA", "b": "CCCC", "d": "GGGG"}
    assert "1 items merged" in out_of(capsys)


# ---- create_alignment ----

@pytest.fixture
def alignment_fasta(write_fasta_file):
    seq = "GCT" * 110
    return write_fasta_file("al.fasta", {"IGHV1-1*01": seq, "IGHV1-1*02": seq[:-3] + "GCC", "IGHV2-1*01": seq})


def test_create_alignment_v(run_script, alignment_fasta, tmp_path, capsys):
    out = tmp_path / "al.txt"
    assert run_script("create_alignment", [alignment_fasta, "V", out]) == 0
    assert "Alignment written to" in out_of(capsys)
    text = out.read_text()
    assert "CDR1" in text and "IGHV1-1*02" in text


def test_create_alignment_filter_wrap_and_coords(run_script, alignment_fasta, tmp_path):
    out = tmp_path / "al.txt"
    assert run_script("create_alignment", [alignment_fasta, "V", out, "-f", "IGHV1-1", "-w", "30", "-c", "2,4,6,8,10"]) == 0
    text = out.read_text()
    assert "IGHV2-1*01" not in text


def test_create_alignment_d_and_j(run_script, alignment_fasta, tmp_path):
    for t in ("D", "J"):
        out = tmp_path / f"{t}.txt"
        assert run_script("create_alignment", [alignment_fasta, t, out]) == 0
        assert "CDR" not in out.read_text()


@pytest.mark.parametrize("coords", ["1,2,3", "1,2,3,4,x", "5,4,6,7,8"])
def test_create_alignment_bad_coords(run_script, alignment_fasta, tmp_path, coords):
    assert run_script("create_alignment", [alignment_fasta, "V", tmp_path / "o.txt", "-c", coords]) == 1


def test_create_alignment_filter_matches_nothing(run_script, alignment_fasta, tmp_path, capsys):
    assert run_script("create_alignment", [alignment_fasta, "V", tmp_path / "o.txt", "-f", "ZZZ"]) == 1
    assert "No sequences matching" in out_of(capsys)


def test_create_alignment_empty_input(run_script, tmp_path, capsys):
    empty = tmp_path / "e.fasta"
    empty.write_text("")
    assert run_script("create_alignment", [empty, "V", tmp_path / "o.txt"]) == 1
    assert "No sequences found" in out_of(capsys)


def test_create_alignment_missing_input(run_script, tmp_path, capsys):
    assert run_script("create_alignment", [tmp_path / "nope.fasta", "V", tmp_path / "o.txt"]) == 1
    assert "Error reading input file" in out_of(capsys)


def test_create_alignment_unwritable_output(run_script, alignment_fasta, tmp_path, capsys):
    assert run_script("create_alignment", [alignment_fasta, "V", tmp_path / "no_dir" / "o.txt"]) == 1
    assert "Error writing output file" in out_of(capsys)


def test_create_alignment_invalid_sequence_type(run_script, alignment_fasta, tmp_path):
    assert run_script("create_alignment", [alignment_fasta, "X", tmp_path / "o.txt"]) == 2


# ---- make_mixcr_json ----

def write_set(tmp_path, alleles=None, **gs_extra):
    gs = {"locus": "IGH", "allele_descriptions": alleles or [v_allele(), d_allele(), j_allele()]}
    gs.update(gs_extra)
    p = tmp_path / "set.json"
    p.write_text(json.dumps({"GermlineSet": [gs]}))
    return p


def test_make_mixcr_json_uses_species_from_file(run_script, tmp_path):
    src = write_set(tmp_path, species={"id": "NCBITAXON:9606", "label": "Homo sapiens"})
    out = tmp_path / "mixcr.json"
    assert run_script("make_mixcr_json", [src, out]) == 0
    data = json.loads(out.read_text())
    assert data[0]["taxonId"] == "9606"
    assert data[0]["speciesNames"] == ["Homo sapiens"]
    assert len(data[0]["genes"]) == 3


def test_make_mixcr_json_command_line_overrides(run_script, tmp_path):
    src = write_set(tmp_path)
    out = tmp_path / "mixcr.json"
    assert run_script("make_mixcr_json", [src, out, "-t", "10090", "-s", "Mus musculus"]) == 0
    data = json.loads(out.read_text())
    assert data[0]["taxonId"] == "10090" and data[0]["speciesNames"] == ["Mus musculus"]


def test_make_mixcr_json_missing_taxon(run_script, tmp_path, capsys):
    src = write_set(tmp_path)
    out = tmp_path / "mixcr.json"
    assert run_script("make_mixcr_json", [src, out]) == 0
    assert "Taxon ID not found" in out_of(capsys)
    assert not out.exists()


def test_make_mixcr_json_missing_species_name(run_script, tmp_path, capsys):
    src = write_set(tmp_path)
    out = tmp_path / "mixcr.json"
    assert run_script("make_mixcr_json", [src, out, "-t", "1"]) == 0
    assert "Species name not found" in out_of(capsys)
    assert not out.exists()


def test_make_mixcr_json_dict_germline_set(run_script, tmp_path):
    p = tmp_path / "set.json"
    p.write_text(json.dumps({"GermlineSet": {"locus": "IGH", "allele_descriptions": [j_allele()],
                                             "species": {"id": "1", "label": "x"}}}))
    out = tmp_path / "mixcr.json"
    assert run_script("make_mixcr_json", [p, out]) == 0
    assert json.loads(out.read_text())[0]["taxonId"] == "1"


# ---- make_leader_fasta (script body runs under __main__ and defines no main()) ----

def run_leader_script(monkeypatch, args):
    monkeypatch.setattr(sys, "argv", ["make_leader_fasta"] + [str(a) for a in args])
    runpy.run_module("scripts.make_leader_fasta", run_name="__main__")


def test_make_leader_fasta(monkeypatch, tmp_path, capsys):
    src = write_set(tmp_path, [v_allele()])
    out = tmp_path / "l.fasta"
    run_leader_script(monkeypatch, [src, out])
    assert "successfully created" in out_of(capsys)
    assert len(simple.read_fasta(str(out))["IGHV1-1*01"]) == 15


def test_make_leader_fasta_aa(monkeypatch, tmp_path):
    src = write_set(tmp_path, [v_allele()])
    out = tmp_path / "l.fasta"
    run_leader_script(monkeypatch, [src, out, "--aa"])
    assert len(simple.read_fasta(str(out))["IGHV1-1*01"]) == 5


def test_make_leader_fasta_none_found(monkeypatch, tmp_path, capsys):
    src = write_set(tmp_path, [j_allele()])
    out = tmp_path / "l.fasta"
    run_leader_script(monkeypatch, [src, out])
    assert "No leader sequences were found" in out_of(capsys)
    assert not out.exists()
