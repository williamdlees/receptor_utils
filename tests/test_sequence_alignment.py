import pytest

from receptor_utils import sequence_alignment as sa

REF = "GCTAAACCCGGGTTTAAA"      # A K P G F K
SEQS = {
    "A*01": REF,
    "A*02": "GCTAAGCCCGGGTTTAAA",     # silent change in codon 2
    "A*03": "GCTAAACCCGGGTCTCAA",     # codon 5 F->S (non-silent), codon 6 K->Q (non-silent)
}


def lines(text):
    return text.split("\n")


def test_empty_input_returns_empty_string():
    assert sa.create_alignment({}, sequence_type="J") == ""


def test_invalid_sequence_type():
    with pytest.raises(ValueError):
        sa.create_alignment(SEQS, sequence_type="X")


def test_unsupported_output_format():
    with pytest.raises(NotImplementedError):
        sa.create_alignment(SEQS, output_format="html")


@pytest.mark.parametrize("coords", [(1, 2, 3), (1, 2, 3, 4, 5, 6), (5, 4, 6, 7, 8), (1, 2, 2, 3, 4)])
def test_invalid_v_coords(coords):
    with pytest.raises(ValueError):
        sa.create_alignment(SEQS, sequence_type="V", v_coords=coords)


def test_j_alignment_shows_reference_and_differences():
    out = sa.create_alignment(SEQS, sequence_type="J")
    assert "GCT AAA CCC GGG TTT AAA" in out           # reference in codons
    assert "A   K   P   G   F   K" in out             # reference translation
    assert "A*02" in out and "A*03" in out
    assert "---" in out
    assert "--g" in out                               # silent change shown as nucleotide only
    assert "-c-" in out and "c--" in out              # codons 5 and 6 differ
    assert "S" in out and "Q" in out                  # non-silent changes show the new amino acid


def test_j_alignment_has_no_cdr_labels():
    out = sa.create_alignment(SEQS, sequence_type="J")
    assert "CDR" not in out


def test_d_alignment_has_no_amino_acids():
    out = sa.create_alignment(SEQS, sequence_type="D")
    body = [ln for ln in lines(out) if ln.startswith("A*")]
    assert len(body) == 3
    assert "A   K" not in out


def test_sequences_sorted_by_name_first_is_reference():
    out = sa.create_alignment({"B": "GCTGCT", "A": "GCTGCA"}, sequence_type="J")
    names = [ln.split()[0] for ln in lines(out) if ln[:1] in "AB" and ln.strip()]
    assert names == ["A", "B"]
    assert "GCT GCA" in out


def test_v_alignment_with_default_coords_labels_cdrs():
    seq = "GCT" * 110
    out = sa.create_alignment({"V1*01": seq, "V1*02": seq[:-3] + "GCC"}, sequence_type="V", codon_wrap=120)
    assert "CDR1" in out and "CDR2" in out and "CDR3" in out


def test_v_alignment_with_custom_coords_and_wrap():
    seq = "GCT" * 30
    out = sa.create_alignment({"V1*01": seq, "V1*02": seq}, sequence_type="V", codon_wrap=10, v_coords=(2, 4, 12, 14, 25))
    assert "CDR1" in out
    assert out.count("1  ") >= 1
    # three wrapped chunks, so three reference lines
    assert sum(1 for ln in lines(out) if ln.startswith("V1*01")) == 3


def test_codon_wrap_chunks_header_numbers():
    seq = "GCT" * 25
    out = sa.create_alignment({"X*01": seq}, sequence_type="J", codon_wrap=10)
    assert sum(1 for ln in lines(out) if ln.startswith("X*01")) == 3
    assert "21" in out and "25" in out


def test_partial_trailing_codon_and_shorter_comparison():
    seqs = {"A": "GCTAAACC", "B": "GCTAA"}
    out = sa.create_alignment(seqs, sequence_type="J")
    assert "GCT AAA CC" in out
    assert "--- --" in out


def test_gaps_are_preserved_in_comparison():
    seqs = {"A": "GCT...AAA", "B": "GCT...AAA"}
    out = sa.create_alignment(seqs, sequence_type="J")
    assert "... " in out


def test_comparison_longer_than_reference():
    out = sa.create_alignment({"A": "GCT", "B": "GCTAAA"}, sequence_type="J")
    assert "aaa" in out


def test_chunk_header_no_cdr_numbering():
    hdr = sa._create_chunk_header_no_cdr(1, 3, 10)
    assert hdr[:2] == ["", ""]
    assert hdr[2].strip().split() == ["1", "2", "3"]


def test_chunk_header_v_regions_complete_and_partial():
    coords = (2, 4, 6, 8, 10)
    whole = sa._create_chunk_header(1, 12, coords, 5)
    assert "CDR1" in whole[0] and "CDR2" in whole[0] and "CDR3" in whole[0]
    start_only = sa._create_chunk_header(1, 3, coords, 5)
    assert start_only[0].rstrip().endswith("CDR1")
    end_only = sa._create_chunk_header(3, 5, (1, 4, 20, 21, 30), 5)
    assert "CDR1" in end_only[0].strip()
    inside = sa._create_chunk_header(3, 4, (1, 10, 20, 21, 30), 5)
    assert set(inside[0].strip()) == {"_"}


def test_chunk_header_triple_digit_numbers():
    hdr = sa._create_chunk_header(99, 101, (2, 4, 6, 8, 200), 5)
    assert "101" in hdr[-1]


def test_reference_and_comparison_lines_directly():
    ref = sa._create_reference_line("R", "GCTAAA", "AK", 0, 2, 6, 1, "J")
    aa_line, nt_line = ref.split("\n")
    assert aa_line.strip().split() == ["A", "K"]
    assert nt_line.startswith("R") and "GCT AAA" in nt_line
    d_ref = sa._create_reference_line("R", "GCTAAA", "AK", 0, 2, 6, 1, "D")
    assert "\n" not in d_ref
    comp = sa._create_comparison_line("Q", "GCTCAA", "GCTAAA", "AK", 0, 2, 6, 1, "J")
    assert "Q" in comp.split("\n")[1] and "c--" in comp
    same = sa._create_comparison_line("Q", "GCTAAA", "GCTAAA", "AK", 0, 2, 6, 1, "J")
    assert same.split("\n")[1].rstrip().endswith("--- ---")
