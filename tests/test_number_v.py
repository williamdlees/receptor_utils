import warnings

import pytest

import receptor_utils
from receptor_utils import number_v
from receptor_utils import simple_bio_seq as simple

# Cases ported from number_v.run_tests(): (query, expected gapped sequence)
NUMBERING_CASES = {
    "1-64": (
        "EMQLVQSEAEVKKPGASVKISCKASGYTFTYRYLHWLRQTPGQGLEWMGWITPYNGNTNYAQKFQDRATITRDRSMSTAYMELSSLRSEDTAVYYCAR",
        "EMQLVQSEA.EVKKPGASVKISCKASGYTF....TYRYLHWLRQTPGQGLEWMGWITPY..NGNTNYAQKFQ.DRATITRDRSMSTAYMELSSLRSEDTAVYYCAR",
    ),
    "2-143": (
        "QVTLKESGPALVKPTQTLTLTCTFSGFSLTTSGMGVGWIRQPPGKALEWLALIYWDDDKRYSTSLKSRLTISKDTSKNQVVLTMTNMDPMDTATYYCARG",
        "QVTLKESGP.ALVKPTQTLTLTCTFSGFSLT..TSGMGVGWIRQPPGKALEWLALIYWD...DDKRYSTSLK.SRLTISKDTSKNQVVLTMTNMDPMDTATYYCARG",
    ),
    "3-8": (
        "EVQLVESGGGLVQPGRSLRPSCAASGFTFSSYGMHWVRQAPEEGLVWVSYIGSSTMYYADSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCVR",
        "EVQLVESGG.GLVQPGRSLRPSCAASGFTF....SSYGMHWVRQAPEEGLVWVSYIGS....STMYYADSVK.GRFTISRDNAKNSLYLQMNSLRAEDTAVYYCVR",
    ),
    "4-86": (
        "QVQLQESGPGLVKPSETLSLTCAVSGGSISSSNWWSWIRQPPGKGLEWIGRISGSGGSTSDNPSLKSRVTISKDTSKNQFSLKLSSVTAADTAVYYCAR",
        "QVQLQESGP.GLVKPSETLSLTCAVSGGSIS...SSNWWSWIRQPPGKGLEWIGRISGS..GGSTSDNPSLK.SRVTISKDTSKNQFSLKLSSVTAADTAVYYCAR",
    ),
    "5-16": (
        "EVQLVQSGAEVKRPGESLKISCKTSGYSFTSYWISWVRQMPGKGLEWMGAIDPSDSDTRYNPSFQGQVTISADKSISTAYLQWSRLKASDTATYYCAK",
        "EVQLVQSGA.EVKRPGESLKISCKTSGYSF....TSYWISWVRQMPGKGLEWMGAIDPS..DSDTRYNPSFQ.GQVTISADKSISTAYLQWSRLKASDTATYYCAK",
    ),
    "6-1": (
        "QVQLQESGPGLVKPSQTLSLTCAISGDSVSSNSATWNWIRQSPSRGLEWLGRTYYRSKWYNDYAQSVQNRISINPDTSKNQFSLQLNSVTPEDMAVYYCAR",
        "QVQLQESGP.GLVKPSQTLSLTCAISGDSVS..SNSATWNWIRQSPSRGLEWLGRTYYRS.KWYNDYAQSVQ.NRISINPDTSKNQFSLQLNSVTPEDMAVYYCAR",
    ),
    "7-108": (
        "QVQLVQSGAEVKQPGASVKVSCKASGYTFTSYGMNWVRQAHGQRLEWMGWINTDTGNPTYAQGFKERFTFSMDTSISTAYLQISSLKAEDTAVYYCAR",
        "QVQLVQSGA.EVKQPGASVKVSCKASGYTF....TSYGMNWVRQAHGQRLEWMGWINTD..TGNPTYAQGFK.ERFTFSMDTSISTAYLQISSLKAEDTAVYYCAR",
    ),
}


@pytest.mark.parametrize("query,expected", NUMBERING_CASES.values(), ids=NUMBERING_CASES.keys())
def test_number_ighv(query, expected):
    gapped, status = number_v.number_ighv(query)
    assert gapped == expected
    assert status == ""
    assert gapped.replace(".", "") == query


def test_number_ighv_no_first_cysteine():
    gapped, status = number_v.number_ighv("A" * 120)
    assert gapped is None
    assert "cysteine" in status


def test_number_ighv_no_trp():
    seq = "A" * 21 + "C" + "A" * 98
    gapped, status = number_v.number_ighv(seq)
    assert gapped is None
    assert "Trp" in status


def test_number_ighv_no_cdr2():
    seq = "A" * 21 + "C" + "A" * 11 + "W" + "A" * 85
    gapped, status = number_v.number_ighv(seq)
    assert gapped is None
    assert "CDR2" in status


def test_nt_diff():
    assert number_v.nt_diff("ACGT", "ACGT") == 0
    assert number_v.nt_diff("ACGT", "ACCT") == 1
    assert number_v.nt_diff("ACNT", "ACGT") == 0
    assert number_v.nt_diff("ACGT", "AC") == 2


def test_match_score():
    assert number_v.match_score("AKLT", [["K"], ["L"], ["T"]], 3) == [(1, 3)]
    assert number_v.match_score("AAAA", [["K"], ["L"]], 1) == []
    assert number_v.match_score("KL", [["K"], ["L"], ["T"]], 0) == []      # target shorter than pattern


def test_distribute_cdr():
    assert number_v.distribute_cdr("ABCD", 8) == "AB....CD"
    assert number_v.distribute_cdr("ABC", 6) == "AB...C"
    assert number_v.distribute_cdr("ABCD", 4) == "ABCD"


def test_insert_space():
    assert number_v.insert_space("ABCD", 2) == "AB CD"
    assert number_v.insert_space("ABCD", 0) == " ABCD"


def test_pretty_gapped(capsys):
    seq = "A" * 104
    number_v.pretty_gapped(seq)
    out = capsys.readouterr().out
    assert "FR1-IMGT" in out
    assert "A" * 26 + " " + "A" * 12 + " " in out


def test_gap_nt_from_aa():
    assert number_v.gap_nt_from_aa("GCTGCTGCT", "A.A.A") == "GCT...GCT...GCT"
    assert number_v.gap_nt_from_aa("GCTGCTGCTGC", "AA") == "GCTGCTGCTGC"       # trailing codon and partial codon retained
    assert number_v.gap_nt_from_aa("GCTGCT", "-A") == "...GCT" + "GCT"


def test_gap_align():
    # leading dots in ref are filled from seq; interior dots become gaps; trailing gaps trimmed
    assert number_v.gap_align("ACGT", "A.CG.T") == "A.CG.T"
    assert number_v.gap_align("ACG", "A.CGTT") == "A.CG"
    assert number_v.gap_align("ACGTTT", "A.CG") == "A.CGTTT"


def test_gap_align_aa():
    assert number_v.gap_align_aa("ACD", "A.CD") == "A.CD"
    assert number_v.gap_align_aa("AC", "A.CDE") == "A.C"
    assert number_v.gap_align_aa("ACDEF", "A.C") == "A.CDEF"


def test_gap_align_aa_from_nt():
    assert number_v.gap_align_aa_from_nt("ACD", "GCT...TGTGAT") == "A.CD"
    assert number_v.gap_align_aa_from_nt("A", "GCT...TGT") == "A."      # a gap is emitted before the residues run out
    assert number_v.gap_align_aa_from_nt("ACD", "...GCTTGTGAT") == "ACD"       # leading gaps consume residues


def test_check_conserved_residues_ok():
    aa = "A" * 22 + "C" + "A" * 17 + "W" + "A" * 62 + "C" + "A" * 5
    assert number_v.check_conserved_residues(aa) == ""


def test_check_conserved_residues_problems():
    aa = "A" * 110
    notes = number_v.check_conserved_residues(aa)
    assert "First cysteine not found" in notes
    assert "Conserved Trp not found" in notes
    assert "Second cysteine not found" in notes


def test_check_conserved_residues_stop_and_truncation():
    assert number_v.check_conserved_residues("A" * 50) == "Sequence truncated before second cysteine"
    aa = "A" * 22 + "C" + "A" * 17 + "W" + "A" * 10 + "*" + "A" * 51 + "C" + "A" * 5
    assert "Stop codon" in number_v.check_conserved_residues(aa)


def test_gap_sequence_with_gapped_reference(tmp_path):
    ungapped = {"IGHV1-1*01": "GCTGCTGCTGCTTGT", "IGHV1-1*02": "GCTGCTGCTGCTAGT"}
    gapped = {"IGHV1-1*01": "GCTGCT...GCTGCTTGT", "IGHV1-1*02": "GCTGCT...GCTGCTAGT"}
    res, aa, notes = number_v.gap_sequence("GCTGCTGCTGCTAGT", gapped, ungapped)
    assert res == "GCTGCT...GCTGCTAGT"
    assert aa == "AA.AAS"
    assert notes == "Sequence truncated before second cysteine"


def test_gap_sequence_non_v_has_no_notes():
    ungapped = {"IGHJ1*01": "TGGGGC"}
    gapped = {"IGHJ1*01": "TGGGGC"}
    res, aa, notes = number_v.gap_sequence("TGGGGC", gapped, ungapped)
    assert (res, aa, notes) == ("TGGGGC", "WG", "")


def test_gap_sequence_with_real_reference(data_dir):
    gapped = simple.read_fasta(str(data_dir / "V_gapped.fasta"))
    ungapped = {k: v.replace(".", "") for k, v in gapped.items()}
    name = "IGHV1-18*01"
    res, aa, notes = number_v.gap_sequence(ungapped[name], gapped, ungapped)
    assert res == gapped[name]
    assert notes == ""
    assert "*" not in aa


def test_package_number_ighv_attribute_is_deprecated():
    with pytest.warns(DeprecationWarning, match="number_v"):
        mod = receptor_utils.number_ighv
    assert mod is number_v


def test_package_unknown_attribute_raises():
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with pytest.raises(AttributeError):
            receptor_utils.does_not_exist
