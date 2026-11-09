import pytest

from receptor_utils import novel_allele_name as nan
from receptor_utils import simple_bio_seq as simple


@pytest.fixture(scope="module")
def v_gapped(data_dir):
    return simple.read_fasta(str(data_dir / "V_gapped.fasta"))


@pytest.fixture(scope="module")
def d_ref(data_dir):
    return simple.read_fasta(str(data_dir / "D.fasta"))


@pytest.fixture(scope="module")
def j_ref(data_dir):
    return simple.read_fasta(str(data_dir / "J.fasta"))


@pytest.fixture(scope="module")
def v_novels(data_dir):
    return simple.read_fasta(str(data_dir / "IGHV3-20.04.fasta"))


# Ported from novel_allele_name.run_tests(): each novel in IGHV3-20.04.fasta must be named by its FASTA id.
# The IGHV4-59 entry is 5'-truncated (leading '.'), which is covered separately below.
V_NOVEL_NAMES = [
    "IGHV1-18*01_g107a_g119a_t167a_c173g_g210a_321aca323",
    "IGHV2-70*15_c1-_c34g_i35ag_g36t",
    "IGHV1-69*01_a162g_t184c_c191t",
    "IGHV3-62*03_a26g_c39a_t116a_c137t_a139g_t151g_a188g_c291g_t296c_g297a_t306c_t310-",
    "IGHV3-20*04_t279c_c300t_t301g_321tgactgac328",
    "IGHV2-70*15_c34g_i35ag_g36t",
    "IGHV2-70*15_c34g_i35ag_g36t_c322-",
    "IGHV3-20*04_t279c_c300t_t301g",
    "IGHV1-18*01_g107a_g119a_t167a_c173g_g210a_g254a",
    "IGHV1-18*01_c66-_t67-_g68-_g107a_g119a_t167a_c173g_i198aa_g210a_g254a_i273a",
    "IGHV1-18*01_c1-_a2-_i58cc_g107a_g119a_t167a_c173g_g210a_g254a_321aca323",
    "IGHV1-2*02_i7g_c44t_g128n",
    "IGHV1-2*02_i19g_g42-_c44t_g128n",
    "IGHV1-2*02_c44t_g128n",
    "IGHV1-2*02_c44t_g128n_i133c",
    "IGHV1-2*02_c44t_g128n_t247-",
]


@pytest.mark.parametrize("expected", V_NOVEL_NAMES)
def test_name_novel_v_matches_expected_name(expected, v_novels, v_gapped):
    seq = v_novels[expected]
    name, gapped, notes = nan.name_novel(seq, v_gapped)
    assert name == expected
    assert isinstance(notes, str)


def test_name_novel_v_exact_match_returns_reference_name(v_gapped):
    ref_name = "IGHV1-18*01"
    name, gapped, notes = nan.name_novel(v_gapped[ref_name].replace(".", ""), v_gapped)
    assert name == ref_name
    assert gapped == v_gapped[ref_name]


def test_name_novel_v_accepts_gapped_input(v_gapped):
    ref_name = "IGHV1-18*01"
    seq = v_gapped[ref_name]
    mutated = seq[:100] + ("A" if seq[100] != "A" else "C") + seq[101:]
    name, gapped, _ = nan.name_novel(mutated, v_gapped)
    assert name.startswith("IGHV1-18*")
    assert name != "IGHV1-18*01"


def test_name_novel_is_case_insensitive(v_gapped):
    seq = v_gapped["IGHV1-18*01"].replace(".", "")
    assert nan.name_novel(seq.lower(), v_gapped)[0] == "IGHV1-18*01"


def test_name_novel_five_prime_truncated_uses_ambiguous_ref(v_novels, v_gapped):
    key = [k for k in v_novels if k.startswith("IGHV4-59*01_11_")][0]
    name, gapped, _ = nan.name_novel(v_novels[key], v_gapped)
    assert name.startswith("IGHV4-59*01")
    assert "g226t_g267t_321gaa323" in name


def test_name_novel_snp_in_d(d_ref):
    base = d_ref["IGHD1-1*01"]
    mutated = base[:3] + ("A" if base[3] != "A" else "C") + base[4:]
    name, seq, notes = nan.name_novel(mutated, d_ref, v_gene=False)
    assert name == f"IGHD1-1*01_{base[3].lower()}4{mutated[3].lower()}"
    assert seq == mutated
    assert notes == ""


def test_name_novel_exact_d_and_j(d_ref, j_ref):
    assert nan.name_novel(d_ref["IGHD1-1*01"], d_ref, v_gene=False)[0] == "IGHD1-1*01"
    assert nan.name_novel(j_ref["IGHJ1*01"], j_ref, v_gene=False)[0] == "IGHJ1*01"


def test_name_novel_d_insertion_and_deletion(d_ref):
    base = d_ref["IGHD1-1*01"]
    ins = base[:5] + "GG" + base[5:]
    name, seq, _ = nan.name_novel(ins, d_ref, v_gene=False)
    assert name.startswith("IGHD1-1*01_i")
    assert seq == ins
    dele = base[:5] + base[6:]
    name, _, _ = nan.name_novel(dele, d_ref, v_gene=False)
    assert name == f"IGHD1-1*01_{base[5].lower()}6-"


def test_name_novel_j_multiple_snps(j_ref):
    base = j_ref["IGHJ1*01"]
    chars = list(base)
    chars[2] = "A" if chars[2] != "A" else "C"
    chars[10] = "A" if chars[10] != "A" else "C"
    name, _, _ = nan.name_novel("".join(chars), j_ref, v_gene=False)
    parts = name.split("_")
    assert parts[0] == "IGHJ1*01"
    assert len(parts) == 3


def test_aligned_diff_identical_and_different():
    a, b, score = nan.aligned_diff("ACGTACGT", "ACGTACGT")
    assert a == b == "ACGTACGT"
    _, _, score2 = nan.aligned_diff("ACGTACGT", "ACGTTCGT")
    assert score2 < score


def test_closest_aligned_ref_returns_best_and_ties():
    ref = {"r1": "ACGTACGTAC", "r2": "TTTTTTTTTT", "r3": "ACGTACGTAC"}
    assert nan.closest_aligned_ref("ACGTACGTAC", ref) == ["r1", "r3"]
    assert nan.closest_aligned_ref("TTTTTTTTTT", ref) == ["r2"]


def test_find_indels_insertion_in_novel():
    ref = "ACGTACGTACGTACGT"
    novel = "ACGTACGTAAACGTACGT"
    fixed, insertions = nan.find_indels(ref, novel)
    assert fixed == ref
    assert insertions and "".join(i[1] for i in insertions) == "AA"


def test_find_indels_deletion_in_novel():
    ref = "ACGTACGTACGTACGT"
    novel = "ACGTACGTCGTACGT"
    fixed, insertions = nan.find_indels(ref, novel)
    assert insertions == []
    assert "-" in fixed
    assert fixed.replace("-", "") == novel


def test_find_indels_identical():
    fixed, ins = nan.find_indels("ACGTACGT", "ACGTACGT")
    assert fixed == "ACGTACGT" and ins == []


def test_build_ambiguous_ref_unique_sequences_keep_names():
    ref = {"IGHV1-1*01": "AAAACCCC", "IGHV1-2*01": "AAAAGGGG"}
    res = nan.build_ambiguous_ref(ref, 2)
    assert res == {"IGHV1-1*01": "..AACCCC", "IGHV1-2*01": "..AAGGGG"}


def test_build_ambiguous_ref_collapses_identical_alleles_of_one_gene():
    ref = {"IGHV1-1*01": "AAAACCCC", "IGHV1-1*02": "GGAACCCC"}
    assert nan.build_ambiguous_ref(ref, 2) == {"IGHV1-1*01_02": "..AACCCC"}


def test_build_ambiguous_ref_collapses_genes_same_family():
    ref = {"IGHV1-1*01": "AAAACCCC", "IGHV1-2*01": "GGAACCCC", "IGHV1-3*02": "TTAACCCC"}
    assert nan.build_ambiguous_ref(ref, 2) == {"IGHV1-1*01_2.01_3.02": "..AACCCC"}


def test_build_ambiguous_ref_collapses_genes_different_family():
    ref = {"IGHV1-1*01": "AAAACCCC", "IGHV2-5*01": "GGAACCCC"}
    res = nan.build_ambiguous_ref(ref, 2)
    assert list(res.values()) == ["..AACCCC"]
    name = list(res)[0]
    assert name.startswith("IGHV1-1*01_") and "01" in name
