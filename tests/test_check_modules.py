import runpy
import sys

import pytest

AUX_HEADER = "#name\tj_codon_frame\tchain_type\tj_cdr3_end\textra_bps\n"
NDM_COLS = ['#name', 'FWR1 start', 'FWR1 stop', 'CDR1 start', 'CDR1 stop', 'FWR2 start', 'FWR2 stop',
            'CDR2 start', 'CDR2 stop', 'FWR3 start', 'FWR3 stop', 'chain type', 'frame start']


def run_module(monkeypatch, name, args):
    monkeypatch.setattr(sys, "argv", [name] + [str(a) for a in args])
    runpy.run_module(f"receptor_utils.{name}", run_name="__main__")


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


# ---- check_aux ----

def test_check_aux_matching(monkeypatch, tmp_path, capsys):
    aux = write(tmp_path, "a.aux", AUX_HEADER + "IGHJ4*01\t0\tJH\t18\t0\n")
    ig = write(tmp_path, "ig.aux", "IGHJ4*01 0 JH 18 0\n")
    run_module(monkeypatch, "check_aux", [aux, ig])
    assert capsys.readouterr().out == ""


def test_check_aux_reports_differences(monkeypatch, tmp_path, capsys):
    aux = write(tmp_path, "a.aux", AUX_HEADER + "IGHJ4*01\t0\tJH\t18\t0\nIGHJ5*01\t1\tJH\t18\t0\n")
    ig = write(tmp_path, "ig.aux", "IGHJ4*01 1 JK 19 2\n")
    run_module(monkeypatch, "check_aux", [aux, ig])
    out = capsys.readouterr().out
    assert "IGHJ4*01 j_codon_frame 0 != 1" in out
    assert "IGHJ4*01 chain_type JH != JK" in out
    assert "IGHJ4*01 j_cdr3_end 18 != 19" in out
    assert "IGHJ4*01 extra_bps 0 != 2" in out
    assert "IGHJ5*01 not found in igblast" in out


def test_check_aux_name_without_allele(monkeypatch, tmp_path, capsys):
    aux = write(tmp_path, "a.aux", AUX_HEADER + "IGHJ-4\t0\tJH\t18\t0\n")
    ig = write(tmp_path, "ig.aux", "IGHJ0-4*00 0 JH 18 0\n")
    run_module(monkeypatch, "check_aux", [aux, ig])
    assert capsys.readouterr().out == ""


def test_check_aux_real_files(monkeypatch, data_dir, capsys):
    run_module(monkeypatch, "check_aux", [data_dir / "Human_IGH_VDJ_rev_2.aux", data_dir / "Human_IGH_VDJ_rev_2.aux"])
    assert capsys.readouterr().out == ""


def test_check_aux_missing_args(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["check_aux"])
    with pytest.raises(SystemExit) as e:
        runpy.run_module("receptor_utils.check_aux", run_name="__main__")
    assert e.value.code == 2


# ---- check_fastas ----

def test_check_fastas(monkeypatch, write_fasta_file, capsys):
    f1 = write_fasta_file("f1.fasta", {"same": "ACGT", "long": "ACGTAA", "short": "AC",
                                       "diff": "ACGA", "new": "AAAA"})
    f2 = write_fasta_file("f2.fasta", {"same": "ACGT", "long": "ACGT", "short": "ACGT", "diff": "ACGT"})
    run_module(monkeypatch, "check_fastas", [f1, f2])
    out = capsys.readouterr().out.splitlines()
    assert "long sequence longer than reference" in out
    assert "short sequence shorter than reference" in out
    assert "diff sequence mismatch" in out
    assert "new not found in reference" in out
    assert not any(l.startswith("same") for l in out)
    assert "long sequence mismatch" not in out


def test_check_fastas_identical_files(monkeypatch, data_dir, capsys):
    run_module(monkeypatch, "check_fastas", [data_dir / "V_gapped.fasta", data_dir / "V_gapped.fasta"])
    assert capsys.readouterr().out == ""


# ---- check_ndm ----

def test_check_ndm_no_errors(monkeypatch, tmp_path, capsys):
    vals = [str(i) for i in range(1, 12)] + ["0"]
    ndm = write(tmp_path, "a.ndm", "\t".join(NDM_COLS) + "\n" + "IGHV1-1*01\t" + "\t".join(vals) + "\n")
    ig = write(tmp_path, "ig.ndm", "IGHV1-1*01 " + " ".join(vals) + "\n")
    run_module(monkeypatch, "check_ndm", [ndm, ig])
    assert "No errors found" in capsys.readouterr().out


def test_check_ndm_errors_and_fallback(monkeypatch, tmp_path, capsys):
    vals = [str(i) for i in range(1, 13)]
    changed = list(vals)
    changed[0] = "99"
    ndm = write(tmp_path, "a.ndm", "\t".join(NDM_COLS) + "\n"
                + "IGHV1-1*01\t" + "\t".join(vals) + "\n"
                + "IGHV1-1*05\t" + "\t".join(vals) + "\n"
                + "IGHV9-9*01\t" + "\t".join(vals) + "\n")
    ig = write(tmp_path, "ig.ndm", "IGHV1-1*01 " + " ".join(changed) + "\n")
    run_module(monkeypatch, "check_ndm", [ndm, ig])
    out = capsys.readouterr().out
    assert "FWR1 start 1 != 99" in out
    assert "IGHV1-1*05 not found in igblast: using IGHV1-1*01" in out
    assert "IGHV9-9*01 not found in igblast. *01 allele also not found: skipping" in out
    assert "2 errors found" in out
