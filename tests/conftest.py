import importlib
import sys
from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent


@pytest.fixture(scope="session")
def data_dir():
    """Directory holding the sample data files shipped in tests/."""
    return DATA_DIR


@pytest.fixture
def write_fasta_file(tmp_path):
    """Return a helper that writes a dict of name -> sequence to a FASTA file in tmp_path."""
    def _write(name, seqs):
        path = tmp_path / name
        with open(path, "w") as f:
            for k, v in seqs.items():
                f.write(f">{k}\n{v}\n")
        return str(path)
    return _write


@pytest.fixture
def run_script(monkeypatch):
    """Run a script's main() with the given argv. Returns the SystemExit code (0 if main returned normally)."""
    def _run(module_name, args):
        mod = importlib.import_module(f"scripts.{module_name}")
        monkeypatch.setattr(sys, "argv", [module_name] + [str(a) for a in args])
        try:
            mod.main()
        except SystemExit as e:
            return e.code if e.code is not None else 0
        return 0
    return _run


# Small synthetic germline references
V_UNGAPPED = {
    "IGHV1-1*01": "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTTTCCTGCAAGGCATCTGGATACACCTTCACC",
    "IGHV1-1*02": "CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTTTCCTGCAAGGCATCTGGATACACCTTCAGC",
}

J_UNGAPPED = {
    "IGHJ4*01": "ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG",
    "IGHJ4*02": "ACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG",
}

D_UNGAPPED = {
    "IGHD1-1*01": "GGTACAACTGGAACGAC",
    "IGHD2-2*01": "AGGATATTGTAGTAGTACCAGCTGCTATGCC",
}


@pytest.fixture
def v_seqs():
    return dict(V_UNGAPPED)


@pytest.fixture
def j_seqs():
    return dict(J_UNGAPPED)


@pytest.fixture
def d_seqs():
    return dict(D_UNGAPPED)
