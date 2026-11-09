import json

import pytest
import requests

import scripts.download_germline_set as dgs
from receptor_utils import simple_bio_seq as simple

from test_aux_formats import germline_set, j_allele, v_allele, d_allele

BASE = "https://ogrdb.test/api"
SPECIES = {"species": [{"label": "Homo sapiens", "id": "NCBITAXON:9606"}, {"label": "Mus musculus", "id": "10090"}]}
SETS = {"germline_species": [
    {"locus": "IGH", "germline_set_name": "IGH_VDJ", "germline_set_id": "S1"},
    {"locus": "IGK", "germline_set_name": "KA", "germline_set_id": "S2"},
    {"locus": "IGK", "germline_set_name": "KB", "germline_set_id": "S3"},
]}
GAPPED = ">IGHV1-1*01\nAAA...CCC\n>IGHD1-1*01\nGGGTTT\n>IGHJ4*01\nTTTGGG\n>IGHG1*01\n" + "CCCC" + "\n"
UNGAPPED = ">IGHV1-1*01\nAAACCC\n>IGHD1-1*01\nGGGTTT\n>IGHJ4*01\nTTTGGG\n"


def ref_set():
    gs = germline_set([v_allele(), d_allele(), j_allele()])
    gs["GermlineSet"][0]["species"] = {"id": "NCBITAXON:9606", "label": "Homo sapiens"}
    return gs


class FakeResponse:
    def __init__(self, payload=None, text="", status=200):
        self.payload, self.text, self.status_code = payload, text, status

    def json(self):
        return self.payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.exceptions.HTTPError(f"{self.status_code} error")


@pytest.fixture
def api(monkeypatch, tmp_path):
    calls = []

    def fake_get(url, *a, **kw):
        calls.append(url)
        assert url.startswith(BASE), "unexpected non-mocked URL"
        path = url[len(BASE):]
        if path == "/germline/species":
            return FakeResponse(SPECIES)
        if path == "/germline/sets/NCBITAXON:9606":
            return FakeResponse(SETS)
        if path == "/germline/set/S1/versions":
            return FakeResponse({"versions": [1, 2]})
        if path.startswith("/germline/set/S1/") and path.endswith("/gapped_ex"):
            return FakeResponse(text=GAPPED)
        if path.startswith("/germline/set/S1/") and path.endswith("/ungapped_ex"):
            return FakeResponse(text=UNGAPPED)
        if path.startswith("/germline/set/S1/"):
            return FakeResponse(ref_set())
        return FakeResponse(status=404, text="not found")

    monkeypatch.setattr(requests, "get", fake_get)
    monkeypatch.chdir(tmp_path)
    return calls


def run(run_script, *args):
    return run_script("download_germline_set", ["Homo sapiens", "IGH", "-u", BASE] + list(args))


# ---- helpers ----

def test_do_request_http_error(monkeypatch, capsys):
    monkeypatch.setattr(requests, "get", lambda url: FakeResponse(status=500, text="boom"))
    with pytest.raises(SystemExit):
        dgs.do_request("http://x")
    assert "HTTP error occurred" in capsys.readouterr().out


def test_do_request_connection_error(monkeypatch, capsys):
    def boom(url):
        raise requests.exceptions.ConnectionError("down")
    monkeypatch.setattr(requests, "get", boom)
    with pytest.raises(SystemExit):
        dgs.do_request("http://x")
    assert "Request error occurred" in capsys.readouterr().out


def test_find_species_id(api):
    assert dgs.find_species_id(BASE, "Homo sapiens") == "NCBITAXON:9606"


def test_find_species_id_unknown(api, capsys):
    assert dgs.find_species_id(BASE, "Dragon") is None
    assert "species not found" in capsys.readouterr().out


def test_find_available_sets_single(api):
    assert dgs.find_available_sets(BASE, "NCBITAXON:9606", "IGH", None) == "S1"


def test_find_available_sets_by_name(api):
    assert dgs.find_available_sets(BASE, "NCBITAXON:9606", "IGK", "KB") == "S3"


def test_find_available_sets_ambiguous(api, capsys):
    assert dgs.find_available_sets(BASE, "NCBITAXON:9606", "IGK", None) is None
    assert "Multiple sets" in capsys.readouterr().out


def test_find_available_sets_wrong_name(api, capsys):
    assert dgs.find_available_sets(BASE, "NCBITAXON:9606", "IGK", "ZZ") is None
    assert "not found for species" in capsys.readouterr().out


def test_find_available_sets_missing_locus(api, capsys):
    assert dgs.find_available_sets(BASE, "NCBITAXON:9606", "TRB", None) is None
    assert "Available loci" in capsys.readouterr().out


def test_get_set_versions(api):
    assert dgs.get_set_versions(BASE, "S1") == ["1", "2"]


def test_request_set_specific_and_missing_version(api, capsys):
    assert dgs.request_set(BASE, "S1", "2", "Homo sapiens", None, "IGH", "AIRRC-JSON")["GermlineSet"]
    assert dgs.request_set(BASE, "S1", "9", "Homo sapiens", "IGH_VDJ", "IGH", "AIRRC-JSON") is None
    assert "Version 9 of set" in capsys.readouterr().out


def test_request_set_ex_suffix_only_for_human(api):
    dgs.request_set(BASE, "S1", "latest", "Homo sapiens", None, "IGH", "SINGLE-FG")
    assert api[-1].endswith("/gapped_ex")
    dgs.request_set(BASE, "S1", "latest", "Mus musculus", None, "IGH", "SINGLE-FG")
    assert api[-1].endswith("/gapped")


def test_write_mixcr_to_disk_requires_species(tmp_path, capsys):
    out = tmp_path / "m.json"
    dgs.write_mixcr_to_disk({"GermlineSet": [{"allele_descriptions": []}]}, str(out), "uri")
    assert not out.exists()
    assert "Taxon ID not found" in capsys.readouterr().out


# ---- main() ----

def test_default_json_download(api, run_script, tmp_path):
    assert run(run_script) == 0
    data = json.loads((tmp_path / "Homo_sapiens_IGH.json").read_text())
    assert len(data["GermlineSet"][0]["allele_descriptions"]) == 3


def test_prefix_and_name(api, run_script, tmp_path):
    assert run(run_script, "-n", "IGH_VDJ", "-p", "pre_") == 0
    assert (tmp_path / "pre.json").exists()


def test_specific_version(api, run_script, tmp_path):
    assert run(run_script, "-v", "2") == 0
    assert any(u.endswith("/germline/set/S1/2") for u in api)


def test_gapped_and_ungapped_fasta(api, run_script, tmp_path):
    assert run(run_script, "-f", "SINGLE-FG") == 0
    assert (tmp_path / "Homo_sapiens_IGH_gapped.fasta").read_text() == GAPPED
    assert run(run_script, "-f", "SINGLE-FU") == 0
    assert (tmp_path / "Homo_sapiens_IGH_ungapped.fasta").read_text() == UNGAPPED


def test_mixcr_format(api, run_script, tmp_path):
    assert run(run_script, "-f", "MIXCR") == 0
    data = json.loads((tmp_path / "Homo_sapiens_IGH_mixcr.json").read_text())
    assert data[0]["taxonId"] == "9606"


def test_leader_formats(api, run_script, tmp_path):
    assert run(run_script, "-f", "LEADER") == 0
    assert len(simple.read_fasta(str(tmp_path / "Homo_sapiens_IGH_leaders.fasta"))["IGHV1-1*01"]) == 15
    assert run(run_script, "-f", "LEADER-AA") == 0
    assert len(simple.read_fasta(str(tmp_path / "Homo_sapiens_IGH_leaders_aa.fasta"))["IGHV1-1*01"]) == 5


def test_multi_f(api, run_script, tmp_path, capsys):
    assert run(run_script, "-f", "MULTI-F") == 0
    for n in ("V", "D", "J", "C", "V_gapped"):
        assert (tmp_path / f"Homo_sapiens_IGH_{n}.fasta").exists()
    assert list(simple.read_fasta(str(tmp_path / "Homo_sapiens_IGH_V_gapped.fasta"))) == ["IGHV1-1*01"]
    assert "FASTA files saved to" in capsys.readouterr().out


def test_multi_igblast(api, run_script, tmp_path):
    assert run(run_script, "-f", "MULTI-IGBLAST") == 0
    assert (tmp_path / "Homo_sapiens_IGH_V.fasta").exists()
    assert (tmp_path / "Homo_sapiens_IGH.ndm").exists()
    assert (tmp_path / "Homo_sapiens_IGH.aux").exists()


def test_multi_igblast_prefix_none(api, run_script, tmp_path):
    assert run(run_script, "-f", "MULTI-IGBLAST", "-p", "NONE") == 0
    assert (tmp_path / "V.fasta").exists()
    assert (tmp_path / ".ndm").exists()


def test_unknown_species_downloads_nothing(api, run_script, tmp_path, capsys):
    assert run_script("download_germline_set", ["Dragon", "IGH", "-u", BASE]) == 0
    assert "species not found" in capsys.readouterr().out
    assert list(tmp_path.iterdir()) == []


def test_ambiguous_set_downloads_nothing(api, run_script, tmp_path):
    assert run_script("download_germline_set", ["Homo sapiens", "IGK", "-u", BASE]) == 0
    assert list(tmp_path.iterdir()) == []


def test_http_failure_exits(monkeypatch, run_script, tmp_path):
    monkeypatch.setattr(requests, "get", lambda url: FakeResponse(status=500, text="boom"))
    monkeypatch.chdir(tmp_path)
    assert run(run_script) == 0     # quit() raises SystemExit(None)
    assert list(tmp_path.iterdir()) == []


def test_invalid_format_rejected(run_script):
    assert run_script("download_germline_set", ["Homo sapiens", "IGH", "-f", "BAD"]) == 2


def test_10x_non_human(api, run_script, capsys):
    assert run_script("download_germline_set", ["Mus musculus", "IGH", "-f", "10X", "-u", BASE]) == 0
    assert "only available for human" in capsys.readouterr().out


def test_10x_unknown_locus(api, run_script, capsys):
    assert run_script("download_germline_set", ["Homo sapiens", "XYZ", "-f", "10X", "-u", BASE]) == 0
    assert "Creating 10x" not in capsys.readouterr().out


@pytest.mark.network
def test_live_ogrdb_species_list():
    resp = requests.get("https://ogrdb.airr-community.org/api_v2/germline/species", timeout=30)
    assert resp.status_code == 200
    assert any(s["label"] == "Homo sapiens" for s in resp.json()["species"])
