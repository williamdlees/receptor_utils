import configparser
import importlib
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def entry_points():
    cfg = configparser.ConfigParser()
    cfg.read(ROOT / "setup.cfg")
    lines = cfg["options.entry_points"]["console_scripts"].strip().splitlines()
    return dict(l.strip().split(" = ") for l in lines if l.strip())


EP = entry_points()


def test_fourteen_console_scripts():
    assert len(EP) == 14


def ep_params():
    return sorted(EP.items())


@pytest.mark.parametrize("name,target", ep_params())
def test_entry_point_importable(name, target):
    module_name, func = target.split(":")
    module = importlib.import_module(module_name)
    assert callable(getattr(module, func))


@pytest.mark.parametrize("name,target", ep_params())
def test_entry_point_help(name, target, monkeypatch, capsys):
    module = importlib.import_module(target.split(":")[0])
    monkeypatch.setattr(sys, "argv", [name, "--help"])
    with pytest.raises(SystemExit) as e:
        module.main()
    assert e.value.code == 0
    assert "usage" in capsys.readouterr().out.lower()


def test_installed_script_runs():
    exe = Path(sys.executable).parent / ("rev_comp.exe" if sys.platform == "win32" else "rev_comp")
    if not exe.exists():
        pytest.skip("package not installed with console scripts")
    out = subprocess.run([str(exe), "AACGT"], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "ACGTT"
