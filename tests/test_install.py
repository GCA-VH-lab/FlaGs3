import gzip
import tarfile
from pathlib import Path

import pytest

from flags3 import home, install
from flags3.tools import Tools
from tests.test_domains import _hmm


@pytest.fixture
def fake_home(tmp_path, monkeypatch):
	root = tmp_path / "home"
	monkeypatch.setattr(home, "HOME", root)
	monkeypatch.setattr(home, "GENOMES", root / "genomes")
	monkeypatch.setattr(home, "DB", root / "db")
	monkeypatch.setattr(home, "TOOLS", root / "tools")
	monkeypatch.setattr(home, "ENVS", root / "envs")
	monkeypatch.setattr(home, "USER_TOOLS_TABLE", root / "tools_table.tsv")
	monkeypatch.setattr(home, "PFAM_HMM", root / "db" / "pfam" / "Pfam-A.hmm")
	monkeypatch.setattr(home, "PFAM_CLANS", root / "db" / "pfam" / "Pfam-A.clans.tsv.gz")
	monkeypatch.setattr(home, "INTERPRO", root / "db" / "interpro" / "interpro_metadata_processed.tsv")
	return root


def test_micromamba_asset_names(monkeypatch):
	monkeypatch.setattr(install.platform, "system", lambda: "Linux")
	monkeypatch.setattr(install.platform, "machine", lambda: "x86_64")
	assert install.Micromamba.asset() == "linux-64"
	monkeypatch.setattr(install.platform, "machine", lambda: "aarch64")
	assert install.Micromamba.asset() == "linux-aarch64"
	monkeypatch.setattr(install.platform, "system", lambda: "Windows")
	with pytest.raises(install.InstallError):
		install.Micromamba.asset()


def test_listing_and_unknown(fake_home):
	installer = install.Installer()
	names = [n for n, _, _ in installer.status()]
	assert names == list(install.BY_NAME) and "deeptmhmm" in names
	assert all(state == "-" for _, state, _ in installer.status())
	with pytest.raises(install.InstallError):
		installer.run(["nothing"], {})


def test_interpro_copies_the_file(fake_home, tmp_path):
	table = tmp_path / "interpro_metadata_processed.tsv"
	table.write_text("accession\tpfam_members\n")
	installer = install.Installer()
	installer.run(["interpro"], {"interpro": table})
	assert home.INTERPRO.read_text() == table.read_text()
	assert installer.components["interpro"].installed() == "installed"
	assert home.USER_TOOLS_TABLE.is_file()
	with pytest.raises(install.InstallError):
		install.Installer().run(["interpro"], {}, force=True)
	install.Installer().run(["interpro"], {})


def test_pfam_from_local_archive_presses(fake_home, tmp_path):
	hmm = _hmm(tmp_path / "tiny.hmm")
	archive = tmp_path / "Pfam-A.hmm.gz"
	with open(hmm, "rb") as src, gzip.open(archive, "wb") as dst:
		dst.write(src.read())
	installer = install.Installer()
	assert installer.components["pfam"].installed() is None
	installer.run(["pfam"], {"pfam": archive})
	assert home.PFAM_HMM.is_file()
	assert home.PFAM_HMM.with_name("Pfam-A.hmm.h3m").is_file()
	assert installer.components["pfam"].installed() == "installed"


def test_defence_hmm_from_local_tarball(fake_home, tmp_path):
	hmm = _hmm(tmp_path / "Cas__Cas9.hmm")
	tar = tmp_path / "defense-finder-models-v2.tar.gz"
	with tarfile.open(tar, "w:gz") as out:
		out.add(hmm, arcname="defense-finder-models-v2/profiles/Cas__Cas9.hmm")
		out.add(hmm, arcname="defense-finder-models-v2/definitions/x.hmm")
	installer = install.Installer()
	installer.run(["defence-hmm"], {"defence-hmm": tar})
	profiles = home.DB / "defensefinder" / "profiles"
	assert (profiles / "Cas__Cas9.hmm").is_file()
	assert installer.components["defence-hmm"].installed() == "installed"


def test_core_builds_missing_and_writes_rows(fake_home, monkeypatch):
	created = []

	def fake_create(self, prefix, packages):
		created.append((prefix, packages))
		(prefix / "bin").mkdir(parents=True, exist_ok=True)
		for name in ("mafft", "trimal", "VeryFastTree", "iqtree", "blastp"):
			(prefix / "bin" / name).write_text("")
		return prefix

	monkeypatch.setattr(install.Micromamba, "create", fake_create)
	monkeypatch.setattr(install.shutil, "which", lambda name: None)
	installer = install.Installer()
	installer.run(["core"], {})
	assert created[0][1] == install.Core.packages
	tools = Tools.load(home.USER_TOOLS_TABLE)
	assert tools["mafft"].command.startswith(str(home.ENVS / "core" / "bin" / "mafft") + " --auto")
	assert tools["blastp"].command.startswith(str(home.ENVS / "core" / "bin" / "blastp") + " -query")
	assert tools["jackhmmer"].engine == "jackhmmer"
	assert install.Installer().components["core"].installed() == "installed"


def test_core_uses_tools_on_path(fake_home, monkeypatch, tmp_path):
	on_path = tmp_path / "usr" / "bin"
	on_path.mkdir(parents=True)
	for name in ("mafft", "trimal", "VeryFastTree", "iqtree", "blastp"):
		(on_path / name).write_text("")
	monkeypatch.setattr(install.shutil, "which", lambda name: str(on_path / name) if (on_path / name).exists() else None)
	monkeypatch.setattr(install.Micromamba, "create", lambda self, prefix, packages: pytest.fail("should not build"))
	installer = install.Installer()
	assert installer.components["core"].installed() == "on PATH"
	installer.run(["core"], {})
	tools = Tools.load(home.USER_TOOLS_TABLE)
	assert tools["trimal"].command.startswith(str(on_path / "trimal"))


def test_mmseqs_from_release_archive(fake_home, tmp_path, monkeypatch):
	monkeypatch.setattr(install.shutil, "which", lambda name: None)
	tar = tmp_path / "mmseqs-linux-avx2.tar.gz"
	with tarfile.open(tar, "w:gz") as out:
		script = tmp_path / "mmseqs"
		script.write_text("#!/bin/sh\n")
		out.add(script, arcname="mmseqs/bin/mmseqs")
	installer = install.Installer()
	assert installer.components["mmseqs"].installed() is None
	installer.run(["mmseqs"], {"mmseqs": tar})
	assert (home.TOOLS / "mmseqs" / "bin" / "mmseqs").is_file()
	tools = Tools.load(home.USER_TOOLS_TABLE)
	assert tools["mmseqs_cluster"].command.startswith(str(home.TOOLS / "mmseqs" / "bin" / "mmseqs") + " easy-search")
	assert tools["mmseqs_cluster"].options["sensitivity"] == "7.5"
	assert install.Installer().components["mmseqs"].installed() == "installed"
	monkeypatch.setattr(install.Mmseqs, "install", lambda self, source=None: pytest.fail("reinstalled"))
	install.Installer().run(["mmseqs"], {})
	monkeypatch.setattr(install.platform, "system", lambda: "Darwin")
	assert install.Mmseqs.asset() == "mmseqs-osx-universal.tar.gz"


def test_micromamba_prefers_existing(fake_home, monkeypatch, tmp_path):
	existing = tmp_path / "micromamba"
	existing.write_text("")
	monkeypatch.setattr(install.shutil, "which", lambda name: str(existing) if name == "micromamba" else None)
	m = install.Micromamba()
	assert m.ensure() == existing
	assert "DYLD_LIBRARY_PATH" not in m.env() and m.env()["MAMBA_ROOT_PREFIX"] == str(home.HOME / "mamba")
