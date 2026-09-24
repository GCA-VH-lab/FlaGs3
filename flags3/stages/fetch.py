from pathlib import Path
from typing import Optional

from flags3 import fasta, mgnify, ncbi, net
from flags3.inputs import InputList
from flags3.log import note
from flags3.schema import MISSING, Failure, GenomeFiles, QueryTarget
from flags3.stage import Stage

SLOTS = ("gff", "faa", "rna", "genome")
EXTENSIONS = {
	"rna": (".rna.fna", ".rna.fna.gz", ".rna.fa", ".rna.fa.gz", "_rna_from_genomic.fna", "_rna_from_genomic.fna.gz"),
	"gff": (".gff", ".gff3", ".gff.gz", ".gff3.gz"),
	"faa": (".faa", ".faa.gz", ".fasta", ".fasta.gz", ".fa", ".fa.gz"),
	"genome": (".fna", ".fna.gz"),
}
INFIX = ("_genomic", "_protein", "_rna_from_genomic", "_cds_from_genomic")


def _basename(name: str, exts) -> Optional[str]:
	for ext in sorted(exts, key=len, reverse=True):
		if name.endswith(ext):
			stem = name[:-len(ext)]
			for infix in INFIX:
				if stem.endswith(infix):
					stem = stem[:-len(infix)]
					break
			return stem
	return None


class GenomeCache:
	def __init__(self, directory: Path):
		self.directory = Path(directory)
		self.directory.mkdir(parents=True, exist_ok=True)
		self.index: dict[str, dict[str, Path]] = {}
		self._protein_index: Optional[dict[str, str]] = None
		self.scan()

	def scan(self) -> None:
		self.index = {}
		claimed = set()
		names = sorted(p.name for p in self.directory.iterdir() if p.is_file() and not p.name.endswith(".part"))
		for slot in ("rna", "gff", "faa", "genome"):
			for name in names:
				if name in claimed:
					continue
				base = _basename(name, EXTENSIONS[slot])
				if base is not None:
					self.index.setdefault(base, {})[slot] = self.directory / name
					claimed.add(name)

	def find(self, assembly: str) -> Optional[str]:
		if assembly in self.index:
			return assembly
		for base in self.index:
			if base.startswith(assembly) or assembly.startswith(base):
				return base
		return None

	def files(self, base: str, slots: list[str]) -> dict[str, Path]:
		good = {}
		for slot in slots:
			path = self.index.get(base, {}).get(slot)
			if path is None:
				continue
			if net.intact(path):
				good[slot] = path
			else:
				path.unlink(missing_ok=True)
				del self.index[base][slot]
		return good

	def by_protein(self, protein: str) -> Optional[str]:
		if self._protein_index is None:
			self._protein_index = {}
			for base, slots in self.index.items():
				if "faa" in slots and net.intact(slots["faa"]):
					for name, _, _ in fasta.read(slots["faa"]):
						self._protein_index.setdefault(name, base)
		return self._protein_index.get(protein)


class Fetch(Stage):
	name = "fetch"

	def run(self, run, config, out: Path) -> None:
		inputs = InputList()
		for name in config.text("inputs").split(","):
			if name:
				inputs.read(run.input_dir / name)
		for warning in inputs.warnings:
			print("Warning: " + warning)
		if run.has("blast"):
			known = {p.split(".")[0] for p, _ in inputs.entries}
			for line in run.stage_file("blast", "accessions.txt").read_text().splitlines():
				accession = line.strip()
				if accession and not accession.startswith("#") and accession.split(".")[0] not in known:
					known.add(accession.split(".")[0])
					inputs.entries.append((accession, None))
					inputs.unpaired.append(accession)
		slots = ["gff", "faa"] + [s for s in config.text("fetch.slots").split(",") if s in ("rna", "genome")]
		offline = config.flag("offline")
		cache = GenomeCache(config.path("genomes"))
		failures: list[Failure] = []
		targets = self._resolve(inputs, config, cache, offline, failures)
		genomes = self._collect(targets, slots, config, cache, offline, failures)
		for i, t in enumerate(targets):
			if t.status == "ok" and not genomes.get(t.assembly, GenomeFiles(t.assembly, "", MISSING, MISSING, MISSING, MISSING)).usable:
				targets[i] = QueryTarget(t.query, t.assembly, t.accessions, "unresolved: genome files unavailable")
		GenomeFiles.write(out / GenomeFiles.FILE, (genomes[a] for a in sorted(genomes)))
		QueryTarget.write(out / QueryTarget.FILE, targets)
		Failure.write(out / Failure.FILE, failures)
		resolved = sum(1 for t in targets if t.status == "ok")
		note("resolved {} of {} queries to {} genomes".format(resolved, len(targets), len(genomes)))
		if targets and not resolved:
			raise RuntimeError("no query could be resolved to a genome; see fetch/{}".format(Failure.FILE))

	def _mapper(self, config) -> ncbi.IpgMapper:
		email = config.text("user_email")
		if not email:
			raise RuntimeError("NCBI resolution needs -u/--user_email; use --offline for genomes you supply yourself")
		return ncbi.IpgMapper(email, api_key=config.text("api_key") or None,
			max_assemblies=config.integer("max_assemblies", 1), cross_db=not config.flag("no_cross_db"))

	def _resolve(self, inputs: InputList, config, cache: GenomeCache, offline: bool,
			failures: list[Failure]) -> list[QueryTarget]:
		targets: list[QueryTarget] = []
		if offline:
			for protein, assembly in inputs.entries:
				base = cache.find(assembly) if assembly else cache.by_protein(protein)
				reason = "unresolved: not in the genome directory" if assembly else "unresolved: not in any genome in the directory"
				targets.append(QueryTarget(protein, base or assembly or MISSING, MISSING, "ok" if base else reason))
			return targets
		local_first = config.flag("local_first")
		found_locally: dict[str, str] = {}
		if local_first:
			for protein in inputs.unpaired:
				base = cache.by_protein(protein)
				if base:
					found_locally[protein] = base
			if found_locally:
				note("{} bare queries found in the genome directory; the rest go to IPG".format(len(found_locally)))
		mapper = self._mapper(config) if any(p not in found_locally for p in inputs.unpaired) or (config.flag("remap") and inputs.paired) else None
		remap = config.flag("remap")
		lookup = [p for p in inputs.unpaired if p not in found_locally] + ([p for p, _ in inputs.paired] if remap else [])
		result = mapper.map(lookup) if (lookup and mapper) else ncbi.Resolution()
		for protein, reason in result.errors.items():
			failures.append(Failure(protein, reason))
		for protein, dropped in result.dropped.items():
			failures.append(Failure(protein, "cross-database assemblies excluded: " + ",".join(sorted(dropped))))
		for protein, assembly in inputs.entries:
			if assembly:
				targets.append(self._paired(protein, assembly, result, remap, failures))
				continue
			if protein in found_locally:
				targets.append(QueryTarget(protein, found_locally[protein], MISSING, "ok"))
				continue
			chosen = result.assemblies.get(protein) or []
			if not chosen:
				reason = result.errors.get(protein, "unresolved: IPG lists no assembly")
				targets.append(QueryTarget(protein, MISSING, MISSING, reason if reason.startswith("unresolved") else "unresolved: " + reason))
				continue
			for choice in chosen:
				targets.append(QueryTarget(protein, choice, ",".join(sorted(result.aliases(protein, choice))) or MISSING, "ok"))
		return targets

	@staticmethod
	def _paired(protein, assembly, result, remap, failures) -> QueryTarget:
		aliases = result.aliases(protein, assembly)
		if remap and protein in result.accessions and not aliases:
			choice = result.assemblies.get(protein) or []
			if choice:
				failures.append(Failure(protein, "remapped from {} to {} by IPG".format(assembly, choice[0])))
				return QueryTarget(protein, choice[0], ",".join(sorted(result.aliases(protein, choice[0]))) or MISSING, "ok")
		return QueryTarget(protein, assembly, ",".join(sorted(aliases)) or MISSING, "ok")
		local_first = config.flag("local_first")
		found_locally: dict[str, str] = {}
		if local_first:
			for protein in inputs.unpaired:
				base = cache.by_protein(protein)
				if base:
					found_locally[protein] = base
			if found_locally:
				note("{} bare queries found in the genome directory; the rest go to IPG".format(len(found_locally)))
		mapper = self._mapper(config) if any(p not in found_locally for p in inputs.unpaired) or (config.flag("remap") and inputs.paired) else None
		remap = config.flag("remap")
		lookup = [p for p in inputs.unpaired if p not in found_locally] + ([p for p, _ in inputs.paired] if remap else [])
		result = mapper.map(lookup) if (lookup and mapper) else ncbi.Resolution()
		for protein, reason in result.errors.items():
			failures.append(Failure(protein, reason))
		for protein, dropped in result.dropped.items():
			failures.append(Failure(protein, "cross-database assemblies excluded: " + ",".join(sorted(dropped))))
		for protein, assembly in inputs.paired:
			aliases = result.aliases(protein, assembly)
			if remap and protein in result.accessions and not aliases:
				choice = result.assemblies.get(protein) or []
				if choice:
					targets.append(QueryTarget(protein, choice[0], ",".join(sorted(result.aliases(protein, choice[0]))) or MISSING,
						"ok"))
					failures.append(Failure(protein, "remapped from {} to {} by IPG".format(assembly, choice[0])))
					continue
			targets.append(QueryTarget(protein, assembly, ",".join(sorted(aliases)) or MISSING, "ok"))
		for protein in inputs.unpaired:
			chosen = result.assemblies.get(protein) or []
			if not chosen:
				reason = result.errors.get(protein, "unresolved: IPG lists no assembly")
				targets.append(QueryTarget(protein, MISSING, MISSING, reason if reason.startswith("unresolved") else "unresolved: " + reason))
				continue
			for assembly in chosen:
				targets.append(QueryTarget(protein, assembly, ",".join(sorted(result.aliases(protein, assembly))) or MISSING, "ok"))
		return targets

	def _collect(self, targets: list[QueryTarget], slots: list[str], config, cache: GenomeCache,
			offline: bool, failures: list[Failure]) -> dict[str, GenomeFiles]:
		found: dict[str, dict[str, Path]] = {}
		source: dict[str, str] = {}
		wanted: dict[str, list[str]] = {}
		for i, t in enumerate(targets):
			if t.status != "ok":
				continue
			base = cache.find(t.assembly)
			if base and base != t.assembly:
				targets[i] = t = QueryTarget(t.query, base, t.accessions, t.status)
			if t.assembly in found:
				continue
			have = cache.files(base, slots) if base else {}
			found[t.assembly] = have
			source[t.assembly] = "local" if offline else "cache"
			missing = [s for s in slots if s not in have and not (mgnify.is_mgnify(t.assembly) and s == "rna")]
			if missing and not offline:
				wanted[t.assembly] = missing
		if wanted:
			self._download(wanted, config, cache, found, source, failures)
		genomes = {}
		for assembly, have in found.items():
			genomes[assembly] = GenomeFiles(assembly, source[assembly],
				*(str(have[s]) if s in have else MISSING for s in SLOTS))
		return genomes

	def _download(self, wanted, config, cache, found, source, failures):
		rate = 10.0 if config.text("api_key") else 5.0
		workers = min(config.workers(), 10)
		groups = {
			"mgnify": {a: s for a, s in wanted.items() if mgnify.is_mgnify(a)},
			"ncbi": {a: s for a, s in wanted.items() if not mgnify.is_mgnify(a)},
		}
		for label, group in groups.items():
			if not group:
				continue
			client = (mgnify.MgnifyGenomes if label == "mgnify" else ncbi.NcbiGenomes)(cache.directory, rate, workers)
			note("downloading {} genomes from {}".format(len(group), label))
			progress = (lambda done, total: note("{} download: {}/{}".format(label, done, total))) if len(group) > 3 else None
			for assembly, got in client.fetch_many(group, progress).items():
				if got:
					found[assembly].update(got)
					source[assembly] = label
			for subject, reason in client.failures.items():
				failures.append(Failure(subject, reason))
		cache.scan()
