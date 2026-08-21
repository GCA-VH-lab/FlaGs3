import os
import shutil
import subprocess
import tempfile
from io import StringIO
from typing import Dict, List, Tuple

from Bio import Phylo, SeqIO



class TreeBuilder:

  GAP_CHARS = "-."

  @staticmethod
  def _debug(message, exc=False):
    try:
      from flags_log import debug
      debug(message, exc=exc)
    except ImportError:
      pass

  @classmethod
  def _run(cls, cmd, **kwargs):
    cls._debug("running: {}".format(" ".join(str(c) for c in cmd)))
    result = subprocess.run(cmd, **kwargs)
    cls._debug("exit {} from {}".format(result.returncode, cmd[0]))
    return result

  def __init__(self, threads: int = 0, engine: str = "veryfasttree",
      gap_threshold: float = 0.1):
    self.threads = threads
    self.engine = engine
    self.gap_threshold = gap_threshold
    self.alignment: Dict[str, str] = {}

  def build(self, sequences: Dict[str, str]) -> Tuple[str, List[str]]:
    names = list(sequences)
    if len(names) < 3:
      return "", names

    try:
      with tempfile.TemporaryDirectory() as tmp:
        fasta = os.path.join(tmp, "q.fasta")
        aln = os.path.join(tmp, "q.aln")
        with open(fasta, "w") as out:
          for name, seq in sequences.items():
            out.write(">{}\n{}\n".format(name, seq))

        with open(aln, "w") as out:
          self._run(["mafft", "--auto", "--anysymbol", "--quiet",
          "--thread", str(self.threads), fasta],
          check=True, stdout=out, stderr=subprocess.DEVNULL)

        raw = self._read_alignment(aln)
        self.alignment = self._trim(raw, self.gap_threshold)
        if raw:
          self._debug("alignment {} cols -> {} after trimming at gt={}".format(
            len(next(iter(raw.values()))),
            len(next(iter(self.alignment.values()))), self.gap_threshold))
        trimmed = os.path.join(tmp, "q.trimmed.aln")
        with open(trimmed, "w") as out:
          for name, seq in self.alignment.items():
            out.write(">{}\n{}\n".format(name, seq))

        if self.engine == "iqtree":
          newick = self._run_iqtree(trimmed, tmp, len(names))
        else:
          newick = self._run(["VeryFastTree", trimmed],
            check=True, capture_output=True, text=True).stdout.strip()
    except FileNotFoundError as e:
      print("Warning: tree building needs mafft and {} on PATH; skipping the tree "
            "({}).".format("iqtree" if self.engine == "iqtree" else "VeryFastTree", e))
      return "", names
    except subprocess.CalledProcessError as e:
      print("Warning: tree building failed, skipping the tree ({}).".format(e))
      return "", names

    leaf_order = [t.name for t in Phylo.read(StringIO(newick), "newick").get_terminals()]
    return newick, leaf_order

  def _run_iqtree(self, aln: str, tmp: str, n_taxa: int) -> str:
    binary = next((b for b in ("iqtree3", "iqtree2", "iqtree")
                   if shutil.which(b)), None)
    if binary is None:
      raise FileNotFoundError("no iqtree binary found")
    cmd = [binary, "-s", aln, "-m", "MFP", "--prefix", os.path.join(tmp, "iq"),
           "-T", str(self.threads) if self.threads else "AUTO", "--quiet"]
    if n_taxa >= 4:
      cmd += ["-B", "1000"]   # ultrafast bootstrap needs at least 4 taxa
    self._run(cmd, check=True, capture_output=True, text=True)
    with open(os.path.join(tmp, "iq.treefile")) as fh:
      return fh.read().strip()

  @staticmethod
  def _read_alignment(path: str) -> Dict[str, str]:
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(path, "fasta")}

  @classmethod
  def _trim(cls, alignment: Dict[str, str], gap_threshold: float) -> Dict[str, str]:
    rows = list(alignment.values())
    if not rows:
      return alignment
    width = len(rows[0])
    need = gap_threshold * len(rows)
    keep = [i for i in range(width)
            if sum(1 for r in rows if r[i] not in cls.GAP_CHARS) >= need]
    if not keep or len(keep) == width:
      return alignment
    return {name: "".join(seq[i] for i in keep) for name, seq in alignment.items()}


def ladderized_leaf_order(newick):
  t = Phylo.read(StringIO(newick), "newick")
  try:
    t.root_at_midpoint()
  except Exception:
    pass
  t.ladderize()
  return [tip.name for tip in t.get_terminals()]