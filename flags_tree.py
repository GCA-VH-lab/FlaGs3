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
      trimal_mode: str = "gt", trimal_value: float = 0.1,
      trimal_extra: str = ""):
    self.threads = threads
    self.engine = engine
    self.trimal_mode = trimal_mode
    self.trimal_value = trimal_value
    self.trimal_extra = trimal_extra
    self.alignment: Dict[str, str] = {}
    self.raw_alignment: Dict[str, str] = {}
    self.commands: List[str] = []

  def _trimal_cmd(self, src, dst):
    import flags_tools
    mode = (self.trimal_mode or "gt").lstrip("-")
    if mode in ("gappyout", "strict", "strictplus", "automated1", "nogaps", "noallgaps"):
      spec = "-" + mode
    else:
      spec = "-{} {}".format(mode, self.trimal_value)
    if self.trimal_extra:
      spec += " " + self.trimal_extra
    cmd, _ = flags_tools.command("trimal", mode=spec, **{"in": src, "out": dst})
    return cmd

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
          import flags_tools
          cmd, wd = flags_tools.command("mafft", threads=self.threads, **{"in": fasta})
          self._run(cmd, check=True, stdout=out, stderr=subprocess.DEVNULL,
                    cwd=wd or None)

        raw = self._read_alignment(aln)
        self.raw_alignment = raw
        trimmed = os.path.join(tmp, "q.trimmed.aln")
        cmd = self._trimal_cmd(aln, trimmed)
        if shutil.which("trimal"):
          self._run(cmd, check=True, capture_output=True, text=True)
          self.commands.append(" ".join(cmd))
          self.alignment = self._read_alignment(trimmed)
        else:
          print("Warning: trimal not found on PATH; using the untrimmed "
                "alignment.")
          self.alignment = raw
          with open(trimmed, "w") as out:
            for name, seq in raw.items():
              out.write(">{}\n{}\n".format(name, seq))
        if raw and self.alignment:
          self._debug("alignment {} cols -> {} after trimal".format(
            len(next(iter(raw.values()))), len(next(iter(self.alignment.values())))))

        if self.engine == "iqtree":
          newick = self._run_iqtree(trimmed, tmp, len(names))
        else:
          import flags_tools
          vft, wd = flags_tools.command("veryfasttree", **{"in": trimmed})
          newick = self._run(vft, check=True, capture_output=True,
                             text=True, cwd=wd or None).stdout.strip()
          self.commands.append(" ".join(vft))
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
    import flags_tools
    cmd, wd = flags_tools.command(
        "iqtree", model="MFP", prefix=os.path.join(tmp, "iq"),
        threads=self.threads if self.threads else "AUTO", **{"in": aln})
    if shutil.which(cmd[0]) is None:
      alt = next((b for b in ("iqtree3", "iqtree2", "iqtree")
                  if shutil.which(b)), None)
      if alt is None:
        raise FileNotFoundError("no iqtree binary found")
      cmd[0] = alt
    if n_taxa >= 4:
      cmd += ["-B", "1000"]
    self._run(cmd, check=True, capture_output=True, text=True, cwd=wd or None)
    self.commands.append(" ".join(cmd))
    with open(os.path.join(tmp, "iq.treefile")) as fh:
      return fh.read().strip()

  @staticmethod
  def _read_alignment(path: str) -> Dict[str, str]:
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(path, "fasta")}

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