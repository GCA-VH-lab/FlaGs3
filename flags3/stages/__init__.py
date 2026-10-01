from flags3.stages.blast import Blast
from flags3.stages.cluster import Cluster, ClusterRna
from flags3.stages.defence import Defence
from flags3.stages.domains import Domains
from flags3.stages.extract import Extract
from flags3.stages.features import Features
from flags3.stages.fetch import Fetch
from flags3.stages.figures import Figures
from flags3.stages.genomad import Genomad
from flags3.stages.report import Report
from flags3.stages.sismis import Sismis
from flags3.stages.tree import Tree

PIPELINE = (Blast, Fetch, Extract, Cluster, ClusterRna, Tree, Domains, Features, Sismis, Genomad, Defence, Report, Figures)
BY_NAME = {stage.name: stage for stage in PIPELINE}
