from pathlib import Path

GENES = [
	("A", 100, 400, "+", "WP_001", "protein one", "cds"),
	("A", 500, 900, "-", "WP_002", "protein two", "cds"),
	("A", 1000, 1080, "+", "rna_001", "tRNA-Ala", "tRNA"),
	("A", 1200, 1500, "+", None, "broken protein", "pseudo"),
	("A", 2000, 2600, "-", "WP_004", "query protein", "cds"),
	("A", 3000, 3500, "+", "WP_005", "protein five", "cds"),
	("A", 4000, 4800, "+", "WP_006", "protein six", "cds"),
	("A", 5000, 5200, "-", "WP_007", "protein seven", "cds"),
	("B", 100, 600, "+", "WP_008", "protein eight", "cds"),
	("B", 700, 1200, "-", "WP_009", "protein nine", "cds"),
]
LENGTHS = {"A": 20000, "B": 5000}
ORGANISM = "Escherichia coli K-12"


def write_genome(directory: Path, name: str = "GCF_TEST") -> None:
	directory.mkdir(parents=True, exist_ok=True)
	gff = ["##gff-version 3"]
	for contig, length in LENGTHS.items():
		gff.append("##sequence-region {} 1 {}".format(contig, length))
	faa, rna, genome = [], [], {}
	for contig, length in LENGTHS.items():
		genome[contig] = ("ACGT" * (length // 4 + 1))[:length]
	for n, (contig, start, end, strand, acc, product, kind) in enumerate(GENES, 1):
		locus = "TEST_{:04d}".format(n)
		gff.append("\t".join((contig, ".", "gene", str(start), str(end), ".", strand, ".",
			"ID=gene-{0};locus_tag={0};gene_biotype={1}".format(locus,
				"pseudogene" if kind == "pseudo" else ("tRNA" if kind == "tRNA" else "protein_coding")))))
		if kind == "cds":
			gff.append("\t".join((contig, ".", "CDS", str(start), str(end), ".", strand, "0",
				"ID=cds-{0};Parent=gene-{1};protein_id={0};locus_tag={1};product={2}".format(acc, locus, product))))
			faa.append(">{} {} [{}]\n{}".format(acc, product, ORGANISM, "M" + "A" * ((end - start + 1) // 3 - 1)))
		elif kind == "pseudo":
			gff.append("\t".join((contig, ".", "CDS", str(start), str(end), ".", strand, "0",
				"ID=cds-{0};Parent=gene-{0};locus_tag={0};product={1};pseudo=true".format(locus, product))))
		else:
			gff.append("\t".join((contig, ".", kind, str(start), str(end), ".", strand, ".",
				"ID=rna-{0};Parent=gene-{1};Name={0};locus_tag={1};product={2}".format(acc, locus, product))))
			rna.append(">{} [locus_tag={}]\n{}".format(acc, locus, "GC" * ((end - start + 1) // 2)))
	(directory / (name + "_genomic.gff")).write_text("\n".join(gff) + "\n")
	(directory / (name + "_protein.faa")).write_text("\n".join(faa) + "\n")
	(directory / (name + "_rna_from_genomic.fna")).write_text("\n".join(rna) + "\n")
	(directory / (name + "_genomic.fna")).write_text(
		"".join(">{}\n{}\n".format(c, s) for c, s in genome.items()))
