#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从 基因组fasta + GFF3 提取最长转录本的 CDS/PEP (低内存, SeqIO.index 随机访问)。

用法: python3 tools/extract_cds_pep.py genome.fa anno.gff3 out.cds.fa out.pep.fa
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from Bio.Seq import Seq

from fake_chrom_gff import parse_gff, representative_transcripts


def main():
    genome_fa, gff_path, cds_out, pep_out = sys.argv[1:5]
    genes, _tx, cds_by_tid, _order = parse_gff(gff_path)
    rep = representative_transcripts(genes, cds_by_tid)
    print(f"基因 {len(genes)}, 代表转录本 {len(rep)}", file=sys.stderr)
    idx = SeqIO.index(genome_fa, "fasta")
    n_skip = 0
    with open(cds_out, "w") as fc, open(pep_out, "w") as fp:
        for gid, tid in rep.items():
            g = genes[gid]
            if g["seqid"] not in idx:
                n_skip += 1
                continue
            chrom = idx[g["seqid"]].seq
            intervals = sorted(cds_by_tid[tid])
            cds = Seq("".join(str(chrom[s - 1:e]) for s, e in intervals))
            if g["strand"] == "-":
                cds = cds.reverse_complement()
            pep = cds[: len(cds) // 3 * 3].translate(table=1)
            fc.write(f">{tid}\n")
            s = str(cds)
            for i in range(0, len(s), 60):
                fc.write(s[i:i + 60] + "\n")
            fp.write(f">{tid}\n")
            s = str(pep).rstrip("*").replace("*", "X")
            for i in range(0, len(s), 60):
                fp.write(s[i:i + 60] + "\n")
    if n_skip:
        print(f"警告: {n_skip} 个转录本所在序列不在基因组中, 已跳过", file=sys.stderr)


if __name__ == "__main__":
    main()
