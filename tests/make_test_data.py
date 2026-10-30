#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成模拟测试数据: 染色体级别的参考基因组 + 近缘种 CDS/PEP。

参考基因组: 2 条染色体, 24 个基因 (单外显子, 随机正负链)。
近缘种: 取 19 个参考基因的同源拷贝 (随机碱基突变), 加 3 个全新随机基因,
        加 1 个仅 CDS 的基因; 转录本 ID 混用 .t1 / -T1 / -RA 后缀,
        记录顺序打乱。其中 4 个参考基因在近缘种中缺失。
"""

import random
import sys
from collections import OrderedDict
from pathlib import Path

from Bio.Seq import Seq

random.seed(42)
OUT = Path(__file__).resolve().parent / "data"
OUT.mkdir(parents=True, exist_ok=True)

COMP = str.maketrans("ACGT", "TGCA")


def rand_seq(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def mutate(seq, rate):
    out = []
    i = 0
    while i < len(seq):
        r = random.random()
        if r < rate:
            out.append(random.choice([b for b in "ACGT" if b != seq[i]]))
            i += 1
        elif r < rate + 0.005 and i + 300 < len(seq):
            i += 3  # 缺失一个密码子
        else:
            out.append(seq[i])
            i += 1
    return "".join(out)


def wrap(seq, fh, header):
    fh.write(f">{header}\n")
    for i in range(0, len(seq), 60):
        fh.write(seq[i:i + 60] + "\n")


# ---------------- 参考基因组 ----------------
chroms = OrderedDict([("chrA", 100_000), ("chrB", 80_000)])
genome = {c: rand_seq(n) for c, n in chroms.items()}

genes = []  # (gid, chrom, start, end, strand, cds_start, cds_end)
for chrom, clen in chroms.items():
    pos = 3_000
    gi = 0
    while pos + 4_000 < clen - 3_000:
        gi += 1
        gid = f"{chrom}g{gi:02d}"
        glen = random.randint(1_200, 2_800)
        start, end = pos, pos + glen
        strand = random.choice(["+", "-"])
        cs, ce = start + 120, end - 120
        ce = cs - 1 + (ce - cs + 1) // 3 * 3  # GFF 闭区间, 长度需为 3 的倍数
        genes.append((gid, chrom, start, end, strand, cs, ce))
        pos = end + random.randint(2_000, 6_000)

with open(OUT / "ref.gff3", "w") as fh:
    fh.write("##gff-version 3\n")
    for c, n in chroms.items():
        fh.write(f"##sequence-region {c} 1 {n}\n")
    for gid, chrom, start, end, strand, cs, ce in genes:
        tid = gid + (".t1" if random.random() < 0.5 else "-T1")
        fh.write(f"{chrom}\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid}\n")
        fh.write(f"{chrom}\ttest\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t"
                 f"ID={tid};Parent={gid}\n")
        fh.write(f"{chrom}\ttest\texon\t{cs}\t{ce}\t.\t{strand}\t.\t"
                 f"ID={tid}.exon1;Parent={tid}\n")
        fh.write(f"{chrom}\ttest\tCDS\t{cs}\t{ce}\t.\t{strand}\t0\t"
                 f"ID={tid}.cds1;Parent={tid}\n")

with open(OUT / "ref.genome.fa", "w") as fh:
    for c, s in genome.items():
        wrap(s, fh, c)

ref_cds, ref_pep = OrderedDict(), OrderedDict()
for gid, chrom, _s, _e, strand, cs, ce in genes:
    cds = genome[chrom][cs - 1:ce]
    if strand == "-":
        cds = cds.translate(COMP)[::-1]
    ref_cds[gid] = cds
    ref_pep[gid] = str(Seq(cds).translate())

with open(OUT / "ref.cds.fa", "w") as fh:
    for gid, chrom, _s, _e, _st, _cs, _ce in genes:
        tid = gid + (".t1" if gid[-1] in "13579" else "-T1")
        wrap(ref_cds[gid], fh, tid)
with open(OUT / "ref.pep.fa", "w") as fh:
    for gid, chrom, _s, _e, _st, _cs, _ce in genes:
        tid = gid + (".t1" if gid[-1] in "13579" else "-T1")
        wrap(ref_pep[gid], fh, tid)

# ---------------- 近缘种 ----------------
kept = random.sample(genes, 19)          # 19 个有同源拷贝
missing = [g for g in genes if g not in kept]  # 27-19=8 个参考基因在近缘种中缺失

suffixes = [".1", "-T1", "-RA"]
tgt_cds, tgt_pep = OrderedDict(), OrderedDict()

for i, (gid, chrom, s, e, st, cs, ce) in enumerate(kept):
    tgid = f"rel_{gid}"
    tid = tgid + suffixes[i % len(suffixes)]
    cds = mutate(ref_cds[gid], 0.04)
    cds = cds[: len(cds) // 3 * 3]
    pep = str(Seq(cds).translate()).rstrip("*").replace("*", "X")
    tgt_cds[tid] = cds
    tgt_pep[tid] = pep

for j in range(3):                          # 3 个全新基因 (应进 ChrUn)
    tid = f"rel_novel{j + 1}.1"
    cds = rand_seq(810)
    tgt_cds[tid] = cds
    tgt_pep[tid] = str(Seq(cds).translate())

cds_only_tid = "rel_cds_only-T1"            # 仅 CDS, 无蛋白 -> ChrUn
tgt_cds[cds_only_tid] = mutate(ref_cds[kept[0][0]], 0.04)[:600]

items = list(tgt_cds.items())
random.shuffle(items)
with open(OUT / "rel.cds.fa", "w") as fh:
    for tid, s in items:
        wrap(s, fh, tid)
items = list(tgt_pep.items())
random.shuffle(items)
with open(OUT / "rel.pep.fa", "w") as fh:
    for tid, s in items:
        wrap(s, fh, tid)

print(f"测试数据已写入 {OUT}", file=sys.stderr)
print(f"  参考基因: {len(genes)} (chrA: {sum(1 for g in genes if g[1] == 'chrA')}, "
      f"chrB: {sum(1 for g in genes if g[1] == 'chrB')})", file=sys.stderr)
print(f"  近缘种转录本: {len(tgt_cds)} CDS / {len(tgt_pep)} PEP", file=sys.stderr)
