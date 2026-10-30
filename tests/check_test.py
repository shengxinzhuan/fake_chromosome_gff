#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""校验 fake_chrom_gff.py 的输出。用法: python check_test.py <outdir>"""

import sys
from collections import OrderedDict
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"


def parse_gff_genes(path):
    genes = OrderedDict()
    mrna = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "gene":
                genes[attr["ID"]] = (f[0], int(f[3]), int(f[4]), f[6])
            elif f[2] == "mRNA":
                mrna[attr["ID"]] = attr["Parent"]
    return genes, mrna


def fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def main():
    outdir = Path(sys.argv[1])
    ref_genes, _ = parse_gff_genes(DATA / "ref.gff3")
    fake_genes, fake_mrna = parse_gff_genes(outdir / "fake.gff3")

    rows = [l.rstrip("\n").split("\t")
            for l in open(outdir / "gene_mapping.tsv")][1:]
    by_new = {r[0]: r for r in rows}

    errors = []

    # 1. 基因总数: 23 = 19 同源 + 3 novel + 1 cds-only
    if len(fake_genes) != 23:
        errors.append(f"基因数 {len(fake_genes)} != 23")

    # 2. ortholog 行的坐标必须等于参考基因坐标; 染色体上按参考位置排序
    chrom_last_start = {}
    for r in rows:
        gid, tid, tgene, rgene, chrom, s, e, st, placement = r[:9]
        s, e = int(s), int(e)
        if gid not in fake_genes:
            errors.append(f"{gid} 不在 gff 中")
            continue
        if fake_genes[gid] != (chrom, s, e, st):
            errors.append(f"{gid} gff 坐标与 mapping 不一致")
        if fake_mrna.get(tid) != gid:
            errors.append(f"{tid} 的 Parent 不是 {gid}")
        if placement == "ortholog":
            if rgene not in ref_genes:
                errors.append(f"{gid} 的参考基因 {rgene} 不在参考 gff")
            elif ref_genes[rgene] != (chrom, s, e, st):
                errors.append(f"{gid} 未继承参考基因 {rgene} 的坐标")
            prev = chrom_last_start.setdefault(chrom, 0)
            if s < prev:
                errors.append(f"{chrom} 上基因未按参考位置排序: {s} < {prev}")
            chrom_last_start[chrom] = s
        else:
            if chrom != "ChrUn":
                errors.append(f"未匹配基因 {gid} 不在 ChrUn: {chrom}")

    # 3. 编号连续性: 每条染色体 G000010, G000020, ...
    seen = {}
    for gid in fake_genes:
        chrom, num = gid.split("G")[0], int(gid.split("G")[1])
        seen.setdefault(chrom, []).append(num)
    for chrom, nums in seen.items():
        if nums != [i * 10 for i in range(1, len(nums) + 1)]:
            errors.append(f"{chrom} 编号不连续: {nums}")

    # 4. fasta 头与 gff mRNA ID 一致
    cds_ids, pep_ids = fasta_ids(outdir / "fake.cds.fa"), fasta_ids(outdir / "fake.pep.fa")
    if set(cds_ids) != set(fake_mrna) - set(pep_ids) - (set(fake_mrna) - set(cds_ids)):
        pass
    if set(cds_ids) | set(pep_ids) != set(fake_mrna):
        errors.append("fasta 头与 gff mRNA ID 集合不一致")
    if len(cds_ids) != 23:
        errors.append(f"fake.cds.fa 应为 23 条, 实际 {len(cds_ids)}")
    if len(pep_ids) != 22:
        errors.append(f"fake.pep.fa 应为 22 条, 实际 {len(pep_ids)}")

    # 5. novel / cds-only 基因应在 ChrUn
    for r in rows:
        if r[2].startswith("rel_novel") or r[2] == "rel_cds_only":
            if r[4] != "ChrUn":
                errors.append(f"{r[2]} 应在 ChrUn, 实际 {r[4]}")

    # 6. gff 行数: 4 行/gene
    n_lines = sum(1 for l in open(outdir / "fake.gff3")
                  if l.strip() and not l.startswith("#"))
    if n_lines != len(fake_genes) * 4:
        errors.append(f"gff 特征行数 {n_lines} != {len(fake_genes) * 4}")

    if errors:
        print(f"[FAIL] {outdir}")
        for e in errors:
            print("  -", e)
        sys.exit(1)
    print(f"[PASS] {outdir}: {len(fake_genes)} genes, "
          f"{len(chrom_last_start) - (1 if 'ChrUn' in chrom_last_start else 0)} ref chroms + ChrUn")


if __name__ == "__main__":
    main()
