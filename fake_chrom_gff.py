#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fake_chrom_gff.py — 为只有 CDS/PEP 的近缘种生成"假染色体" GFF。

流程:
  1. 从参考基因组取得基因位置与蛋白序列 (GFF+FA 提取, 或直接给 GFF+CDS+PEP)
  2. 近缘种 PEP vs 参考 PEP 做蛋白比对 (diamond 优先, 可 --aligner 指定), 每个基因取最佳命中 (可选 --rbh 互逆过滤)
  3. 匹配上的近缘种基因按参考同源基因的染色体/位置排序, 直接继承其坐标
  4. 未匹配上的基因放到 ChrUn 伪染色体上顺序排布
  5. 重新编号 (如 Chr01G000010), 输出 fake.gff3 / fake.cds.fa / fake.pep.fa / gene_mapping.tsv

仅供需要 GFF 坐标框架的工具使用 (共线性分析、绘图等), 坐标不代表真实位置。
"""

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict, defaultdict

from Bio import SeqIO
from Bio.Seq import Seq

ISOFORM_RE = re.compile(
    r"(\.mRNA\d+|\.[Tt]\d+|[-_]T\d+|-R[A-Z]|[-_.]seq\d+|\.p\d+|\.\d+)$")


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_fasta(path):
    """读 fasta; 重复 ID 时保留最长序列并警告 (真实数据常有 ID 碰撞)。"""
    seqs = OrderedDict()
    n_dup = 0
    with open_text(path) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            if rec.id in seqs:
                n_dup += 1
                if len(rec.seq) > len(seqs[rec.id]):
                    seqs[rec.id] = rec.seq
            else:
                seqs[rec.id] = rec.seq
    if n_dup:
        print(f"  警告: {path} 中有 {n_dup} 个重复序列 ID, 已保留最长的一条",
              file=sys.stderr)
    return seqs


def gene_id_of(tid):
    """剥离转录本/蛋白后缀得到基因 ID, 循环剥离以兼容叠加后缀
    (如 Trinity + TransDecoder 的 comp101_c0_seq1.p1 -> comp101_c0)。"""
    gid = tid
    for _ in range(5):
        new = ISOFORM_RE.sub("", gid)
        if new == gid or not new:
            break
        gid = new
    return gid or tid


def parse_gff(path):
    """解析 GFF3, 返回 (genes, transcripts, cds_by_tid, seqid_order)。

    genes: OrderedDict gid -> {seqid, start, end, strand, transcripts: [tid]}
    transcripts: tid -> gid
    cds_by_tid: tid -> [(start, end), ...]
    """
    genes_raw = OrderedDict()
    mrnas = []
    cds_by_tid = defaultdict(list)
    seqid_order = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seqid, _src, ftype, start, end, _score, strand, _phase, attrs = f
            if seqid not in seqid_order:
                seqid_order.append(seqid)
            attr = {}
            for kv in attrs.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr[k] = v
            start, end = int(start), int(end)
            if ftype == "gene":
                gid = attr.get("ID")
                if gid:
                    genes_raw[gid] = {
                        "seqid": seqid, "start": start, "end": end,
                        "strand": strand, "transcripts": [],
                    }
            elif ftype in ("mRNA", "transcript"):
                tid = attr.get("ID")
                if tid:
                    gid = attr.get("Parent", "").split(",")[0] or gene_id_of(tid)
                    mrnas.append((tid, gid, seqid, start, end, strand))
            elif ftype == "CDS":
                for parent in attr.get("Parent", "").split(","):
                    if parent:
                        cds_by_tid[parent].append((start, end))

    genes = OrderedDict()
    for gid, info in genes_raw.items():
        genes[gid] = dict(info)
    transcripts = {}
    for tid, gid, seqid, start, end, strand in mrnas:
        if tid in transcripts:
            continue
        if gid not in genes:
            genes[gid] = {"seqid": seqid, "start": start, "end": end,
                          "strand": strand, "transcripts": []}
        else:
            genes[gid]["start"] = min(genes[gid]["start"], start)
            genes[gid]["end"] = max(genes[gid]["end"], end)
        transcripts[tid] = gid
        genes[gid]["transcripts"].append(tid)
        if tid in cds_by_tid:
            genes[gid]["start"] = min(genes[gid]["start"],
                                      min(s for s, _e in cds_by_tid[tid]))
            genes[gid]["end"] = max(genes[gid]["end"],
                                    max(e for _s, e in cds_by_tid[tid]))
    # CDS 直接挂在 gene 上的情形
    for tid, intervals in cds_by_tid.items():
        if tid in transcripts:
            continue
        if tid in genes:
            transcripts[tid] = tid
            if tid not in genes[tid]["transcripts"]:
                genes[tid]["transcripts"].append(tid)
    return genes, transcripts, dict(cds_by_tid), seqid_order


def representative_transcripts(genes, cds_by_tid):
    """每个基因选 CDS 总长度最长的转录本。"""
    rep = {}
    for gid, g in genes.items():
        best, best_len = None, 0
        for tid in g["transcripts"]:
            clen = sum(e - s + 1 for s, e in cds_by_tid.get(tid, []))
            if clen > best_len:
                best, best_len = tid, clen
        if best is not None:
            rep[gid] = best
    return rep


def extract_from_genome(genome, genes, cds_by_tid, rep):
    """从基因组 fasta 按 GFF 提取代表转录本的 CDS 与蛋白序列。"""
    cds_out, pep_out = OrderedDict(), OrderedDict()
    skipped = []
    for gid, tid in rep.items():
        g = genes[gid]
        if g["seqid"] not in genome:
            skipped.append(gid)
            continue
        intervals = sorted(cds_by_tid[tid])
        seq = Seq("".join(str(genome[g["seqid"]][s - 1:e]) for s, e in intervals))
        if g["strand"] == "-":
            seq = seq.reverse_complement()
        cds_out[gid] = seq
        pep = seq[: len(seq) // 3 * 3].translate(table=1)
        pep_out[gid] = Seq(str(pep).rstrip("*").replace("*", "X"))
    if skipped:
        print(f"  警告: {len(skipped)} 个参考基因在基因组 fasta 中找不到序列, 已跳过",
              file=sys.stderr)
    return cds_out, pep_out


def assign_to_genes(seqs, genes, transcripts):
    """把 fasta 序列按 ID 归属到基因 (精确转录本 ID > 精确基因 ID > 去转录本后缀),
    同一基因取最长序列。返回 (gid -> seq, 未匹配 ID 列表)。"""
    out = OrderedDict()
    unmatched = []
    for sid, seq in seqs.items():
        gid = None
        if sid in genes:
            gid = sid
        elif sid in transcripts:
            gid = transcripts[sid]
        else:
            cand = gene_id_of(sid)
            if cand in genes:
                gid = cand
        if gid is None:
            unmatched.append(sid)
        elif gid not in out or len(seq) > len(out[gid]):
            out[gid] = seq
    return out, unmatched


def assign_target(seqs):
    """近缘种 fasta 全部保留, 转录本后缀合并为基因, 同基因取最长。"""
    out = OrderedDict()
    for sid, seq in seqs.items():
        gid = gene_id_of(sid)
        if gid not in out or len(seq) > len(out[gid]):
            out[gid] = seq
    return out


def emit_fasta(fh, rid, seq):
    fh.write(f">{rid}\n")
    s = str(seq)
    for i in range(0, len(s), 60):
        fh.write(s[i:i + 60] + "\n")


def write_fasta(path, records):
    with open(path, "w") as fh:
        for rid, seq in records:
            emit_fasta(fh, rid, seq)


def resolve_aligner(choice):
    """返回实际使用的比对器名称。choice: blastp / diamond / auto。"""
    if choice != "auto":
        if shutil.which(choice) is None:
            sys.exit(f"错误: 找不到 {choice}, 请先安装或改用 --aligner auto")
        return choice
    for cand in ("diamond", "blastp"):
        if shutil.which(cand):
            return cand
    sys.exit("错误: 找不到 diamond 或 blastp, 请先安装其一")


def run_blastp(query_fa, db_fa, tmpdir, tag, threads, evalue):
    db = os.path.join(tmpdir, f"db_{tag}")
    subprocess.run(["makeblastdb", "-in", db_fa, "-dbtype", "prot",
                    "-out", db, "-logfile", os.devnull], check=True)
    out = os.path.join(tmpdir, f"blast_{tag}.out")
    fmt = "6 qseqid sseqid pident length evalue bitscore"
    with open(out, "w") as fh:
        subprocess.run(["blastp", "-query", query_fa, "-db", db,
                        "-outfmt", fmt, "-evalue", str(evalue),
                        "-max_target_seqs", "10",
                        "-num_threads", str(threads)],
                       check=True, stdout=fh)
    return parse_tabular(out)


def run_diamond(query_fa, db_fa, tmpdir, tag, threads, evalue, sensitive=False):
    db = os.path.join(tmpdir, f"db_{tag}.dmnd")
    subprocess.run(["diamond", "makedb", "--in", db_fa, "-d", db,
                    "--quiet"], check=True)
    out = os.path.join(tmpdir, f"diamond_{tag}.out")
    cmd = ["diamond", "blastp", "--query", query_fa, "--db", db,
           "--outfmt", "6", "qseqid", "sseqid", "pident",
           "length", "evalue", "bitscore",
           "--evalue", str(evalue),
           "--max-target-seqs", "10",
           "--threads", str(threads), "--quiet"]
    if sensitive:
        cmd.append("--more-sensitive")
    with open(out, "w") as fh:
        subprocess.run(cmd, check=True, stdout=fh)
    return parse_tabular(out)


def parse_tabular(out):
    hits = []
    with open(out) as fh:
        for line in fh:
            q, s, pident, length, ev, bs = line.split()
            hits.append((q, s, float(pident), int(length), float(ev), float(bs)))
    return hits


def best_hit_per_query(hits):
    """每个 query 取 (bitscore, pident, length) 最高的命中。"""
    best = {}
    for q, s, pident, length, ev, bs in hits:
        score = (bs, pident, length)
        if q not in best or score > best[q][0]:
            best[q] = (score, s, ev)
    return {q: {"subject": v[1], "bitscore": v[0][0], "pident": v[0][1],
                "evalue": v[2]} for q, v in best.items()}


def build_placements(target_genes, best, ref_genes, chrun, gap):
    """返回 (matched, unplaced)。
    matched: chrom -> [(target_gid, ref_gid, hit)], 按 ref 基因位置排序。
    unplaced: [target_gid] 保持输入顺序。"""
    matched = defaultdict(list)
    unplaced = []
    for tg in target_genes:
        hit = best.get(tg)
        if hit is None:
            unplaced.append(tg)
            continue
        rg = hit["subject"]
        if rg not in ref_genes:
            unplaced.append(tg)
            continue
        matched[ref_genes[rg]["seqid"]].append((tg, rg, hit))
    for chrom in matched:
        matched[chrom].sort(key=lambda x: (ref_genes[x[1]]["start"], x[0]))
    return matched, unplaced


def main():
    ap = argparse.ArgumentParser(
        description="根据参考基因组的基因顺序, 为近缘种 CDS/PEP 生成假染色体 GFF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--ref-gff", required=True, help="参考基因组 GFF3 (染色体级别)")
    ap.add_argument("--ref-fa", help="参考基因组 fasta (与 --ref-cds/--ref-pep 二选一)")
    ap.add_argument("--ref-cds", help="参考 CDS fasta (需同时给 --ref-pep)")
    ap.add_argument("--ref-pep", help="参考 PEP fasta (需同时给 --ref-cds)")
    ap.add_argument("--cds", required=True, help="近缘种 CDS fasta")
    ap.add_argument("--pep", required=True, help="近缘种 PEP fasta")
    ap.add_argument("--outdir", required=True, help="输出目录")
    ap.add_argument("--evalue", type=float, default=1e-5, help="蛋白比对 e-value 阈值")
    ap.add_argument("--aligner", choices=["auto", "blastp", "diamond"],
                    default="auto",
                    help="蛋白比对器; auto 优先使用 diamond (快), 其次 blastp")
    ap.add_argument("--sensitive", action="store_true",
                    help="diamond 使用 --more-sensitive (更敏感稍慢, 仅 diamond 有效)")
    ap.add_argument("--min-pident", type=float, default=0.0,
                    help="最佳命中 pident 低于该值的基因按未匹配处理 (0 = 不过滤)")
    ap.add_argument("--shift-duplicates", action="store_true",
                    help="共享参考坐标的基因向右平移, 保证任意基因不重叠 "
                         "(坐标不再忠实参考, 推荐 MCScanX 等工具使用)")
    ap.add_argument("--rbh", action="store_true", help="仅保留 reciprocal best hit")
    ap.add_argument("--threads", type=int, default=4, help="比对线程数")
    ap.add_argument("--tag", default="G", help="基因编号标签, 如 Chr01G000010")
    ap.add_argument("--chrun", default="ChrUn", help="未匹配基因所在伪染色体名")
    ap.add_argument("--gap", type=int, default=1000, help="ChrUn 上基因间间隔 (bp)")
    ap.add_argument("--pad", type=int, default=10000, help="染色体末端 padding (bp)")
    ap.add_argument("--make-genome", action="store_true",
                    help="额外输出由 N 组成的 fake.genome.fa")
    args = ap.parse_args()

    if bool(args.ref_fa) == bool(args.ref_cds or args.ref_pep):
        sys.exit("错误: 请二选一: --ref-fa, 或 --ref-cds + --ref-pep")
    if (args.ref_cds or args.ref_pep) and not (args.ref_cds and args.ref_pep):
        sys.exit("错误: --ref-cds 与 --ref-pep 需要同时提供")

    os.makedirs(args.outdir, exist_ok=True)

    print("[1/4] 解析参考 GFF ...", file=sys.stderr)
    ref_genes, ref_tx, ref_cds_by_tid, seqid_order = parse_gff(args.ref_gff)
    if not ref_genes:
        sys.exit(f"错误: {args.ref_gff} 中没有可用的 gene 记录")
    print(f"      参考基因: {len(ref_genes)} 条, 参考染色体/序列: {len(seqid_order)} 条",
          file=sys.stderr)

    if args.ref_fa:
        print("[2/4] 从参考基因组提取 CDS/PEP ...", file=sys.stderr)
        genome = load_fasta(args.ref_fa)
        rep = representative_transcripts(ref_genes, ref_cds_by_tid)
        if not rep:
            sys.exit("错误: GFF 中没有任何 CDS 记录, 无法从基因组提取序列; "
                     "请改用 --ref-cds/--ref-pep 模式")
        ref_cds, ref_pep = extract_from_genome(genome, ref_genes, ref_cds_by_tid, rep)
    else:
        print("[2/4] 读取参考 CDS/PEP ...", file=sys.stderr)
        ref_cds_raw = load_fasta(args.ref_cds)
        ref_pep_raw = load_fasta(args.ref_pep)
        ref_cds, un1 = assign_to_genes(ref_cds_raw, ref_genes, ref_tx)
        ref_pep, un2 = assign_to_genes(ref_pep_raw, ref_genes, ref_tx)
        if un1 or un2:
            print(f"      警告: {len(un1) + len(un2)} 条参考序列无法对应到 GFF 基因, 已忽略",
                  file=sys.stderr)
        ref_cds = OrderedDict((g, s) for g, s in ref_cds.items() if g in ref_pep)
        ref_pep = OrderedDict((g, s) for g, s in ref_pep.items() if g in ref_cds)
    ref_pep = OrderedDict((g, s) for g, s in ref_pep.items() if g in ref_genes)
    if not ref_pep:
        sys.exit("错误: 没有任何参考蛋白序列与 GFF 基因对应")
    print(f"      可用于比对的参考蛋白: {len(ref_pep)} 条", file=sys.stderr)

    print("[3/4] 读取近缘种 CDS/PEP 并构建代表序列 ...", file=sys.stderr)
    tgt_cds_raw = load_fasta(args.cds)
    tgt_pep_raw = load_fasta(args.pep)
    n_tx = max(len(tgt_cds_raw), len(tgt_pep_raw))
    tgt_cds = assign_target(tgt_cds_raw)
    tgt_pep = assign_target(tgt_pep_raw)
    tgt_genes = list(OrderedDict.fromkeys(list(tgt_cds) + list(tgt_pep)))
    n_cds_only = sum(1 for g in tgt_genes if g not in tgt_pep)
    n_pep_only = sum(1 for g in tgt_genes if g not in tgt_cds)
    if n_cds_only or n_pep_only:
        print(f"      警告: {n_cds_only} 个基因仅有 CDS, {n_pep_only} 个基因仅有 PEP "
              f"(仅 CDS 的基因将无法参与比对, 放到 {args.chrun})", file=sys.stderr)
    print(f"      近缘种基因: {len(tgt_genes)} 个 (由 {n_tx} 条转录本合并), "
          f"其中 {len(tgt_pep)} 个有蛋白可用于比对", file=sys.stderr)

    with tempfile.TemporaryDirectory() as tmpdir:
        tgt_pep_fa = os.path.join(tmpdir, "tgt.pep.fa")
        ref_pep_fa = os.path.join(tmpdir, "ref.pep.fa")
        write_fasta(tgt_pep_fa, tgt_pep.items())
        write_fasta(ref_pep_fa, ref_pep.items())
        aligner = resolve_aligner(args.aligner)

        def run_one(query_fa, db_fa, tag):
            if aligner == "diamond":
                return run_diamond(query_fa, db_fa, tmpdir, tag,
                                   args.threads, args.evalue, args.sensitive)
            return run_blastp(query_fa, db_fa, tmpdir, tag,
                              args.threads, args.evalue)

        print(f"[4/4] {aligner} 比对 (近缘种 -> 参考) ...", file=sys.stderr)
        fwd = best_hit_per_query(
            run_one(tgt_pep_fa, ref_pep_fa, "fwd"))
        best = {tg: {"subject": h["subject"], "bitscore": h["bitscore"],
                     "pident": h["pident"], "evalue": h["evalue"]}
                for tg, h in fwd.items()}
        if args.rbh:
            print(f"      {aligner} 比对 (参考 -> 近缘种, RBH 反向) ...", file=sys.stderr)
            rev = best_hit_per_query(
                run_one(ref_pep_fa, tgt_pep_fa, "rev"))
            rev_best = {rg: h["subject"] for rg, h in rev.items()}
            before = len(best)
            best = {tg: h for tg, h in best.items()
                    if rev_best.get(h["subject"]) == tg}
            print(f"      RBH 过滤: {before} -> {len(best)}", file=sys.stderr)
        if args.min_pident > 0:
            before = len(best)
            best = {tg: h for tg, h in best.items()
                    if h["pident"] >= args.min_pident}
            print(f"      pident >= {args.min_pident} 过滤: "
                  f"{before} -> {len(best)} (被过滤的进入 {args.chrun})",
                  file=sys.stderr)

    matched, unplaced = build_placements(tgt_genes, best, ref_genes, args.chrun, args.gap)

    print("编号并输出 ...", file=sys.stderr)
    chroms_used = [c for c in seqid_order if c in matched]
    mapping_rows = []
    placed = []
    counters = defaultdict(int)
    n_shifted = 0
    for chrom in chroms_used:
        cursor = 0
        for tg, rg, hit in matched[chrom]:
            counters[chrom] += 1
            new_gid = f"{chrom}{args.tag}{counters[chrom] * 10:06d}"
            g = ref_genes[rg]
            s, e = g["start"], g["end"]
            if args.shift_duplicates:
                if s <= cursor:
                    s = cursor + args.gap + 1
                    e = s + (g["end"] - g["start"])
                    n_shifted += 1
                cursor = e
            placed.append({
                "new_gene": new_gid, "target": tg, "ref": rg,
                "chrom": chrom, "start": s, "end": e,
                "strand": g["strand"], "hit": hit,
            })
    pos = 1
    for tg in unplaced:
        counters[args.chrun] += 1
        new_gid = f"{args.chrun}{args.tag}{counters[args.chrun] * 10:06d}"
        glen = len(tgt_cds.get(tg) or tgt_pep[tg] * 3)
        start, end = pos, pos + max(glen, 300) + 199
        placed.append({
            "new_gene": new_gid, "target": tg, "ref": None,
            "chrom": args.chrun, "start": start, "end": end,
            "strand": "+", "hit": None,
        })
        pos = end + 1 + args.gap

    # 染色体长度: 参考染色体按基因组 fasta, 否则按基因末端 + padding
    chrom_len = OrderedDict()
    if args.ref_fa:
        for c in seqid_order:
            chrom_len[c] = len(genome[c])
    else:
        for c in seqid_order:
            ends = [ref_genes[g]["end"] for g in ref_genes
                    if ref_genes[g]["seqid"] == c]
            chrom_len[c] = (max(ends) + args.pad) if ends else args.pad
    chrom_len[args.chrun] = pos + args.pad if unplaced else 0
    if args.shift_duplicates:
        for p in placed:
            if p["end"] >= chrom_len[p["chrom"]]:
                chrom_len[p["chrom"]] = p["end"] + 1

    gff_path = os.path.join(args.outdir, "fake.gff3")
    with open(gff_path, "w") as fh:
        fh.write("##gff-version 3\n")
        for c in seqid_order:
            fh.write(f"##sequence-region {c} 1 {chrom_len[c]}\n")
        if unplaced:
            fh.write(f"##sequence-region {args.chrun} 1 {chrom_len[args.chrun]}\n")
        for p in placed:
            gid, tid = p["new_gene"], p["new_gene"] + ".t1"
            row = [p["chrom"], "fakeGFF", "gene", p["start"], p["end"],
                   ".", p["strand"], ".", f"ID={gid};Name={gid}"]
            fh.write("\t".join(map(str, row)) + "\n")
            row[2], row[8] = "mRNA", f"ID={tid};Parent={gid}"
            fh.write("\t".join(map(str, row)) + "\n")
            row[2], row[8] = "exon", f"ID={tid}.exon1;Parent={tid}"
            fh.write("\t".join(map(str, row)) + "\n")
            row[2], row[7], row[8] = "CDS", "0", f"ID={tid}.cds1;Parent={tid}"
            fh.write("\t".join(map(str, row)) + "\n")

    cds_path = os.path.join(args.outdir, "fake.cds.fa")
    pep_path = os.path.join(args.outdir, "fake.pep.fa")
    n_cds = n_pep = 0
    with open(cds_path, "w") as fc, open(pep_path, "w") as fp:
        for p in placed:
            tid = p["new_gene"] + ".t1"
            tg = p["target"]
            if tg in tgt_cds:
                emit_fasta(fc, f"{tid} {tg}", tgt_cds[tg])
                n_cds += 1
            if tg in tgt_pep:
                emit_fasta(fp, f"{tid} {tg}", tgt_pep[tg])
                n_pep += 1

    map_path = os.path.join(args.outdir, "gene_mapping.tsv")
    with open(map_path, "w") as fh:
        fh.write("new_gene_id\tnew_mrna_id\ttarget_gene_id\tref_gene_id\t"
                 "chrom\tstart\tend\tstrand\tplacement\tpident\tevalue\tbitscore\n")
        for p in placed:
            tid = p["new_gene"] + ".t1"
            if p["hit"] is not None:
                fh.write(f"{p['new_gene']}\t{tid}\t{p['target']}\t{p['ref']}\t"
                         f"{p['chrom']}\t{p['start']}\t{p['end']}\t{p['strand']}\t"
                         f"ortholog\t{p['hit']['pident']:.2f}\t{p['hit']['evalue']:.3g}"
                         f"\t{p['hit']['bitscore']:.1f}\n")
            else:
                fh.write(f"{p['new_gene']}\t{tid}\t{p['target']}\t.\t"
                         f"{p['chrom']}\t{p['start']}\t{p['end']}\t{p['strand']}\t"
                         f"unplaced\t.\t.\t.\n")

    if args.make_genome:
        genome_path = os.path.join(args.outdir, "fake.genome.fa")
        with open(genome_path, "w") as fh:
            for c, length in chrom_len.items():
                if length == 0:
                    continue
                fh.write(f">{c}\n")
                for i in range(0, length, 80):
                    fh.write("N" * min(80, length - i) + "\n")

    n_matched = len(placed) - len(unplaced)
    pids = [p["hit"]["pident"] for p in placed if p["hit"] is not None]
    pident_stat = ""
    if pids:
        b = [sum(1 for x in pids if x >= 90),
             sum(1 for x in pids if 70 <= x < 90),
             sum(1 for x in pids if 50 <= x < 70),
             sum(1 for x in pids if x < 50)]
        pident_stat = (f"\n  命中一致性: mean={sum(pids) / len(pids):.1f}%  "
                       f">=90%:{b[0]}  70-90%:{b[1]}  50-70%:{b[2]}  <50%:{b[3]}")
    dup = len(matched and [g for c in matched for _tg, g, _h in matched[c]]) \
        - len({g for c in matched for _tg, g, _h in matched[c]})
    print(f"""
完成! 输出目录: {args.outdir}
  fake.gff3        假染色体 GFF (gene/mRNA/exon/CDS)
  fake.cds.fa      按假染色体顺序重排的近缘种 CDS ({n_cds} 条)
  fake.pep.fa      按假染色体顺序重排的近缘种 PEP ({n_pep} 条)
  gene_mapping.tsv 新 ID -> 原 ID / 参考基因 对照表
""" + ("  fake.genome.fa   N 填充的假基因组序列\n" if args.make_genome else "")
          + f"""统计:
  匹配到参考基因(继承其坐标): {n_matched}
  未匹配 -> {args.chrun}: {len(unplaced)}
  多个近缘种基因命中同一参考基因(坐标重复): {max(dup, 0)}
  每染色体匹配基因数: """ + ", ".join(
              f"{c}:{len(matched[c])}" for c in chroms_used) + pident_stat
          + (f"\n  坐标平移(消除重叠): {n_shifted}" if args.shift_duplicates else ""),
          file=sys.stderr)


if __name__ == "__main__":
    main()
