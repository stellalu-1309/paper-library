# -*- coding: utf-8 -*-
"""
從論文 PDF 拆出參考文獻清單，並統計每一條在內文被引用幾次、出現在哪幾頁。

需要 pdftotext（poppler-utils 或 xpdf）。

用法：
    python tools/extract_refs.py paper.pdf
    python tools/extract_refs.py paper.pdf --json out.json

設計原則（和整個專案一致）：拆不出來的就留白，不猜、不補。
每一筆都帶頁碼，方便回原文核對。
"""

import argparse
import collections
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

# 參考文獻清單的起始標題，依序嘗試
REF_HEADINGS = [
    "References and Notes",
    "\nReferences\n",
    "\nREFERENCES\n",
    "\nLiterature Cited\n",
]

# 內文裡看起來像引用、其實不是的前綴（腦區、病毒、試劑、單位…）
# 用來擋掉 CA1、AAV9、GCaMP6、5xFAD 這類命名中的數字
NOT_A_CITATION = re.compile(
    r"(CA|Fig|Figure|Table|Extended|Data|Movie|AAV|hM|GCaMP|ArchT|TetO|ChR|"
    r"ZT|NS|SD|RS|PRE|POST|Hz|cm|mm|ms|mV|nm|kHz|v|p|n)$",
    re.I,
)


def pdf_to_pages(pdf_path):
    """
    用 pdftotext 抽文字，回傳 {頁碼: 該頁文字}。頁碼從 1 開始，對應 PDF 實際頁數。

    注意：Windows 上的 xpdf pdftotext 打不開含非 ASCII 字元的路徑
    （例如資料夾名稱有中文時，會直接回 I/O Error）。
    所以先把 PDF 複製到純 ASCII 的暫存路徑再處理。
    """
    workdir = tempfile.mkdtemp(prefix="refx_")
    src = os.path.join(workdir, "in.pdf")
    tmp = os.path.join(workdir, "out.txt")
    try:
        shutil.copyfile(pdf_path, src)
    except OSError as e:
        sys.exit("讀不到 PDF：%s（%s）" % (pdf_path, e))
    pdf_path, _orig = src, pdf_path
    try:
        # 掃描版或壓縮流有損的 PDF 常讓 pdftotext 回傳非零，但文字其實抽得出來。
        # 所以不看 exit code，只看有沒有抽到東西。
        subprocess.run(
            ["pdftotext", "-enc", "UTF-8", pdf_path, tmp],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        text = io.open(tmp, encoding="utf-8", errors="replace").read()
        if not text.strip():
            sys.exit("這份 PDF 抽不到文字，可能是掃描檔，需要先做 OCR：%s" % _orig)
    except FileNotFoundError:
        sys.exit("找不到 pdftotext。請先安裝 poppler-utils 或 xpdf。")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    # pdftotext 用 \f 分頁
    return {i + 1: t for i, t in enumerate(text.split("\f")) if t.strip()}


def split_reference_list(full_text):
    """切出參考文獻區塊。找不到就回 None —— 不去猜它在哪。"""
    for h in REF_HEADINGS:
        if h in full_text:
            return full_text.split(h, 1)[1]
    return None


def parse_refs(blob):
    """把 '1. Author ... (2019). [PubMed: 123]' 逐條拆開，回傳 {編號: 原文字串}。"""
    # 在「數字 + 點 + 空白 + 大寫字母」處切開
    parts = re.split(r"(?:(?<=\])|(?<=\.)|^)\s*(\d{1,3})\.\s+(?=[A-Z])", blob)
    refs = {}
    i = 1
    while i < len(parts) - 1:
        n = int(parts[i])
        body = re.sub(r"\s+", " ", parts[i + 1]).strip()
        if n not in refs:          # 同一編號只留第一次出現的
            refs[n] = body
        i += 2
    return refs


def ref_fields(raw):
    pmid = re.search(r"PubMed:\s*(\d+)", raw)
    doi = re.search(r"\b(10\.\d{4,9}/[^\s,;]+)", raw)
    year = re.search(r"\((\d{4})\)", raw)
    return {
        "raw": raw,
        "pmid": pmid.group(1) if pmid else None,
        "doi": doi.group(1).rstrip(".") if doi else None,
        "year": year.group(1) if year else None,
    }


def expand(group, max_refs, max_range=12):
    """把 '12-16' 或 '1,3' 展開成編號清單。過長的區間視為誤判，丟掉。"""
    out = []
    for chunk in group.split(","):
        chunk = chunk.strip()
        rng = re.match(r"^(\d+)\s*[–—-]\s*(\d+)$", chunk)
        if rng:
            a, b = int(rng.group(1)), int(rng.group(2))
            if 0 < b - a < max_range:
                out += list(range(a, b + 1))
        elif chunk.isdigit():
            out.append(int(chunk))
    return [n for n in out if 1 <= n <= max_refs]


def count_parenthetical(pages, refs, body_pages):
    """括號式引用，例如 (17)、(12-16)。這種格式可靠。"""
    cites = collections.defaultdict(list)
    for p in body_pages:
        t = pages.get(p, "")
        for m in re.finditer(r"\((\d{1,3}(?:\s*[,–—-]\s*\d{1,3})*)\)", t):
            for n in expand(m.group(1), max(refs)):
                if n in refs:
                    ctx = re.sub(r"\s+", " ", t[max(0, m.start() - 90):m.start()]).strip()
                    cites[n].append({"page": p, "context": ctx[-85:]})
    return cites


def count_superscript(pages, refs, body_pages):
    """
    上標式引用，例如 memory^5,6。

    上標抽成純文字後會和統計數字混在一起（'0.57 Hz' 看起來像 '引用第 57 條'），
    所以規則收得很緊：數字前面必須緊接小寫字母，且前 5 個字元內不能有數字、
    等號或乘號。寧可漏抓，不要錯抓。

    真正可靠的做法是讀 PDF 的字型大小與基線（PyMuPDF page.get_text("dict")），
    上標字級較小、基線較高。這個函式只是沒有該資訊時的退路。
    """
    cites = collections.defaultdict(list)
    for p in body_pages:
        t = pages.get(p, "")
        for m in re.finditer(
            r"(?<=[a-z])(\d{1,3}(?:[,–—-]\d{1,3})*)(?![\d.–—)%a-zA-Z-])", t
        ):
            pre = t[max(0, m.start() - 5):m.start()]
            if re.search(r"[\d=×]", pre):
                continue
            word = re.search(r"([A-Za-z]+)$", t[max(0, m.start() - 14):m.start()])
            if word and NOT_A_CITATION.search(word.group(1)):
                continue
            for n in expand(m.group(1), max(refs)):
                if n in refs:
                    ctx = re.sub(r"\s+", " ", t[max(0, m.start() - 90):m.start()]).strip()
                    cites[n].append({"page": p, "context": ctx[-85:]})
    return cites


def detect_style(pages, refs, body_pages):
    """看哪一種格式命中的多，就當作是哪一種。"""
    par = count_parenthetical(pages, refs, body_pages)
    sup = count_superscript(pages, refs, body_pages)
    n_par = sum(len(v) for v in par.values())
    n_sup = sum(len(v) for v in sup.values())
    if n_par >= n_sup:
        return "parenthetical", par, True
    return "superscript", sup, False


def main():
    ap = argparse.ArgumentParser(description="拆解論文參考文獻並統計引用次數")
    ap.add_argument("pdf", help="論文 PDF")
    ap.add_argument("--json", help="把完整結果寫成 JSON")
    ap.add_argument("--top", type=int, default=10, help="列出前幾名（預設 10）")
    args = ap.parse_args()

    pages = pdf_to_pages(args.pdf)
    full = "\f".join(pages.get(i, "") for i in range(1, max(pages) + 1))

    blob = split_reference_list(full)
    if blob is None:
        sys.exit("找不到參考文獻區塊的標題。這份 PDF 可能用了別的寫法，請自己看一下。")

    refs = {n: ref_fields(r) for n, r in parse_refs(blob).items()}
    if not refs:
        sys.exit("參考文獻區塊找到了，但一條都拆不出來。")

    # 參考文獻清單開始的那一頁之前，才算內文
    ref_start_page = next(
        (p for p in sorted(pages) if any(h.strip() in pages[p] for h in REF_HEADINGS)),
        max(pages),
    )
    body_pages = [p for p in sorted(pages) if p < ref_start_page]

    style, cites, reliable = detect_style(pages, refs, body_pages)
    total = max(refs)
    missing = [i for i in range(1, total + 1) if i not in refs]
    with_pmid = sum(1 for r in refs.values() if r["pmid"])

    print("檔案：%s" % os.path.basename(args.pdf))
    print("總頁數：%d　內文頁：1-%d　參考文獻從第 %d 頁開始" % (max(pages), len(body_pages), ref_start_page))
    print()
    print("參考文獻編號到：%d" % total)
    print("成功拆出：     %d" % len(refs))
    print("有 PMID：      %d" % with_pmid)
    if missing:
        print("拆不出來：     第 %s 條" % "、".join(map(str, missing)))
        print("               （通常是跨頁被切斷。讀到多少留多少，沒有補。）")
    print()
    print("引用格式：%s" % ("括號式，例如 (17)" if style == "parenthetical" else "上標式，例如 memory^5,6"))
    print("引用位置：%s" % ("可靠" if reliable else "僅供參考——上標會和統計數字混淆，這裡是保守估計，一定有漏"))
    print("共抓到 %d 處引用標記" % sum(len(v) for v in cites.values()))
    print()

    if cites:
        print("被引用最多次的 %d 條：" % min(args.top, len(cites)))
        for n, lst in sorted(cites.items(), key=lambda x: -len(x[1]))[:args.top]:
            pgs = sorted({c["page"] for c in lst})
            pm = refs[n]["pmid"]
            print("  [%3d] %2d 次　p.%s%s" % (n, len(lst), "、p.".join(map(str, pgs)),
                                              "　PMID %s" % pm if pm else ""))
            print("        %s" % refs[n]["raw"][:100])
            if lst:
                print("        例：…%s" % lst[0]["context"][-55:])
    else:
        print("內文裡沒抓到任何引用標記。")

    if args.json:
        data = {
            "file": os.path.basename(args.pdf),
            "style": style,
            "reliable": reliable,
            "total": total,
            "parsed": len(refs),
            "missing": missing,
            "with_pmid": with_pmid,
            "refs": refs,
            "citations": {str(k): v for k, v in cites.items()},
        }
        with io.open(args.json, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=1)
        print("\n已寫入 %s" % args.json)


if __name__ == "__main__":
    main()
