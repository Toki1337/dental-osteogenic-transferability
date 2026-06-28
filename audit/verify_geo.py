#!/usr/bin/env python3
"""Verify GEO accessions via NCBI E-utilities (esearch -> esummary, db=gds).
Honest data audit: existence, n_samples, platform, organism, submission date.
"""
import json, time, sys, urllib.request, urllib.parse

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# (accession, planned_role, planned_note)
ACCS = [
    # success / dental osteogenic (RRA pool)
    ("GSE99958",  "success/dental", "PDLSC osteo RNA-seq+miRNA"),
    ("GSE163354", "success/dental", "PDLSC SDF-1/Exendin-4"),
    ("GSE226347", "success/dental", "DPSC-ECM Wnt/Hippo"),
    ("GSE296018", "success/dental", "DPSC cuscuta flavonoid"),
    ("GSE271641", "success/dental", "DPSC CBD odonto/osteo"),
    ("GSE299041", "success/dental", "DPSC shear/p38"),
    ("GSE286540", "success/dental", "SCAP simvastatin"),
    ("GSE236009", "success/dental", "PDLSC Cald1"),
    ("GSE266150", "success/dental", "PDLSC naringenin"),
    ("GSE159507", "success/dental", "PDLSC osteo mRNA/lncRNA"),
    ("GSE159508", "success/dental", "PDLSC osteo miRNA"),
    ("GSE266257", "success/dental", "SHED inorganic phosphate"),
    ("GSE49007",  "success/dental", "DFC dental follicle D7 array"),
    ("GSE316449", "success/dental-ext", "PDLF a-KG RNA-seq (ext valid)"),
    ("GSE316447", "success/dental-ext", "PDLF a-KG ATAC (ext valid)"),
    ("GSE160451", "neg-control", "DPSC ETV2 reprogramming (NEG ctrl, exclude RRA)"),
    ("GSE105145", "anchor", "DPSC vs iliac BMSC baseline"),
    # success / jaw regeneration positive reference
    ("GSE104473", "success/jaw", "mouse mandible distraction osteogenesis"),
    ("GSE223778", "success/jaw", "extraction socket healing Mertk"),
    # failure / MRONJ-ORNJ
    ("GSE7116",   "failure", "human BRONJ peripheral blood (2007 array; MM confound)"),
    ("GSE303003", "failure", "human MRONJ sequestrum granulation scRNA"),
    ("GSE269255", "failure", "mouse ORNJ mandible BMSC scRNA+metabolome (taurine)"),
    ("GSE295106", "failure", "mouse BRONJ mandible marrow scRNA"),
    ("GSE296096", "failure", "mouse osteoclast bisphosphonate toxicity"),
    ("GSE306512", "failure", "mouse macrophage mitochondria"),
    # jaw position-specific background
    ("GSE58474",  "position", "human mandible vs iliac osteoblasts"),
    ("GSE30167",  "position", "mouse jaw/alveolar vs long bone"),
]

def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "geo-audit/1.0 (academic)"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read().decode("utf-8", "replace")

def esearch_uid(acc):
    q = urllib.parse.urlencode({"db": "gds", "term": f"{acc}[ACCN] AND gse[ETYP]", "retmode": "json", "retmax": 5})
    try:
        d = json.loads(fetch(f"{EUTILS}/esearch.fcgi?{q}"))
        ids = d.get("esearchresult", {}).get("idlist", [])
        # prefer the 200....... uid that ends with the GSE number
        num = acc.replace("GSE", "")
        for i in ids:
            if i.endswith(num.zfill(len(i) - 3)) or i.endswith(num):
                return i
        return ids[0] if ids else None
    except Exception as e:
        return f"ERR:{e}"

def esummary(uid):
    q = urllib.parse.urlencode({"db": "gds", "id": uid, "retmode": "json"})
    d = json.loads(fetch(f"{EUTILS}/esummary.fcgi?{q}"))
    res = d.get("result", {})
    rec = res.get(uid, {})
    return rec

rows = []
for acc, role, note in ACCS:
    uid = esearch_uid(acc)
    time.sleep(0.4)
    if not uid:
        rows.append({"acc": acc, "role": role, "note": note, "status": "NOT_FOUND", "title": "", "n": "", "gpl": "", "taxon": "", "pdat": ""})
        print(f"{acc:11s} NOT_FOUND")
        continue
    if str(uid).startswith("ERR"):
        rows.append({"acc": acc, "role": role, "note": note, "status": uid, "title": "", "n": "", "gpl": "", "taxon": "", "pdat": ""})
        print(f"{acc:11s} {uid}")
        continue
    try:
        rec = esummary(uid)
        time.sleep(0.4)
        title = rec.get("title", "")
        n = rec.get("n_samples", "")
        gpl = rec.get("gpl", "")
        taxon = rec.get("taxon", "")
        pdat = rec.get("pdat", "")
        gdstype = rec.get("gdstype", "")
        rows.append({"acc": acc, "role": role, "note": note, "status": "FOUND",
                     "title": title, "n": n, "gpl": "GPL" + str(gpl) if gpl else "",
                     "taxon": taxon, "pdat": pdat, "gdstype": gdstype})
        print(f"{acc:11s} n={str(n):>3s} {taxon[:22]:22s} {pdat:10s} {gdstype[:28]:28s} | {title[:60]}")
    except Exception as e:
        rows.append({"acc": acc, "role": role, "note": note, "status": f"SUMMARR:{e}", "title": "", "n": "", "gpl": "", "taxon": "", "pdat": ""})
        print(f"{acc:11s} SUMMARY_ERR {e}")

out = sys.argv[1] if len(sys.argv) > 1 else "geo_audit.json"
with open(out, "w", encoding="utf-8") as f:
    json.dump(rows, f, ensure_ascii=False, indent=2)
print(f"\nSaved -> {out}  ({sum(1 for r in rows if r['status']=='FOUND')}/{len(rows)} FOUND)")
