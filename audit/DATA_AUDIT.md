# 数据审计报告 — 颌骨干细胞再生成败对照项目

> 核验日期:2026-06（NCBI E-utilities 实时查询）
> 方法:`audit/verify_geo.py`（esearch→esummary, db=gds/sra），原始结果 `audit/geo_audit_raw.json`
> 原则:**只报告真实查询所得**；样本数 = GEO esummary `n_samples`（系列总样本，含分组，非各臂净样本）；标注与方案 v1 的差异。

## 0. 一句话结论

**27/27 个 GEO accession 全部真实存在；DRA006607 在 NCBI SRA 真实存在；三处"白地/旧范式"声明全部经 PubMed 复核成立。** 审计同时暴露 3 个必须修正的执行风险(失败侧人源数据比预想薄),已在下方逐条给出处置。**成功侧极稳,失败侧须改为"以鼠下颌特异 scRNA 为主、人源为辅、定位假设生成"。**

---

## 1. 成功侧 · 牙源成骨(RRA-meta 候选池)

| accession | 细胞 | n | 物种 | 上线 | 平台类型 | 处置 |
|---|---|---|---|---|---|---|
| GSE99958 | PDLSC | 8 | 人 | 2018 | RNA-seq | ✅ 纳入 |
| GSE163354 | PDLSC (SDF-1/Exendin-4) | 20 | 人 | 2021 | RNA-seq | ✅ 纳入(统计锚点) |
| GSE226347 | DPSC-ECM (Wnt/Hippo) | 12 | 人 | 2024 | RNA-seq | ⚠️ 纳入前**逐臂确认受体细胞**(勿混 ECM 诱导的非 DPSC 受体) |
| GSE296018 | DPSC (菟丝子黄酮) | 6 | 人 | 2026 | RNA-seq | ✅ 纳入(药物×成骨臂,取成骨向对照) |
| GSE271641 | DPSC (CBD) | 6 | 人 | 2025 | RNA-seq | ✅ 纳入(odonto/osteo,裁出成骨臂) |
| GSE299041 | DPSC (剪切力/p38) | 6 | 人 | 2025 | RNA-seq | ✅ 纳入(力学-成骨) |
| GSE286540 | SCAP (辛伐他汀) | 6 | 人 | 2025 | RNA-seq | ✅ 纳入(SCAP 稀缺) |
| GSE236009 | PDLSC (Cald1) | 6 | 人 | 2025 | RNA-seq | ✅ 纳入 |
| GSE266150 | PDLSC (柚皮素) | 8 | 人 | 2024 | RNA-seq | ✅ 纳入(验证集) |
| GSE159507 | PDLSC 成骨 mRNA/lncRNA | 6 | 人 | 2020 | 芯片 | ✅ 纳入(跨平台,RRA 权重) |
| GSE159508 | PDLSC 成骨 miRNA | 6 | 人 | 2020 | 芯片(ncRNA) | ➖ miRNA,作 miRNA-mRNA 旁证,不入 mRNA RRA |
| GSE266257 | **SHED**(无机磷酸盐) | 4 | 人 | 2026 | RNA-seq | ✅ 纳入(已确认 SHED+茜素红成骨/矿化) |
| GSE49007 | DFC 牙囊 D7 | 4 | 人 | 2013 | 芯片 | ✅ 纳入(DFC 稀缺,跨平台权重低) |
| GSE316449 | PDLF (α-KG) | 8 | 人 | 2026 | RNA-seq | ✅ **外部独立成骨验证集**(非 RRA 训练) |
| GSE316447 | PDLF (α-KG) | 8 | 人 | 2026 | ATAC/占位 | ✅ 外部表观验证(配套 GSE316449) |
| GSE160451 | DPSC (ETV2 过表达) | 6 | 人 | 2022 | RNA-seq | ⛔ **阴性对照,不纳入 RRA**(谱系重编程,非成骨诱导) |
| GSE105145 | DPSC vs 髂骨 BMSC | 6 | 人 | 2019 | RNA-seq | ✅ 跨源锚点(非 RRA,作 §4.6 跨源对照) |

**小结**:成功侧牙源数据**全部真实**,覆盖 PDLSC/DPSC/SCAP/SHED/DFC 五类。注意多数 2024-2026 集为"药物/扰动 × 成骨"设计 → **Step0 必须逐臂裁出"成骨诱导 vs 未诱导"纯对比**,真正进入 mRNA RRA 训练的纯臂约 **11–12 套**(剔除 GSE159508 miRNA、GSE160451 阴对照、GSE316449/447 留作外部验证、GSE105145 跨源锚点)。**入选臂最终数按裁剪结果如实报告**。

## 2. 成功侧 · 颌骨再生正参照

| accession | 内容 | n | 物种 | 上线 | 处置 |
|---|---|---|---|---|---|
| GSE104473 | **下颌牵张成骨 RNA-seq+ATAC** | 51 | 鼠 | 2018 | ✅ 成功侧主力,样本量罕见大,跨情境验证核心 |
| GSE223778 | 拔牙窝愈合(Mertk) | 12 | 人+鼠 | 2023 | ✅ 第二成功侧(含人源,价值高) |

## 3. 失败侧 · MRONJ/ORNJ ⚠️(审计重点)

| accession | 内容 | n | 物种 | 上线 | 处置 |
|---|---|---|---|---|---|
| GSE7116 | **MM 患者 ONJ 表征**(双膦酸盐) | 26 | 人 | 2007 | ⚠️ **重混杂**:11 例 MM+ONJ,无健康对照、含 MM/化疗背景。**仅作支持性证据,不作主诊断队列** |
| GSE303003 | 人 MRONJ 死骨周围肉芽 scRNA | **2** | 人 | 2025 | ⚠️ **样本极小**(4 患者 MRONJ vs 根尖囊肿对照,GEO 计 2 个处理组)。**仅机制示意,不做定量队列推断**;且对应 2025-11 deposit,**有抢跑**,本研究须切割于"成败框架"而非复述其 detached-osteoclast 机制 |
| GSE269255 | **下颌 ORNJ 放疗基质细胞 scRNA+代谢组** | 3 | 鼠 | 2024 | ✅ 部位确认为下颌(irradiated mandibular stromal cells);牛磺酸轴,作 ORNJ 失败侧+已发表机制旁证 |
| GSE295106 | **下颌骨髓 BRONJ vs 对照 scRNA** | 4 | 鼠 | 2025 | ✅ **失败侧主力(下颌特异、有对照)** |
| GSE296096 | GGOH 逆转 N-BP 毒性 | 9 | 鼠 | 2026 | ✅ 救援轴正交支持(GGOH) |
| GSE306512 | 巨噬线粒体救援防 ONJ | 6 | 鼠 | 2026 | ✅ 救援轴正交支持(线粒体) |

**关键修正(影响方案)**:人源失败侧数据**比 v1 预想薄得多**(GSE303003 n=2;GSE7116 为 MM 混杂)。⇒
- 失败侧主轴改为**鼠下颌特异 scRNA**(GSE295106 对照 vs BRONJ + GSE269255 ORNJ);
- 人源(GSE7116/GSE303003)**降为佐证**,全程定位"假设生成";
- 救援轴(GGOH/牛磺酸/线粒体)作为**靶点方向的正交文献支持**,强化 Discussion。

## 4. 颌骨位置特异背景

| accession | 内容 | n | 物种 | 处置 |
|---|---|---|---|---|
| GSE58474 | 人下颌 vs 髂骨成骨细胞 | 12 | 人 | ✅ 位置特异主证(成骨弱) |
| GSE30167 | 鼠颌/牙槽 vs 长骨 | 4 | 鼠 | ✅ 跨物种佐证 |
| DRA006607 | 人 maxilla/mandible/ilium MSC 位置记忆 | — | 人 | ✅ 真实存在于 NCBI SRA(uid 10220584);**原始 fastq,需自行比对定量**;DDBJ/SRA 非 GEO |

## 5. 遗传层 / 药物层(公开成熟资源,未逐一再核,均为领域标准)

- eBMD GWAS — Morris JA & Kemp JP et al. *Nat Genet* 2019(GEFOS,GWAS Catalog 可取汇总统计)
- FinnGen(R10+)骨质疏松/骨折端点;GEFOS DXA-BMD(复制)
- eQTLGen 血液 cis-eQTL(n=31,684);GTEx v8
- Finan 2017 druggable genome(*Sci Transl Med*);Open Targets / DGIdb / ChEMBL
- CMap/LINCS L1000(clue.io);L1000CDS2;L1000FWD
- MSigDB Hallmark/GO-BP osteoblast differentiation;WikiPathways;DrugBank

## 6. 三处白地/旧范式 — PubMed 复核结果(论文新颖性命门)

| 声明 | 检索 | 结果 | 结论 |
|---|---|---|---|
| **MRONJ 特异孟德尔随机化 = 白地** | `(MRONJ OR osteonecrosis of the jaw) AND Mendelian randomization` | 全库 **1** 命中,唯一命中 PMID 39749391 = "特应性皮炎与痴呆"(词义爆炸误命中,与颌骨坏死无关) | ✅ **真白地,实质 0 篇** |
| **跨牙源成骨 RRA-meta = 白地** | `dental (pulp/PDL) stem cell AND osteogenic AND robust rank aggregation` | 全库 **0** 命中 | ✅ **真白地** |
| **CMap→成骨→老药重定位 = 旧范式** | `parbendazole AND osteogenic` | PMID **26420877**(Brey, *PNAS* 2015, "Connectivity Map-based discovery of parbendazole...")+ PMID 29194609 跟进 | ✅ **确为占据领地 → CMap 降级为 coda 正确** |

> 创新点据此最终锁定为:**① 再生"成功 vs 失败"统一对照框架(以下颌牵张成骨为成功正参照、MRONJ/ORNJ 为失败)+ ② 泛牙源 RRA"再生胜任程序"(0 先例)+ ③ MRONJ 语境的 druggable-genome MR(0 先例)**;CMap 仅作与 MR 取交集的"双证据"探索性 coda。

## 7. 对方案 v1 的净修正清单

1. **失败侧重心** 人源→鼠下颌特异 scRNA(GSE295106 主力 + GSE269255 ORNJ);人源 GSE7116/GSE303003 降为佐证。
2. **GSE303003** n=2、且为潜在抢跑源 → 仅机制示意,创新切割于"成败框架"。
3. **GSE7116** 为 MM 混杂、无健康对照 → 不作主诊断队列;若需人源失败 DE,只能在 MM 内做 ONJ vs non-ONJ 且重声明混杂。
4. **GSE266257=SHED 成骨**、**GSE269255=下颌 ORNJ** 两处疑虑解除(均确认)。
5. **RRA 纯训练臂** 经裁剪后约 11–12 套(非 14),最终数按 Step0 结果报告。
6. **DRA006607** 为 SRA 原始数据,需纳入"需自行比对"的工作量预算。
7. 三处白地坐实 → 可放心写进 Introduction 的 gap statement(附本审计为证)。
