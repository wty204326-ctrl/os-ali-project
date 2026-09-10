# 基于氧化应激相关基因的急性肝损伤机器学习诊断模型

> 计算生物学与人工智能课程项目 · 截止日期 2026-09-20
> 全流程在 Linux (WSL2 Ubuntu) 下完成,一键可复现
> **v1.1 方法学修订版**(修正特征泄漏、跨平台尺度、循环论证与选择性报告;数字以本版为准)

## 📊 结果速览

| 指标 | 结果 |
|---|---|
| 差异表达基因(GSE17649,strain 校正) | 710 个(上调 429 / 下调 281) |
| OS-DEGs(与 GO:0006979 交集) | 13 探针 / 9 基因:Btg1, Gclc, Hmox1, Itch, Rcan1, Rcan2, Srxn1, Txnrd1, Xpa |
| 富集(全部 DEGs,非循环论证) | KEGG 31 条 p.adjust<0.05:MAPK 信号、FoxO 信号、IL-17 信号、凋亡、TNF 信号 |
| 训练集 AUC(**嵌套特征选择**,无泄漏) | randomForest 0.958 / LASSO 0.926 / XGBoost 0.917 |
| 外部验证 AUC(GSE111828,独立 RNA-seq) | LASSO 0.925 / RF 0.888 / XGBoost 0.813 —— **仅代表区分度(排序),概率未校准** |
| 跨物种方向参考(GSE120652,人,n=3v3) | MGI 同源映射下 6/9 一致(GCLC/RCAN1/BTG1 相反);人源无一基因过 FDR |

完整报告见 [`report/06_report.html`](report/06_report.html)(一键渲染:`rmarkdown::render("report/06_report.Rmd")`)。

### 与 v1.0 的方法学差异(重要)

| 问题 | v1.0 | v1.1 |
|---|---|---|
| 特征泄漏 | 全样本选特征后 CV,AUC 偏乐观 | limma 特征选择嵌入每折训练集内部重跑(嵌套 CV) |
| 跨平台尺度 | 直接用芯片模型预测 RNA-seq 计数 | 用训练集均值/标准差对验证集 z-score 标准化 |
| 循环论证 | 对 OS-DEGs 报"氧化应激通路显著" | 富集改在**全部 DEGs** 上做;OS-DEGs 交集不再当证据 |
| 校准 | 称"诊断模型" | 明确区分**区分度(AUC)**与**校准**;概率塌缩如实报告 |
| 跨物种 | 符号大小写硬配,只报一致项 | MGI HOM_MouseHumanSequence 官方同源映射,一致与不一致全列 |
| 品系混杂 | 未纳入设计矩阵 | strain 作为协变量纳入 limma 设计 |
| 环境 | 含 CUDA 构建的锁定文件 | 精简可移植规格 + 依赖自检脚本 |



## 1. 研究背景

对乙酰氨基酚(Acetaminophen, APAP)过量是欧美及国内急性肝衰竭(ALF)最主要的药物性病因,氧化应激(oxidative stress)是其核心致病机制:过量 NAPQI 耗竭谷胱甘肽、线粒体活性氧(ROS)爆发,最终导致肝细胞大面积坏死。

本研究拟基于 GEO 公共转录组数据:

1. 筛选急性肝损伤中的**差异表达基因(DEGs)**(品系协变量校正);
2. 与氧化应激相关基因集取交集,锁定**氧化应激相关差异基因(OS-DEGs)**作为候选特征;
3. 用 LASSO / 随机森林 / XGBoost 构建**急性肝损伤分类模型**(嵌套交叉验证,无特征泄漏);
4. 在独立外部队列上评估**区分度(AUC)与校准**,并如实区分二者;
5. 对全部 DEGs 做富集分析提供独立的生物学解释,并以 MGI 同源映射检验跨物种方向一致性。

## 2. 数据来源

> 注:急性肝损伤非肿瘤,TCGA 无相应数据,全部数据来自 GEO。

| 角色 | 数据集 | 物种 | 设计 | n |
|---|---|---|---|---|
| 训练集 | GSE17649 | 小鼠 | 多品系 APAP 中毒 vs 对照(12:24,均衡) | 36 |
| 外部验证集 | GSE111828 | 小鼠 | C57BL/6J APAP vs 对照(4:20) | 24 |
| 跨物种方向验证 | GSE120652 | 人 | APAP-ALF vs 正常肝(FFPE,GPL6244) | 6 |

氧化应激基因集:GO:0006979 (response to oxidative stress),来源 EBI GOA 注释数据库(UniProtKB 行,排除 NOT 注释;脚本自动下载解析)。跨物种映射使用 **MGI HOM_MouseHumanSequence** 官方同源报告。

**设计说明**:训练集选择样本量最大且组间均衡的 GSE17649,品系作为协变量纳入设计矩阵;GSE111828 为独立实验、不同技术(RNA-seq vs 芯片),用于外部验证并评估跨平台稳健性;人类数据样本少,仅作关键基因方向的跨物种参考。

**对照定义说明**:GSE17649 的对照为同一实验 0h 基线(非独立健康肝),因此识别的是"给药后早期转录变化"而非"健康 vs 损伤";该局限在报告讨论中说明。

## 3. 分析流程

```
GEO series matrix ──→ 预处理(log2) ──→ limma 差异表达(strain 协变量校正,|log2FC|>1, adj.P<0.05)
                                            │
                    GO:0006979 氧化应激基因集 ─┴──→ OS-DEGs 交集(用于建模特征候选)
                                            │
                         全部 DEGs 的 GO/KEGG 富集(clusterProfiler,非循环论证)
                                            │
    嵌套交叉验证:每折训练集内部 limma 选特征 → LASSO/randomForest/XGBoost → 折外预测
                                            │
                      外部验证(z-score 用训练集统计量)+ 区分度/校准分别报告
```

## 4. 目录结构

```
├── README.md                     项目说明(本文件)
├── env/rna.yml                   精简可移植环境规格(不含 CUDA 锁定)
├── scripts/
│   ├── 00_install_deps.sh        依赖自检与补装(含 org.Mm/org.Hs.eg.db 分段安装)
│   ├── 01_download.sh            数据与基因集下载(带完整性校验)
│   ├── 02_de_analysis.R          strain 校正的差异表达 + 探针注释 + OS-DEGs
│   ├── 03_ml_modeling.R          嵌套交叉验证建模(折内选特征)
│   ├── 04_enrichment.R           全部 DEGs 的 GO/KEGG 富集
│   ├── 05_validation.R           外部验证(z-score 标准化 + 校准报告)
│   ├── 06b_crossspecies.R        MGI 同源映射的跨物种方向检验
│   └── run_all.sh                一键复现入口
├── data/raw/                     原始数据(不提交,.gitignore 已排除)
├── results/                      分析结果(表格与图,提交)
└── report/06_report.Rmd|html     课程报告
```

## 5. 复现方法

```bash
# 需先安装 Miniforge(conda/mamba)于 ~/miniforge3
mamba env create -f env/rna.yml
source ~/miniforge3/bin/activate rna
bash scripts/run_all.sh        # 含依赖自检、数据下载、全部分析与报告渲染
```

脚本不写死路径:每个 R 脚本由自身位置推导仓库根目录,可克隆到任意路径运行。
若网络受限导致 org.Mm/org.Hs.eg.db 安装失败,`00_install_deps.sh` 会自动改用
分段并行下载源码包安装。

## 6. 里程碑

- [x] v0.1-data 仓库初始化、数据集核验与下载(2026-09-09)
- [x] v0.2-de 差异表达分析 + OS-DEGs + 富集(2026-09-09)
- [x] v0.3-ml 机器学习建模与外部验证(2026-09-09)
- [x] v1.0 报告、跨物种验证与一键复现脚本(2026-09-10)
- [x] **v1.1 方法学修订**:嵌套 CV 消除特征泄漏、验证集 z-score 标准化、
      富集改全部 DEGs(去循环论证)、MGI 同源映射、strain 协变量、精简环境(2026-09-10)

## 7. 环境与工具

Ubuntu 26.04 (WSL2) · R 4.5.3 · limma/clusterProfiler/org.Mm.eg.db ·
glmnet/randomForest/xgboost(CPU)/pROC · pheatmap / readxl / rmarkdown

## 8. 已知局限(如实声明)

1. 训练集对照为同批实验 0h 基线,非独立健康肝;识别的是给药后早期转录变化。
2. 外部验证对照组仅 4 例,AUC 置信区间宽。
3. LASSO 在外部数据上预测概率整体塌缩(~0),AUC 仅代表排序能力,不可作诊断概率。
4. 人源数据集 n=3v3,无一基因过 FDR,跨物种结论仅为弱参考。
5. 未纳入 ALT/AST 等临床协变量。

