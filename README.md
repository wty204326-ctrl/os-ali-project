# 基于氧化应激相关基因的急性肝损伤机器学习诊断模型

> 课程项目:计算生物学与人工智能课程项目 · 截止日期 2026-09-20
> 全流程在 Linux (WSL2 Ubuntu) 下完成,一键可复现

## 结果速览(2026-09-13 重跑,分层 CV + z-score 版)

| 指标 | 结果 |
|---|---|
| 差异表达基因(GSE17649) | 566 个(上调 361 / 下调 205) |
| OS-DEGs(与 GO:0006979 交集) | 12 探针 / 8 基因:Hmox1, Srxn1, Gclc, Txnrd1, Rcan1, Btg1, Itch, Xpa |
| 训练集 AUC(分层 5×5 折 CV,折内 z-score) | randomForest 0.954 / XGBoost 0.953 / LASSO 0.927 |
| **外部验证 AUC(GSE111828,小鼠 RNA-seq)** | **LASSO 0.975** / randomForest 0.956 / XGBoost 0.919 |
| 外部验证阈值 0.5(LASSO) | 灵敏度 0.70 / 特异度 1.00(4 对照全部判对) |
| 跨物种方向一致性(GSE120652,人,n=6) | 5/8 方向一致,均未过 FDR(探索性) |

完整报告见 [`report/06_report.html`](report/06_report.html)(一键渲染:`rmarkdown::render("report/06_report.Rmd")`)。

## 1. 研究背景

对乙酰氨基酚(Acetaminophen, APAP)过量是欧美及国内急性肝衰竭(ALF)最主要的药物性病因,氧化应激(oxidative stress)是其核心致病机制:过量 NAPQI 耗竭谷胱甘肽、线粒体活性氧(ROS)爆发,最终导致肝细胞大面积坏死。

本研究拟基于 GEO 公共转录组数据:

1. 筛选急性肝损伤中的**差异表达基因(DEGs)**;
2. 与氧化应激相关基因集取交集,锁定**氧化应激相关差异基因(OS-DEGs)**;
3. 用机器学习方法(LASSO / 随机森林 / XGBoost)构建**急性肝损伤诊断模型**;
4. ROC/AUC 评估并在独立外部数据集上验证;
5. 对关键基因做富集分析与文献佐证,提供生物学解释。

## 2. 数据来源

> 注:急性肝损伤非肿瘤,TCGA 无相应数据,全部数据来自 GEO。

| 角色 | 数据集 | 物种 | 设计 | n |
|---|---|---|---|---|
| 训练集 | GSE17649 | 小鼠 | 多品系 APAP 后 0h vs 3–6h(12:24) | 36 |
| 外部验证集 | GSE111828 | 小鼠 | C57BL/6J 未处理 vs APAP 12–72h(4:20) | 24 |
| 跨物种方向验证 | GSE120652 | 人 | APAP-ALF vs 正常肝(FFPE,GPL6244) | 6 |

氧化应激基因集:GO:0006979 (response to oxidative stress),来源 EBI GOA 直接注释(脚本自动下载解析)。

**设计说明**:训练集对照为给药后 0h 基线,并含抗性品系 SJL。GSE111828 为**同物种、不同平台**的独立实验,用于跨平台验证。人类数据样本少,仅用于关键基因表达方向的探索性讨论。

## 3. 分析流程

```
GEO series matrix ──→ 预处理 / 必要时 log2 ──→ limma(|log2FC|>1, adj.P<0.05)
GPL annot.gz ──→ 探针注释 ─────────────────→ OS-DEGs(与 GO:0006979 交集)
                                                        │
                              GO/KEGG 富集(clusterProfiler;GO:0006979 不作独立发现)
                                                        │
              LASSO ─┐   分层 CV + 折内 z-score         │
        randomForest ─┼── 7 个 OS-DEG 基因 ─────────────→ 诊断模型
             XGBoost ─┘                                  ↓
                         验证集独立 z-score → ROC/AUC → GSE111828
```

## 4. 目录结构

```
├── README.md            项目说明(本文件)
├── env/rna.yml          conda 锁文件(Linux-64)
├── scripts/             全部分析脚本(按编号顺序执行)
├── data/raw/            原始数据(不提交,.gitignore 已排除)
├── results/             分析结果(表格与图,提交)
└── report/              R Markdown 最终报告
```

## 5. 复现方法

```bash
# 1. 创建环境(Ubuntu, conda/mamba)
mamba env create -f env/rna.yml
conda activate rna

# 2. 依次运行脚本(会自动设置项目根目录)
bash scripts/run_all.sh
```

`env/rna.yml` 从本机 Linux 环境导出,含 CUDA 相关构建。若 `mamba env create` 失败,在已有 R 4.x 中安装 limma、clusterProfiler、org.Mm.eg.db、glmnet、randomForest、xgboost、pROC、ggplot2、rmarkdown、readxl 后直接运行脚本即可。

## 6. 里程碑

- [x] v0.1-data 仓库初始化、数据集核验与下载(2026-09-09)
- [x] v0.2-de 差异表达分析 + OS-DEGs + 富集(2026-09-09)
- [x] v0.3-ml 机器学习建模与外部验证(2026-09-09)
- [x] v1.0-final 报告、跨物种验证与一键复现脚本(2026-09-09)
- [x] v1.1-fix 可复现路径、自动探针注释、分层 CV、跨平台 z-score、报告表述(2026-09-13)

## 7. 环境与工具

Ubuntu (WSL2) · R 4.x · limma/clusterProfiler · glmnet/randomForest/xgboost/pROC
