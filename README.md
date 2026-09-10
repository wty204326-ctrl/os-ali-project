# 基于氧化应激相关基因的急性肝损伤机器学习诊断模型

> 课程项目:计算生物学与人工智能课程项目 · 截止日期 2026-09-20
> 全流程在 Linux (WSL2 Ubuntu) 下完成,一键可复现

## 📊 结果速览

| 指标 | 结果 |
|---|---|
| 差异表达基因(GSE17649) | 566 个(上调 361 / 下调 205) |
| OS-DEGs(与 GO:0006979 交集) | 9 探针 / 7 基因:Hmox1, Srxn1, Gclc, Txnrd1, Rcan1, Btg1, Itch |
| 最显著富集通路 | GO:0006979 response to oxidative stress(FE=38, adj.P=3.3e-09) |
| 训练集 AUC(5×5 折 CV) | randomForest 0.957 / XGBoost 0.944 / LASSO 0.924 |
| **外部验证 AUC(GSE111828,独立 RNA-seq)** | **LASSO 0.938** / XGBoost 0.913 / RF 0.831 |
| 跨物种方向一致性(GSE120652,人) | Nrf2 核心模块(HMOX1/SRXN1/TXNRD1)方向全部一致 |

完整报告见 [`report/06_report.html`](report/06_report.html)(一键渲染:`rmarkdown::render("report/06_report.Rmd")`)。


## 1. 研究背景

对乙酰氨基酚(Acetaminophen, APAP)过量是欧美及国内急性肝衰竭(ALF)最主要的药物性病因,氧化应激(oxidative stress)是其核心致病机制:过量 NAPQI 耗竭谷胱甘肽、线粒体活性氧(ROS)爆发,最终导致肝细胞大面积坏死。

本研究拟基于 GEO 公共转录组数据:

1. 筛选急性肝损伤中的**差异表达基因(DEGs)**;
2. 与氧化应激相关基因集取交集,锁定**氧化应激相关差异基因(OS-DEGs)**;
3. 用机器学习方法(LASSO / 随机森林 / XGBoost)构建**急性肝损伤诊断模型**;
4. ROC/AUC 评估并在独立外部数据集上验证;
5. 对关键基因(hub genes)做富集分析与文献佐证,提供生物学解释。

## 2. 数据来源

> 注:急性肝损伤非肿瘤,TCGA 无相应数据,全部数据来自 GEO。

| 角色 | 数据集 | 物种 | 设计 | n |
|---|---|---|---|---|
| 训练集 | GSE17649 | 小鼠 | 多品系 APAP 中毒 vs 对照(12:24,均衡) | 36 |
| 外部验证集 | GSE111828 | 小鼠 | C57BL/6J APAP vs 对照(4:20) | 24 |
| 跨物种方向验证 | GSE120652 | 人 | APAP-ALF vs 正常肝(FFPE,GPL6244) | 6 |

氧化应激基因集:GO:0006979 (response to oxidative stress),来源 EBI GOA 注释数据库(脚本自动下载解析)。

**设计说明**:训练集选择样本量最大且组间均衡的 GSE17649;GSE111828 为独立实验、不同平台(GPL19057),用于外部验证并评估跨平台稳健性;人类数据样本少,用于关键基因表达方向的跨物种一致性讨论。

## 3. 分析流程

```
GEO series matrix ──→ 预处理/标准化(替代log2) ──→ limma 差异表达(|log2FC|>1, adj.P<0.05)
                                                        │
                    GO:0006979 氧化应激基因集 ────────→ OS-DEGs 交集
                                                        │
                              GO/KEGG 富集(clusterProfiler)
                                                        │
              LASSO ─┐                                  │
        randomForest ─┼── 特征筛选 → hub genes → 合并特征  │
             XGBoost ─┘                                  ↓
                                        诊断模型 → ROC/AUC → 外部验证集
```

## 4. 目录结构

```
├── README.md            项目说明(本文件)
├── env/rna.yml          conda 环境定义(可复现)
├── scripts/             全部分析脚本(按编号顺序执行)
├── data/raw/            原始数据(不提交,.gitignore 已排除)
├── results/             分析结果(表格与图,提交)
└── report/              R Markdown 最终报告
```

## 5. 复现方法

```bash
# 1. 创建环境(Ubuntu, conda/mamba)
mamba env create -f env/rna.yml
source ~/miniforge3/bin/activate rna

# 2. 依次运行脚本
bash scripts/run_all.sh
```

## 6. 里程碑

- [x] v0.1-data 仓库初始化、数据集核验与下载(2026-09-09)
- [x] v0.2-de 差异表达分析 + OS-DEGs + 富集(2026-09-09,提前 4 天)
- [x] v0.3-ml 机器学习建模与外部验证(2026-09-09)
- [x] v1.0-final 报告、跨物种验证与一键复现脚本(2026-09-09)

## 7. 环境与工具

Ubuntu 26.04 (WSL2) · R 4.5.3 · DESeq2/limma/clusterProfiler · glmnet/randomForest/xgboost/pROC
