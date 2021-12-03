import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from bioinfokit import analys, visuz

from rpy2 import robjects
from rpy2.robjects import Formula

from rpy2.robjects import pandas2ri
pandas2ri.activate()

from rpy2.robjects.packages import importr

base = importr("base")
stats = importr("stats")
DESeq2 = importr("DESeq2")

SRR3414629 = pd.read_csv("SRR3414629.counts", sep="\t", names=['r1'], index_col=0)
SRR3414630 = pd.read_csv("SRR3414630.counts", sep="\t", names=['r2'], index_col=0)
SRR3414631 = pd.read_csv("SRR3414631.counts", sep="\t", names=['r3'], index_col=0)
SRR3414635 = pd.read_csv("SRR3414635.counts", sep="\t", names=['c1'], index_col=0)
SRR3414636 = pd.read_csv("SRR3414636.counts", sep="\t", names=['c2'], index_col=0)
SRR3414637 = pd.read_csv("SRR3414637.counts", sep="\t", names=['c3'], index_col=0)
counts = pd.read_csv("SRR3414629.counts", sep="\t", names=['r1'], index_col=0)
counts['r2'] = SRR3414630['r2']
counts['r3'] = SRR3414631['r3']
counts['c1'] = SRR3414635['c1']
counts['c2'] = SRR3414636['c2']
counts['c3'] = SRR3414637['c3']
counts.to_csv("ALL.counts", sep="\t")
# counts = pd.read_csv("ALL.counts", sep="\t", index_col=0)

# Define meta
meta = pd.DataFrame({"Type": ["Sample"]*3 + ["Control"]*3}, index=counts.columns)
meta["Type"] = stats.relevel(robjects.vectors.FactorVector(meta["Type"]), ref="Control")
print(meta)

# Calculate normalization factors
dds = DESeq2.DESeqDataSetFromMatrix(countData=counts, colData=meta, design=Formula("~ Type"))
dds = DESeq2.DESeq(dds)
print(dds)

res = DESeq2.results(dds, name="Sample_vs_Control")
res = DESeq2.lfcShrink(dds, coef="Sample_vs_Control", type="apeglm")
res = pd.DataFrame(base.as_data_frame(res))
res.index = counts.index
res = res.sort_values("padj")
res = res.loc[res["padj"] < 0.01]
# res = res.loc[res["log2FoldChange"].abs() >= 1]
print(res)

res.to_csv("differentially_expressed_genes.txt", sep="\t")
'''
sns.heatmap(res)
plt.savefig("HW3_heatmap.pdf")
plt.close()

# visuz.GeneExpression.ma(df=res, lfc='log2FC', ct_count='value1', st_count='value2', pv='p-value')
'''