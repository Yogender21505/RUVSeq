# : Remove Unwanted Variation from RNA-Seq Data Using R Programming
In this document, we show how to conduct a differential expression (DE) analysis that
controls for “unwanted variation”, e.g., batch, library preparation, and other nuisance
effects, using the between-sample normalization methods proposed in . We call this
approach RUVSeq for remove unwanted variation from RNA-Seq data.
Briefly, RUVSeq works as follows. 
For n samples and J genes
consider the following:
generalized linear model (GLM), where the RNA-Seq read counts are regressed on both
the known covariates of interest and unknown factors of unwanted variation,
# log E[Y |W, X, O] = W α + Xβ + O. 1
Here, Y is the n × J matrix of observed gene-level read counts, W is an n × k matrix
corresponding to the factors of “unwanted variation” and α its associated k × J matrix
of nuisance parameters, X is an n × p matrix corresponding to the p covariates of
interest/factors of “wanted variation” (e.g., treatment effect) and β its associated p × J
matrix of parameters of interest, and O is an n×J matrix of offsets that can either be set
to zero or estimated with some other normalization procedure (such as upper-quartile
normalization).
The matrix X is a random variable, assumed to be known a priori. For instance, in the
usual two-class comparison setting (e.g., treated vs. control samples), X is an n × 2
design matrix with a column of ones corresponding to an intercept and a column of
indicator variables for the class of each sample (e.g., 0 for control and 1 for treated)
. The simultaneous estimation of W, α, β, and k is infeasible. For a given k, we consider
instead the following three approaches to estimate the factors of unwanted variation W:
• RUVg uses negative control genes, assumed to have constant expression across
samples;
• RUVs uses centered (technical) replicate/negative control samples for which the
covariates of interest are constant;
• RUVr uses residuals, e.g., from a first-pass GLM regression of the counts on the
covariates of interest.
The resulting estimate of W can then be plugged into Equation 1 , for the full set of genes
and samples, and α and β estimated by GLM regression. Normalized read counts can
be obtained as residuals from ordinary least squares (OLS) regression of log Y − O on
the estimated W.
Note that although here we illustrate the RUV approach using the GLM implementation
of edgeR and DESeq2, all three RUV versions can be readily adapted to work with any
DE method formulated within a GLM framework.
See Below for full details and algorithms for each of the three RUV procedures
[Normalized_PCA.pdf](https://github.com/Yogender21505/RUVSeq/files/10133826/Normalized_PCA.pdf)
[Normalized_RLE.pdf](https://github.com/Yogender21505/RUVSeq/files/10133827/Normalized_RLE.pdf)
[spike_PCA_plot.pdf](https://github.com/Yogender21505/RUVSeq/files/10133828/spike_PCA_plot.pdf)
[spike_RLE_plot.pdf](https://github.com/Yogender21505/RUVSeq/files/10133829/spike_RLE_plot.pdf)
[variability_PCA.pdf](https://github.com/Yogender21505/RUVSeq/files/10133831/variability_PCA.pdf)
[Varibility.pdf](https://github.com/Yogender21505/RUVSeq/files/10133832/Varibility.pdf)
[VolcanoPlot.pdf](https://github.com/Yogender21505/RUVSeq/files/10133833/VolcanoPlot.pdf)
