#!/usr/bin/Rscript --vanilla

#########################################
## Program: Differential_gene_expression_analysis_on_counts.R
## Purpose: Perform differential gene expression analysis on a feature x sample count data matrix which might be originated from STAR-HTSeq/STAR-Salmon/STAR-RSEM pipeline.
##	This tool supports only count data.
##  It performs the following major tasks:
#   1. Data filter based on abundance
#   2. Normalization and Transformation
#   3. Sample-level QC
#	4. Gene-level QC
#   5. Identification of Differentially Expressed Genes (DEGs)
##
## Usage: /path/to/Rscript --vanilla /path/to/Differential_gene_expression_analysis_on_counts.R --DGE_tool <Tool for DGE; DESeq2/edgeR> --count_file <path/to/expression data> --gene_id_column <column number for gene identifiers, e.g., 1> --count_column <starting column number for expression data, e.g., 3> --use_CPM_cutoff <whether to filter count data based on CPM calculation, e.g., TRUE> --sample_coverage_threshold <sample coverage threshold, e.g., 0.5> --normalization_method <method to normalize data such as Median/DESeq2-median-of-ratios/edgeR-TMM> --pseudocounts <prior count, e.g., 0> --sample_info_file <path/to/sample metadata> --sample_id_column <column number for sampleids, e.g., 1> --primary_factor_column <column number for sample conditions, e.g., 5> --other_known_factors_columns <column numbers for other known factors which might act as sources of variation in gene expression> --contrasts <list of comma separated contrasts in the format: cond2_vs_cond1> --counts_threshold <threshold for count data filter, e.g., 10> --DESeq2_dispersion_method <method to estimate dispersion such as parametric/local/mean/glmGamPoi> --DESeq2_test_type <type of statistical test such as Wald/LRT> --stat_sig_measure <type of measurement for statistical significance, e.g., pvalue, padj or svalue> --stat_cutoff <cutoff for calling differentially expressed features, e.g., 0.05> --log2_fold_change_cutoff <log2-fold-change cutoff for calling differentially expressed fatures, e.g., 1> --top_n_for_viz <n most variable genes will be considered for generating following plots: clustered heatmap, MDS, e.g., 500> --output_dir <path/to/output directory> --output_prefix <prefix for output files such as figures> --log_file <filename for log>
##
## Developer: Srikant Verma (srikant_verma@persistent.com)
## Version: 1.0: Created on February 28, 2022
## Version: 1.1: Updated on March 9, 2022
## Version: 1.2: Updated on March 11, 2022
## Version: 1.3: Updated on April 4, 2022
## Version: 1.4: Updated on May 9, 2022
## Version: 1.5: Updated on June 26, 2022
#########################################

## Version: 1.0: Created on February 28, 2022
#	1. Data filtering based on abundance
#	2. Normalization (DESeq2's median-of-ratio) and transformation (VST/rlog)
#	3. Sample-level QC: Density plot, MDS, Clustered Heatmap
#	4. Gene-level QC: Dispersion estimation plot
#	5. DGE analysis using DESeq2 only: pairwise comparision using Wald statistics (DESeq2_test_type: Wald)
#	6. log2FC estimation using 'apeglm'
#	7. Support for MAplot, Volcano plot
#	

## Version: 1.1: Updated on March 9, 2022
#	1. Support for DESeq2_test_type: LRT
#	2. User can choose padj or svalue as measurement for statistical significance via --stat_sig_measure
#	3. Adjust the effect of other known sources of variation such as batch, gender etc. via --other_known_factors_columns
#	4. No Support for interaction term in design formula. Users are requested to follow a much simpler approach to achieve their goal which includes combining the factors of interest into a single factor with all combinations of the original factors, and then changing the design to include just this factor, e.g., ~ group. [Reference: 25]. However, this approach could be only usefull in pairwise comparison, such as DESeq2's Wald test.
#	5. --condition_column is replaced with --primary_factor_column


## Version: 1.2: Updated on March 11, 2022
#	1. Remove --full_design and --reduced_design parameters as these can easily created based on --primary_factor_column and --other_known_factors_columns
#	2. Add --gene_id_column by which user can specify the column number (in counts data) to be used as identifier during DGE analysis. This is usefull when the counts data file contains multiple columns for gene annotation, for example, in case of pipeline involving Salmon and tximport.
#	3. 'pvalue' is added as an option for --stat_sig_measure, though it is not recommended.
#	4. drawIHeatmap() updated to handle cases where fixed column or row order is required.


## Version: 1.3: Updated on April 4, 2022
#	1. Added new parameter --use_CPM_cutoff. Now, pre-filtering on counts data will be based on combination of --counts_threshold and --use_CPM_cutoff. If --use_CPM_cutoff is FALSE and --counts_threshold is 10, then only --counts_threshold will be used to keep genes that have at least 10 reads in total. However, if --use_CPM_cutoff is TRUE and --counts_threshold is 10 and library size of least sequenced sample is 10,000,000, then biological meaningful expression level will be (10/10000000)*1000000 = 1 CPM, and therefore along with --sample_coverage_threshold, it will be decided which genes will be kept for downstram analysis.
#	2. Include --other_known_factors_columns into column annotation of heatmaps
#	3. Re-introduced --full_design and --reduced_design parameters to support the study of interaction terms in case of DESeq2's LRT. This will help to address a biological question like - whether a condition induces a change in gene expression at any time point after the reference level time point. [References: 25, 26]
#	4. fitType = "glmGamPoi", which performs faster dispersion and parameter estimation, is no longer mandatory for LRT.
#	5. Identifying gene clusters exhibiting particular patterns across samples (on DESeq2's LRT output) [inspired by DEGreport::degPatterns(), references: 24, 27]
#	6. To create a factor, following code is used: factor(x, levels = unique(x))

## Version: 1.4: Updated on May 9, 2022
#	1. Output files such as list of DEGs and various plots which are generated based on values of --stat_sig_measure, --stat_cutoff and --log2_fold_change_cutoff, will now be stored in sub-directories, for example: log2FC2_padj0.01.


## Version: 1.5: Updated on June 26, 2022
#	Enhancements:
#		1. Removed edgeR from list of supported DE tools from description of the parameter --DGE_tool.
#		2. Validate --gene_id_column, --count_column
#		3. Handle spaces (around comma for more than one input values) in following parameters: --contrasts/--other_known_factors_columns/--stat_cutoff/--log2_fold_change_cutoff/--top_n_for_viz
#		4. Send more error messages to STDOUT (which were previously directed to log files only).
#
#	Bug fixes:
#		1. Explicite use of '--other_known_factors_columns NULL' on command-line gives error "Not all other_known_factors_columns exist in sample_info_file". How NULL went so deep inside the script?  Because make_option(.., type = "character", ...) is converting the command line value (in this case NULL) to a character and assigning it to the variable. Thus it is no longer a NULL and is passing through 'if(!is.null(arguments$other_known_factors_columns)) {other_known_factors_columns <- ...}'. Similar issues might appear for other variables wherever default is NULL and we are using the parameter explicitely on ther command line as argument. The current solution to this problem is to add a check like 'if(other_known_factors_columns == 'NULLL')'.
#		2. Following error encountered in a run with --DESeq_test_type = LRT:
#			<simpleError in normalized_transformed_mat[rownames(res_sig_genes), ]: subscript out of bounds>
#			Issue was rownames of 'res' object got removed at 'res <- cbind(gene_info_df[rownames(res),], res)' where res was a DESeqDataSet object. The solution is to first convert 'res' to a data.frame object using 'res <- data.frame(res)'
#
#		3. Remove 'save(list = ls(), file = …)' from warning/error blocks of tryCatch() because it only saves warning/error message.


## Future releases will have the following features:
#	2. Identify unknown sources of variation using SVA, and adjust their effects in differential expression
#	3. edgeR as DGE_tool
#	4. support for other normalization methods: Median, edgeR-TMM

## Known issues:
#	1. color legends for samples are not consistent across all plots. glimmaMDSplot() uses its own color scheme. Unable to set sample.cols in glimmaXY().



## Notes:
## Usually, all required libraries are called the the beginning of the script. However, this is not the case here as we want the user to quickly check the help page, instead of waiting to see all libraries getting loaded.

## Filtering counts data by abundance:
#	Genes that are not expressed at a biologically meaningful level in any condition should be discarded to reduce the subset of gene to those that are of interest, and to reduce the number of tests carried out downstream when looking at differential expression. Although any sensible value can be used as the expression cutoff, typically a CPM value of 1 is used in analyses as it separates expressed genes from unexpressed genes well for most datasets. Here, a CPM value of 1 means that a gene is expressed if it has at least 20 counts in the sample with library size approx. 20 million or at least 76 counts in the sample with library size approx. 76 million. If sequence reads are summarised by exons rather than genes and/or experiments have low sequencing depth, a lower CPM cutoff may be considered. Additionally, a gene must be expressed in at least 3 samples across the entire experiment to be kept for downstream analysis. [Reference: 1]
#	However, the cpm_cutoff = 1 is a data-independent threshold which may remove genes with potential signal when sequencing depth is high. An alternate is: "The filtering rule is to keep only those genes with n or more samples with CPM greater than the CPM value for a raw count of 10 for the least sequenced sample." For an experiment with the least library size, 15313570, we can set cpmCutoff to (10/15313570) * 1000000 = 0.6530156 [Reference: 2]. Instead of number, we can define a percentage cutoff for the samples. For example, filter out genes which are lowly expressed in less than 80% of a sample group.
#	cpm cutff calculation - Reference: 3

# Filering by adundance should be preferably performed before normalization. [References: 4, 5, 6]

## Normalizaion:
#	Popular methods used for normalization of Counts data: DESeq (RLE/median-of-ratios), TMM-edgeR, FPKM-CuffDiff, total count (TC), Med, UQ and full quantile (FQ) [Reference: 7]. Median normalization has been found as most commonly observed component in top-performing RNA-seq pipelines for various applications [Reference: 8]. CPM has been used as data normalization method in MammaPrint and BluePrint molecular diagnostics using NGS [Reference: 9]. Which method is suitable for between sample comparisons or differential expression and which are for within sample comparions [Reference: 10]. 
#	Some people call CPM, FPKM as data normalization methods [Reference: 11, 12] while some call these data transformation methods [Reference: 13].
#	Popular method for count data transformation is: log2(n + pseudocounts). Two alternative approaches are rlog and VST which offer more rational way of choosing parameters equivalent to pseudocounts. Both transformations produce transformed data on the log2 scale which has been normalized with respect to library size or other normalization factors. If size facors or normalization factors are not present, rlog() and varianceStabilizingTransformation() will call DESeq2::estimateSizeFactors() internally. [Reference: 14]

## Transformation:
#	When the expected amount of variance is approximately the same across different mean values, the data is said to be homoskedastic.
# 	For RNA-seq counts, however, the expected variance grows with the mean. As a solution, DESeq2 offers two transformations for count data that stabilize the variance across the mean: the regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014), and the variance stabilizing transformation (VST) for negative binomial data with a dispersion-mean trend (Anders and Huber 2010), implemented in the vst function.
#	Note that neither rlog transformation nor the VST are used by the differential expression estimation in DESeq, which always occurs on the raw count data
#	The rlog() function transforms the count data to the log2 scale in a way which minimizes differences between samples for rows with small counts, and which normalizes with respect to library size.
#	The varianceStabilizingTransformation() function calculates a variance stabilizing transformation(VST) from the fitted dispersion-mean relation(s) and then transforms the count data (normalized by division by the size factors or normalization factors), yielding a matrix of values which are now approximately homoskedastic (having constant variance along the range of mean values). The transformation also normalizes with respect to library size.
#	How are both rld and VST better than ordinary log2 transformation? - For miRNAs with high counts, the rlog and VST will give similar result to the ordinary log2 transformation of normalized counts. For miRNAs with lower counts, however, the values are shrunken towards the miRNAs’ averages across all samples. The rlog-transformed or VST data then becomes approximately homoskedastic, and can be used directly for computing distances between samples, making PCA plots, or as input to downstream methods which perform best with homoskedastic data.
#	Which transformation to choose? - The rlog tends to work well on small datasets (n < 30), sometimes outperforming the VST when there is a large range of sequencing depth across samples (an order of magnitude difference). The VST is much faster to compute and is less sensitive to high count outliers than the rlog. Therefore the recommended choice is the VST for large datasets (hundreds of samples).
#	So what does the transformed data contain? - In all cases (rlog or VST), the transformation is scaled such that for large counts, it becomes asymptotically (for large values) equal to the logarithm to base 2 of normalized counts.
#	The parameter blind in varianceStabilizingTransformation() or rlog() is logical. It is used to mark whether to blind the transformation to the          experimental design (which is typically done while creating dds object usually using design parameter in DESeqDataSetFromMatrix function). blind=TRUE should be used for comparing samples in a manner unbiased by prior information on samples, for example to perform sample QA (quality assurance). blind=FALSE should be used for transforming data for downstream analysis, where the full use of the design information should be made. blind=FALSE will skip re-estimation of the dispersion trend, if this has already been calculated. If many of genes have large differences in counts due to the experimental design, it is important to set blind=FALSE for downstream analysis. [Reference: ?DESeqDataSetFromMatrix]

## Mean-Variance relationship: Plot the SD of the transformed expression data across samples (of each gene), against the mean. Note that the vertical axis in such plots is the square root of the variance over all samples. The ideal situation: the standard deviation should be roughly constant along the whole dynamic range; the red line depicting the median estimator should be approximately horizontal. However, this may be unreasonable in the case of datasets with many true differences due to the experimental conditions.

## prior.count in edgeR::cpm() only comes into picture when log = TRUE.

## base::as.name is a way to refer to R objects by name, rather than value. It helps in case we want to use non-standard column names, such as those containing hyphen.

## Sample-level QC:
#	Density plot: [Reference: 19]
#		It is a curve on a graph that represents the distribution of values in a dataset. It is useful for three reasons:
#			It gives us good idea of the "shape" of a distribution, including whether or not a distribution has one or more peaks of frequently occurring values and whether or not the distribution is skewed to the left or right.
#			It lets us visually see whether the mean and median of a distribution are located.
#			It lets us visually see what percentage of data points are falling between different values. 

#	PCA plot:
#		It shows unsupervised clustering of the samples in 2D as obtained by PCA. It helps in determining intra- and inter-group variability and outliers. When assessing variability within the dataset, it is preferable that the intergroup variability, representing differences between experimental conditions in comparison with control conditions, is greater than the intragroup variability, representing technical or biological variability. A global overview of the data allows for the characterization of variation between replicates and whether investigator-defined experimental groups show actual differences between groups (a group being a set of replicates from the same condition or of the same cell type).

#	Hierarchical clustering with heatmap:
#		Unsupervised clustering of samples based on pairwise distance. It helps to determine whether the main factor of interest (say, condition) is the most prominent source of variation in the data. It is typically done on Z-scaled transformed data. Z-scaling is done per gene across all samples. This Z-transformed data is also suitable for gene clustering. 
#		"The Euclidean distance of the z-scores is the same as correlation distance ... The most popular methods for gene expression data are to use log2(expression + 0.25), correlation distance and complete linkage clustering." [Reference: 20]

#	MDS plot: [Reference: 17]
#		The multidimensional scaling (MDS) plot is frequently used to explore differences in samples. When data has been MDS transformed, the first two dimensions explain the greatest variance between samples, and the amount of variance decreases monotonically with increasing dimension.
#		The Glimma MDS contains two main components:
#			1)a plot showing two MDS dimensions
#			2)a plot of the eigenvalues of each dimension which indicates how much variation the dimension can explain
#		By default, Glimma uses 500 most variable genes selected based on pairwise distance between samples.


## Gene-level QC: 
#	Gene-wise dispersion plot: [Reference: 15]
#		To examine whether the data is a good fit for the DESeq2 model. You expect your data to generally scatter around the (fitted) curve, with the dispersion decreasing with increasing mean expression levels. If you see a cloud or different shapes, then you might want to explore your data more to see if you have contamination (mitochondrial, etc.) or outlier samples.


## Identification of DE genes:
#	What are the different types of hypothesis tests? "DESeq2 offers two kinds of hypothesis tests: the Wald test, where we use the estimated standard error of a log2 fold change to test if it is equal to zero, and the likelihood ratio test (LRT). The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model, in which some of the terms of the full model are removed. The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero. The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor at once, or all interactions between two variables. The LRT for count data is conceptually similar to an analysis of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a reduced model" [Reference: 23]
#
#	What do the functions nbinomWaldTest() and nbinomLRT() from DESeq2 do?
#	nbinomWaldTest - "This function tests for significance of coefficients in a Negative Binomial GLM, using previously calculated ‘sizeFactors’ (or ‘normalizationFactors’) and dispersion estimates.  See ‘DESeq’ for the GLM formula." [?nbinomWaldTest]
#	nbinomLRT - "This function tests for significance of change in deviance between a full and reduced model which are provided as ‘formula’. Fitting uses previously calculated ‘sizeFactors’ (or ‘normalizationFactors’) and dispersion estimates." [?nbinomLRT]
#
#	What does the function results() do?
#	"‘results’ extracts a result table from a DESeq analysis giving base means across samples, log2 fold changes, standard errors, test statistics, p-values and adjusted p-values; ‘resultsNames’ returns the names of the estimated effects (coefficents) of the model; ‘removeResults’ returns a ‘DESeqDataSet’ object with     results columns removed." [?results]
#
#	Shrinking log2-fold-change associated with a specific constrast/coef using 'apeglm' in lfcShrink() for visualization and ranking of genes [Reference: 16, 21]
#	Why is it required? "The moderated log fold changes proposed by Love, Huber, and Anders (2014) use a normal prior distribution, centered on zero and with a scale that is fit to the data. The shrunken log fold changes are useful for ranking and visualization, without the need for arbitrary filters on low count genes. The normal prior can sometimes produce too strong of shrinkage for certain datasets. In DESeq2 version 1.18, we include two additional adaptive shrinkage estimators, available via the type argument of lfcShrink. [Reference: 16]. apeglm stabds for "Approximate Posterior Estimation for GLM”".
#
#	What is s-value, and how to interpret it? "The adjusted p-values and s-values are similar but with a different definition of error. One focuses on falsely rejecting what are truly null genes, and the other on getting the sign of the LFC wrong. You can use whichever you prefer or feel more comfortable with. Because we were computing posterior distributions in apeglm, and because we were adding ashr at the same time to lfcShrink, we decided it made sense to provide s-values as output. I think there is one clear benefit to outputting s-values, which is during method development, we can assess and benchmark power and error control simultaneously with real data. It's very difficult to find real datasets in which there are non-null genes and null genes, if we have a point null hypothesis of LFC=0". [Reference: 22]
#	
#	What is LRT? 
#	"LRT is used to identify any genes that show change in expression across the different levels. This type of test can be especially useful in analyzing time course experiments." [Reference: 24]
#	"Likelihood ratio test (chi-squared test) for GLMs: This function tests for significance of change in deviance between a full and reduced model which are provided as ‘formula’. Fitting uses previously calculated ‘sizeFactors’ (or ‘normalizationFactors’) and dispersion estimates. The difference in deviance is compared to a chi-squared distribution with df = (reduced residual degrees of freedom - full residual degrees of freedom)." [?nbinomLRT]
#
#	How to interpret p-value in LRT? 
#	"The p-values are determined solely by the difference in deviance between the ‘full’ and ‘reduced’ model formula (not log2 fold changes). Essentially the LRT test is testing whether the term(s) removed in the ‘reduced’ model explains a significant amount of variation in the data?" [Reference: 24]
#	"even though there are fold changes present they are not directly associated with the actual hypothesis test. Thus, when filtering significant genes from the LRT we use only the FDR as our threshold. How many genes are significant at padj < 0.05?. We are unable to set a fold change criteria here since the statistic is not generated from any one pairwise comparison. This list includes genes that can be changing in any number of combinations across the three factor levels. It is advisable to instead increase the stringency on our criteria and lower the FDR threshold." [Reference: 24]
#	"The likelihood ratio test p values therefore represent a test of all the variables and all the levels of factors which are among these variables. However, the results table only has space for one column of log fold change, so a single variable and a single comparison is shown (among the potentially multiple log fold changes which were tested in the likelihood ratio test). This is indicated at the top of the results table with the text, e.g., log2 fold change (MLE): condition C vs A, followed by, LRT p-value: ‘~ batch + condition’ vs ‘~ batch’. This indicates that the p value is for the likelihood ratio test of all the variables and all the levels, while the log fold change is a single comparison from among those variables and levels." [Reference: 23]
#
#
#	How to analysis time series data: There are multiple ways depending on biological question of interest. [Reference: 25]
#	1. In order to test for any differences (in expression of a gene) over multiple time points, one can use a design including the time factor (e.g., ~ condition + time), and then test using the likelihood ratio test where the time factor is removed in the reduced formula (e.g., ~ condition). 
#	2. For a control and treatment time series, one can use a design formula containing the condition factor, the time factor, and the interaction of the two (e.g., ~ condition + time + condition:time). In this case, using the likelihood ratio test with a reduced model which does not contain the interaction terms (e.g., ~ condition + time) will test whether the condition induces a change in gene expression at any time point after the reference level time point (time 0).

#	On interaction terms:
#	Currently, the design for full and redcued models are internally created based on --primary_factor_column and --other_known_factors_columns. The script does not support interaction terms in design formula for full or reduced model. Therefore, users are requested to follow a much simpler approach to achieve their goal which includes combining the factors of interest into a single factor with all combinations of the original factors, and then changing the design to include just this factor, e.g., ~ group. Please refer this: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions [Reference: 25]
#	Update as in version 1.3 on 18th March 2022: user can directly specify --full_design and --reduced_design with interaction terms.
#
#	MA plot: [Reference: 17]
#		It helps to visualize log-fold-change between experimental groups (M) against the mean expression across all the samples (A) for each gene.
#		The Glimma MA plot has two components:
#			1) A plot of summary statistics across all genes that have been tested	[nn Left Panel]
#			2) A plot of gene expression from individual samples for a given gene	[on Right Panel]
#		Below the two plots, we have a table wherein user can apply a filter
#		Note: transform.counts = "logcpm" in Glimma::glimmaMA corresponds to edgeR::cpm(counts, log=TRUE)
#
#	Volcano plot: [Reference 18]
#		It helps to visualize statistical significance (on Y-axis) against log-fold-change (on X-axis) between experimental groups for each gene.
#		It has two components:
#			1) A plot of gene-wise log-fold-change on x-axis versus -log10(svalue)	[on Left Panel]
#			2) A plot of gene expression from individual samples for a given gene	[on Right Panel]
#		Below the two plots, we have a table wherein user can apply a filter
#

##	Identification of gene clusters exhibiting particular patterns across samples [Reference: 24]
#	"Often we are interested in genes that have particular patterns across the sample groups (levels) of our condition. For example, with the MOV10 dataset, we may be interested in genes that exhibit the lowest expression for the Mov10_KD and highest expression for the Mov10_OE sample groups (i.e. KD < CTL < OE). To identify genes associated with these patterns we can use a clustering tool, degPatterns from the ‘DEGreport’ package, that groups the genes based on their changes in expression across sample groups." 
#	"The function 'degPatterns' determines sets of genes that exhibit similar expression patterns across sample groups."
#	"The degPatterns tool uses a hierarchical clustering approach based on pair-wise correlations, then cuts the hierarchical tree to generate groups of genes with similar expression profiles. The tool cuts the tree in a way to optimize the diversity of the clusters, such that the variability inter-cluster > the variability intra-cluster."
#	There is some error (Error: $ operator is invalid for atomic vectors) working with this function. Therefore, had to write custom function to achieve similar result following instructions here in Reference: 27, 29, 32. The instruction mentioned a step where correlation was calculated using cor.test() function. However, cor.test() on all pairwise combinations of a long list of genes seems very time-consuming. Therefore, I have used cor() function. Additionally, as an alternative approach a soft clustering using fuzzy c-means algorithm as part of Mfuzz R package can be used [Reference: 28].
#	As a metric for goodness of a clustering technique, 'Silhouette Coefficient' can be used. Its value ranges from -1 to 1, with 1 means clusters are well apart from each other and clearly distinguished, 0 means clusters are indifferent, or we can say that distance between clusters is not significant, and -1 means clustered are assigned in the wrong way [Reference: 30, 31]. 


## References:
# 1. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/
# 2. https://academic.oup.com/bioinformatics/article/35/12/2084/5159452
# 3. https://doi.org/10.1093/bioinformatics/bty895
# 4. https://www.biostars.org/p/356032/
# 5. https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#removing-genes-that-are-lowly-expressed
# 6. https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# 7. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411036/
# 8. https://www.nature.com/articles/s41598-020-74567-y#Sec5
# 9.https://doi.org/10.1016/j.jmoldx.2019.04.007
# 10. https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
# 11. https://www.frontiersin.org/articles/10.3389/fgene.2020.00594/full
# 12. https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-021-02936-w#Sec2
# 13. https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#transformations-from-the-raw-scale
# 14. https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#count-data-transformations
# 15. https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html
# 16. http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#methods-changes-since-the-2014-deseq2-paper
# 17. Glimma: http://bioconductor.org/packages/release/bioc/vignettes/Glimma/inst/doc/DESeq2.html
# 18. Glimma: https://bioconductor.org/packages/release/bioc/manuals/Glimma/man/Glimma.pdf
# 19. Density plot: https://www.statology.org/density-curves/
# 20. clustering on gene expression data: https://online.stat.psu.edu/stat555/node/85/
# 21. http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#moreshrink
# 22. https://support.bioconductor.org/p/113664/
# 23. LRT: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#likelihood-ratio-test
# 24. LRT: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
# 25. Interaction terms: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions
# 26. http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# 27. degPatterns: https://rdrr.io/github/lpantano/DEGreport/man/degPatterns.html
# 28. Mfuzz: https://www.bioconductor.org/packages/release/bioc/vignettes/Mfuzz/inst/doc/Mfuzz.pdf
# 29: diana: https://www.rdocumentation.org/packages/cluster/versions/2.1.2/topics/diana
# 30. Silhouette Coefficient: https://towardsdatascience.com/silhouette-coefficient-validating-clustering-techniques-e976bb81d10c
# 31. Silhouette Coefficient: https://medium.com/codesmart/r-series-k-means-clustering-silhouette-794774b46586
# 32. hierarchical clustering: agglomerative vs divisive: https://bradleyboehmke.github.io/HOML/hierarchical.html
# 33. 

##########################################

cat("For help: /path/to/Rscript --vanilla /path/to/Differential_gene_expression_analysis_on_counts.R --help", "\n", sep = "");


## =========================================================
## Take user specified parameters
## =========================================================
library("optparse");	# command line option parser

args_list <- list(
	optparse::make_option(opt_str = c("--DGE_tool"), type = "character", default = "DESeq2", help = "Supported tools are: DESeq2. Default: DESeq2"),
    optparse::make_option(opt_str = c("--count_file"), type = "character", default = NULL, help = "Feature (rows) x Sample (columns) expression file in CSV or TSV format. Default: NULL"),
	optparse::make_option(opt_str = c("--gene_id_column"), type = "integer", default = 1, help = "Column number in --counts_file which contains gene identifiers. Default: 1"),
    optparse::make_option(opt_str = c("--count_column"), type = "integer", default = 3, help = "Starting column number containing the expression data values. Default: 3"),
	optparse::make_option(opt_str = c("--counts_threshold"), type = "integer", default = 10, help = "A raw count value used for pre-filtering genes; If --use_CPM_cutoff is set FALSE, then pre-filtering will be done only using --counts_threshold. In that case, only those genes which have total counts across all samples more than this value will get selected for downstream analysis. If --use_CPM_cutoff is set TRUE, a gene will considered to be expressed at biologically meaningful level in a sample if its abundance in Counts-Per-Million (CPM) unit is greater than the CPM value calculated on this raw count for the least sequenced sample. For example, if --counts_threshold is 10 and library size of least sequenced sample is 10,000,000, then biological meaningful expression level will be (10/10000000)*1000000 = 1 CPM. In such case, --sample_coverage_threshold will also be considered, e.g., if --sample_coverage_threshold is set to 0.5, then all genes with abundance value of at least 1 CPM in at least 50% samples will be considered for downstream analysis. Default: 10"),
	optparse::make_option(opt_str = c("--use_CPM_cutoff"), type = "logical", default = TRUE, help = "Whether a cutoff calculated on --counts_threshold and library size of the least sequenced sample will be used for pre-filtering on counts data. Already described above. Default: TRUE"),
	optparse::make_option(opt_str = c("--sample_coverage_threshold"), type = "double", default = 0.5, help = "A feature will be discarded if not abundant at biologically meaningful level in at least this fraction of samples. This is applicable only when --use_CPM_cutoff is set TRUE. Range: 0.0 - 1.0; Default: 0.5"),
	optparse::make_option(opt_str = c("--normalization_method"), type = "character", default = 'DESeq2-median-of-ratios', help = "Supported methods are:  DESeq2-median-of-ratios. Default: DESeq2-median-of-ratios"),
	optparse::make_option(opt_str = c("--pseudocounts"), type = "double", default = 0, help = "It is a positive constant to avoid taking log of 0 or any negative values. Useful across data filter, normalization and transformation steps. Default: 0"),
	optparse::make_option(opt_str = c("--sample_info_file"), type = "character", default = NULL, help = "A sample metadata file (in TSV/CSV) with list of samples of interest. Default: NULL"),
	optparse::make_option(opt_str = c("--sample_id_column"), type = "integer", default = 1, help = "Column number for sample ids. Default: 1"),
	optparse::make_option(opt_str = c("--primary_factor_column"), type = "integer", default = 5, help = "Column number for the primary factor of interest such as sample conditions/groups. Default: 5"),
	optparse::make_option(opt_str = c("--contrasts"), type = "character", default = "all", help = "Comma separated list of contrasts for DGE in the defined format: cond2_vs_cond1 where cond1 and cond2 should be values from primary_factor_column. Applicable only for pairwise DGE analysis such as Wald test from DESeq2. Default: all, means all pairwise combinations will be considered"),
	optparse::make_option(opt_str = c("--other_known_factors_columns"), type = "character", default = NULL, help = "Comma separated list of column numbers of sample_info_file which can act as other possible source of variation in gene expression. Default: NULL"),
	optparse::make_option(opt_str = c("--full_design"), type = "character", default = NULL, help = "Design formula for full model. If defined, it will get preference over --primary_factor_column and --other_known_factors_columns. If not provided, it will be defined based on both parameters, e.g., '~ genotype + condition'. Interaction terms are supported, so user can specify something like '~ condition + time + condition:time'. Default: NULL"),
	optparse::make_option(opt_str = c("--reduced_design"), type = "character", default = NULL, help = "Design formula for reduced model. If defined, it will get preference over --other_known_factors_columns. If not provided, it will be defined based on --other_known_factors_columns, e.g., '~ genotype'. Default: NULL"),
	optparse::make_option(opt_str = c("--DESeq2_dispersion_method"),	type = "character",	default = "parametric",					help = "Method to estimate dispersion such as parametric/local/mean/glmGamPoi. 'glmGamPoi' dispersion estimator should be used in combination with a LRT and not a Wald test. Default: parametric."),
	optparse::make_option(opt_str = c("--DESeq2_test_type"), type = "character", default = "Wald", help = "Type of statistical test such as Wald/LRT. Default: Wald"),
	optparse::make_option(opt_str = c("--stat_sig_measure"), type = "character", default = "padj", help = "Type of measurement for statistical significance level, pvalue, padj or svalue. svalue is only applicable for --DESeq2_test_type = Wald. Default: padj"),
	optparse::make_option(opt_str = c("--stat_cutoff"),	type = "character",	default = "0.05", help = "Cutoff for calling differentially expressed genes. Comma-separated if multiple. Default: 0.05"),
	optparse::make_option(opt_str = c("--log2_fold_change_cutoff"),	type = "character",	default = "1", help = "log2-fold-change cutoff for calling differentially expressed genes. Applicable only for DESeq2_test_type = Wald. Comma-separated values if multiple. Default: 1"),
	optparse::make_option(opt_str = c("--top_n_for_viz"), type = "character", default = NULL, help = "n most variable genes will be considered for generating following plots: clustered heatmaps, MDS. Comma-separated values if multiple. Default: NULL, means all genes will be considered."),
	optparse::make_option(opt_str = c("--output_dir"), type = "character", default = "./", help = "Output directory path. Default: ./"),
	optparse::make_option(opt_str = c("--output_prefix"), type = "character", default = NULL, help = "Prefix for output files. Default: NULL"),
	optparse::make_option(opt_str = c("--log_file"), type = "character", default = "DGE.log", help = "Log file. Default: DGE.log")
)
arguments <- optparse::parse_args(optparse::OptionParser(option_list = args_list));


## =========================================================
## Load required libraries
## =========================================================
library("typed");		# to set types for a variable, and for arguments of a function
library("log4r");		# logging library
library("data.table");	# read input data fast
library("utils");		# for combination calculation
library("DESeq2");		# normalize count data by median-of-ratios or RLE; apply VST/rlog transformation; perform DGE analysis
library("edgeR");		# normalize count data by TMM; apply CPM transformation
library("apeglm");		# for effect size shrinkage
library("glmGamPoi");	# for faster dispersion and parameter estimation in LRT, particularly useful for single-cell data
library("ggplot2");		# for drawing plots
library("ggsci");		# to support Nature Publishing Group (NPG) color palettes
library("gplots");		# to generate smoothly varying set of colors
library("plotly");		# for interactive visualization
library("htmlwidgets");	# for saving html widgets
library("SummarizedExperiment");	# for extracting transformed matrix from DESeqTransform object
library("vsn");			# for generating mean-sd plot
library("hexbin");		# for generating bins during mean-sd plot
library("ComplexHeatmap");	# for generating static cluster heatmaps
library("magrittr");	# brings %>% operator
library("iheatmapr");	# for generating interactve cluster heatmaps
library("Glimma");		# for generating interactive MDS and MA plots
library("DEGreport");	# identify gene clusters exhibiting particular patterns across samples (based on DESeq2's LRT output); uses hard partition for clustering.
library("ggpubr");		# combine multiple ggplot on one page
#library("Mfuzz");		# identify gene clusters exhibiting particular patterns across samples (based on DESeq2's LRT output); uses soft partition for clustering.

# To avoid following error, Error in .External2(C_X11, paste0("png::", filename), g$width, g$height,
options(bitmapType = 'cairo');

# set highcharter options
#options(highcharter.theme = highcharter::hc_theme_smpl(tooltip = list(valueDecimals = 2)));


## =========================================================
## Set imp. variables
## =========================================================
declare("DGE_tool", Character(length = 1, null_ok = FALSE), value = arguments$DGE_tool, const = TRUE);

declare("count_file", Character(length = 1, null_ok = TRUE), value = arguments$count_file, const = TRUE);
declare("gene_id_column", Integer(length = 1, null_ok = FALSE), value = arguments$gene_id_column, const = TRUE);
declare("count_column", Integer(length = 1, null_ok = FALSE), value = arguments$count_column, const = TRUE);
declare("counts_threshold", Integer(length = 1, null_ok = FALSE), value = arguments$counts_threshold, const = TRUE);
declare("use_CPM_cutoff", Logical(length = 1, null_ok = TRUE), value = arguments$use_CPM_cutoff, const = TRUE);
declare("sample_coverage_threshold", Double(length = 1, null_ok = FALSE), value = arguments$sample_coverage_threshold, const = TRUE);
declare("normalization_method", Character(length = 1, null_ok = FALSE), value = arguments$normalization_method, const = TRUE);
declare("pseudocounts", Double(length = 1, null_ok = FALSE), value = arguments$pseudocounts, const = TRUE);

declare("sample_info_file", Character(length = 1, null_ok = TRUE), value = arguments$sample_info_file, const = TRUE);
declare("sample_id_column", Integer(length = 1, null_ok = FALSE), value = arguments$sample_id_column, const = TRUE);
declare("primary_factor_column", Integer(length = 1, null_ok = FALSE), value = arguments$primary_factor_column, const = TRUE);

if(!is.null(arguments$contrasts))
{
	arguments$contrasts <- gsub(x = arguments$contrasts, pattern = "\\s", replacement = "");
}
declare("contrasts", Character(null_ok = TRUE), value = arguments$contrasts, const = FALSE);

if(!is.null(arguments$other_known_factors_columns))
{
	arguments$other_known_factors_columns <- gsub(x = arguments$other_known_factors_columns, pattern = "\\s", replacement = "");
	if(arguments$other_known_factors_columns == 'NULL')
	{
		arguments$other_known_factors_columns <- NULL;
	}
}
declare("other_known_factors_columns", Integer(null_ok = TRUE), const = FALSE);
if(!is.null(arguments$other_known_factors_columns))
{
	other_known_factors_columns <- as.integer( unlist(strsplit(arguments$other_known_factors_columns, split = ",")) );
	if(any(is.na(other_known_factors_columns)))
	{
		cat("Error: --other_known_factors_columns does not have all integer values.", "\n", sep = "");
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed.", "\n", sep = "");
		q(save = "no", status = 1);
	}
}

if(!is.null(arguments$full_design) && (arguments$full_design == 'NULL'))
{
	arguments$full_design <- NULL;
}else if(!is.null(arguments$full_design) && (arguments$full_design != 'NULL'))
{
	if(!grepl(arguments$full_design, pattern = "^~"))
	{
		cat("Error: If defined, --full_design should be in following format: ~ gentotype + condition1 + condition2 + .. where all factors are appropriate columns from --sample_info_file", "\n", sep = "");
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed.", "\n", sep = "");
		q(save = "no", status = 1);
	}
}
declare("full_design", Any(), value = arguments$full_design, const = FALSE);

if(!is.null(arguments$reduced_design) && (arguments$reduced_design == 'NULL'))
{
	arguments$reduced_design <- NULL;
}else if(!is.null(arguments$reduced_design) && (arguments$reduced_design != 'NULL'))
{
	if(!grepl(arguments$reduced_design, pattern = "^~"))
	{
		cat("Error: If defined, --reduced_design should be in following format: ~ 1 or ~ condition1 + condition2 where all factors are appropriate columns from --sample_info_file", "\n", sep = "");
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed.", "\n", sep = "");
		q(save = "no", status = 1);
	}
}
declare("reduced_design", Any(), value = arguments$reduced_design, const = FALSE);

declare("DESeq2_dispersion_method", Character(length = 1, null_ok = FALSE), value = arguments$DESeq2_dispersion_method, const = FALSE);
declare("DESeq2_test_type", Character(length = 1, null_ok = FALSE), value = arguments$DESeq2_test_type, const = TRUE);
declare("stat_sig_measure", Character(length = 1, null_ok = FALSE), value = arguments$stat_sig_measure, const = FALSE);

arguments$stat_cutoff <- gsub(x = arguments$stat_cutoff, pattern = "\\s", replacement = "");
if(arguments$stat_cutoff == 'NULL')
{
	cat("Error: --stat_cutoff can not be a NULL.", "\n", sep = "");
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed.", "\n", sep = "");
	q(save = "no", status = 1);
}
declare("stat_cutoff", Double(null_ok = FALSE), value = as.double(unlist(strsplit(arguments$stat_cutoff, split = ","))), const = TRUE);
if(any(is.na(stat_cutoff)))
{
	cat("Error: --stat_cutoff does not have all double values.", "\n", sep = "");
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed.", "\n", sep = "");
	q(save = "no", status = 1);
}

arguments$log2_fold_change_cutoff <- gsub(x = arguments$log2_fold_change_cutoff, pattern = "\\s", replacement = "");
if(arguments$log2_fold_change_cutoff == 'NULL')
{
	cat("Error: --log2_fold_change_cutoff can not be a NULL.", "\n", sep = "");
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed.", "\n", sep = "");
	q(save = "no", status = 1);
}
declare("log2_fold_change_cutoff", Double(null_ok = FALSE), value = as.double(unlist(strsplit(arguments$log2_fold_change_cutoff, split = ","))), const = TRUE);
if(any(is.na(log2_fold_change_cutoff)))
{
	cat("Error: --log2_fold_change_cutoff does not have all double values.", "\n", sep = "");
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed.", "\n", sep = "");
	q(save = "no", status = 1);
}

if(!is.null(arguments$top_n_for_viz))
{
	arguments$top_n_for_viz <- gsub(x = arguments$top_n_for_viz, pattern = "\\s", replacement = "");
	if(arguments$top_n_for_viz == 'NULL')
	{
		arguments$top_n_for_viz <- NULL;
	}
}
declare("top_n_for_viz", Integer(null_ok = TRUE), const = FALSE);
if(!is.null(arguments$top_n_for_viz))
{
	top_n_for_viz <- as.integer( unlist(strsplit(arguments$top_n_for_viz, split = ",")) );
	if(any(is.na(top_n_for_viz)))
	{
		cat("Error: --top_n_for_viz does not have all integer values.", "\n", sep = "");
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed.", "\n", sep = "");
		q(save = "no", status = 1);
	}
}

declare("output_dir", Character(length = 1, null_ok = FALSE), value = arguments$output_dir, const = TRUE);

if(!is.null(arguments$output_prefix) && (arguments$output_prefix == 'NULL'))
{
	arguments$output_prefix <- NULL;
}
declare("output_prefix", Character(length = 1, null_ok = TRUE), value = arguments$output_prefix, const = FALSE);
output_prefix <- ifelse(test = is.null(output_prefix), yes = "", no = paste(output_prefix, "_", sep = ""));

declare("log_file", Character(length = 1, null_ok = FALSE), value = arguments$log_file, const = FALSE);


declare("supported_DGE_tools", Character(), value = c("DESeq2"), const = TRUE);
declare("supported_normalization_methods", Character(), value = c("DESeq2-median-of-ratios"), const = TRUE);

declare("expression_mat", Matrix(null_ok = TRUE), const = FALSE);				# count data - matrix
declare("conditions", Character(), const = FALSE);								# sample conditions - vector
declare("condition_col", Character(), const = FALSE);							# labels for sample conditions - vector
declare("dds", Any(), const = FALSE);											# DESeqDataSet object
declare("normalized_transformed_mat", Matrix(null_ok = TRUE), const = FALSE);	# normalized and transformed data - matrix 
declare("DGE_result", Data.frame());											# DGE results for all comparisions - data frame
declare("effect_size_measure", Character(length = 1, null_ok = TRUE));			# Measure for effect size; log2FoldChange/log2FoldChange_apeglm/NULL based on DESeq2_test_type and stat_sig_measure


## Create output directory ----------------------------------------
if(!dir.exists(output_dir))
{
	dir.create(path = output_dir, recursive = TRUE, mode = "0777");
}

declare("log_dir", Character(length = 1), value = paste(output_dir, "/", "log", sep = ""), const = FALSE);
if(!dir.exists(log_dir))
{
	dir.create(path = log_dir, recursive = TRUE, mode = "0777");
}

## Set file logger -------------------------------------------------
declare("log_file", Character(length = 1), value = paste(log_dir, "/", log_file, sep = ""), const = TRUE);
declare("file_logger", Any(), value = log4r::logger(threshold = "DEBUG", appenders = log4r::file_appender(file = log_file, append = TRUE)), const = TRUE);
info(file_logger, "Processing starts");
# Additionally, for capturing stdouts 
base::sink(file = log_file, append = TRUE);

## Record user selected parameters ---------------------------------
info(file_logger, "Parameters are: ");
info(file_logger, paste(head(x = names(arguments), n = -1), head(x = as.character(arguments), n = -1), sep = ":", collapse = ", "));
info(file_logger, "Note: Above, all parameters of the scripts have been shown including those which were default and not set by the User. Few might not be relevant as per the options chosen by the User, and those will be internally ignored.");

## Check validity of some parameters: DGE_tool, count_file, sample_info_file, normalization_method, DESeq2_test_type, stat_sig_measure, contrasts -------
if(!(DGE_tool %in% supported_DGE_tools))
{
	fatal(file_logger, paste("A non-supported value for --DGE_tool: ", DGE_tool, " has been selected. Please choose one of the following: DESeq2", "\n", "For more help: /path/to/Rscript --vanilla /path/to/Differential_gene_expression_analysis_on_counts.R --help", sep = ""));
	base::sink();
	cat("A non-supported value for --DGE_tool: ", DGE_tool, " has been selected. Please choose one of the following: DESeq2", "\n", sep = "");
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed. Check log for more details: ", log_file, "\n", sep = "");
	q(save = "no", status = 1);
}else if(is.null(count_file))
{
	fatal(file_logger, paste("--count_file should have a valid file. User has supplied NULL.", "\n", "For more help: /path/to/Rscript --vanilla /path/to/Differential_gene_expression_analysis_on_counts.R --help", sep = ""));
	base::sink();
	cat("--count_file should have a valid file. User has supplied NULL.", "\n", sep = "")
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed. Check log for more details: ", log_file, "\n", sep = "");
	q(save = "no", status = 1);
}else if(is.null(sample_info_file))
{
	fatal(file_logger, paste("--sample_info_file should have a valid file. User has supplied NULL.", "\n", "For more help: /path/to/Rscript --vanilla /path/to/Differential_gene_expression_analysis_on_counts.R --help", sep = ""));
	base::sink();
	cat("--sample_info_file should have a valid file. User has supplied NULL.", "\n", sep = "")
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed. Check log for more details: ", log_file, "\n", sep = "");
	q(save = "no", status = 1);
}else if(!(normalization_method %in% supported_normalization_methods))
{
	fatal(file_logger, paste("A non-supported value for --normalization method: ", normalization_method," has been selected. Please choose one of the following: ", paste(supported_normalization_methods, collapse = ","), "\n", "For more help: /path/to/Rscript --vanilla /path/to/Differential_gene_expression_analysis_on_counts.R --help", sep = ""));
	base::sink();
	cat("A non-supported value for --normalization method: ", normalization_method," has been selected. Please choose one of the following: ", paste(supported_normalization_methods, collapse = ","), "\n", sep = "")
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed. Check log for more details: ", log_file, "\n", sep = "");
	q(save = "no", status = 1);
}else if(normalization_method == "DESeq2-median-of-ratios" && DGE_tool != "DESeq2")
{
	fatal(file_logger, paste("Currently we only support normalization method - ", normalization_method, " with DESeq2", ", but you have selected ", DGE_tool, sep = ""));
	base::sink();
	cat("Currently we only support normalization method - ", normalization_method, " with DESeq2", ", but you have selected ", DGE_tool, "\n", sep = "");
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed. Check log for more details: ", log_file, "\n", sep = "");
	q(save = "no", status = 1);
}else if(!(DESeq2_dispersion_method %in% c("parametric", "local", "mean", "glmGamPoi")))
{
	fatal(file_logger, paste("A non-supported value for --dispersion_method: ", DESeq2_dispersion_method, " has been selected. Please choose one of the following: parametric,local,mean,glmGamPoi", "\n", "For more help: /path/to/Rscript --vanilla /path/to/Differential_gene_expression_analysis_on_counts.R --help", sep = ""));
	base::sink();
	cat("A non-supported value for --dispersion_method: ", DESeq2_dispersion_method, " has been selected. Please choose one of the following: parametric,local,mean,glmGamPoi", "\n", sep = "");
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed. Check log for more details: ", log_file, "\n", sep = "");
	q(save = "no", status = 1);
}else if(!(DESeq2_test_type %in% c("Wald", "LRT")))
{
	fatal(file_logger, paste("A non-supported value for --DESeq2_test_type: ", DESeq2_test_type, " has been selected. Please choose one of the following: Wald, LRT", "\n", "For more help: /path/to/Rscript --vanilla /path/to/Differential_gene_expression_analysis_on_counts.R --help", sep = ""));
	base::sink();
	cat("A non-supported value for --DESeq2_test_type: ", DESeq2_test_type, " has been selected. Please choose one of the following: Wald, LRT", "\n", sep = "");
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed. Check log for more details: ", log_file, "\n", sep = "");
	q(save = "no", status = 1);
}else if(DESeq2_dispersion_method == "glmGamPoi" && DESeq2_test_type != "LRT")
{
	fatal(file_logger, "--DESeq2_dispersion_method glmGamPoi should be used in combination with --DESeq2_test_type LRT");
	base::sink();
	cat("--DESeq2_dispersion_method = glmGamPoi should be used in combination with --DESeq2_test_type = LRT", "\n", sep = "");
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed. Check log for more details: ", log_file, "\n", sep = "");
	q(save = "no", status = 1);
}else if(!(stat_sig_measure %in% c("pvalue", "padj", "svalue")))
{
	fatal(file_logger, paste("A non-supported value for --stat_sig_measure: ", stat_sig_measure, " has been provided. Please choose one of the following: pvalue, padj, svalue", "\n", "For more help: /path/to/Rscript --vanilla /path/to/Differential_gene_expression_analysis_on_counts.R --help", sep = ""));
	base::sink();
	cat("A non-supported value for --stat_sig_measure: ", stat_sig_measure, " has been provided. Please choose one of the following: pvalue, padj, svalue", "\n", sep = "");
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed. Check log for more details: ", log_file, "\n", sep = "");
	q(save = "no", status = 1);
}

if( DESeq2_test_type == "Wald" && (is.null(contrasts) || contrasts == 'NULL') )
{
	warn(file_logger, "Since, --DESeq2_test_type is set to Wald, --contrasts can not be NULL. Therefore, all pairwise combinations of the primary factor of interest (which is based on --primary_factor_column) will be considered.");
	contrasts <- "all";
}

if( DESeq2_test_type == "LRT" && stat_sig_measure != "padj")
{
	warn(file_logger, "Since LRT has been selected for --DESeq2_test_type, --stat_sig_measure will be set to 'padj'");
	stat_sig_measure <- "padj";
}

## Set effect_size_measure based on stat_sig_measure -------------------------
if(DESeq2_test_type == "Wald" && stat_sig_measure %in% c("pvalue", "padj"))
{
	effect_size_measure <- "log2FoldChange";
	info(file_logger, paste("Since --stat_sig_measure was set to pvalue/padj, we are going to use log2-fold-change estimated by differential expression tool as effect size measure. The column 'log2FoldChange' from the main result will be used for various plots inluding MA, Volcano etc.", "\n", sep = ""));
}else if(DESeq2_test_type == "Wald" && stat_sig_measure == "svalue")
{
	effect_size_measure <- "log2FoldChange_apeglm";
	info(file_logger, paste("Since --stat_sig_measure was set to svalue, we are going to use log2-fold-change as estimated by 'apeglm' R package as effect size measure. The column 'log2FoldChange_apeglm' from the main result will be used for various plots including MA, Volcano etc.", "\n", sep = ""));
}else if(DESeq2_test_type == "LRT")
{
	info(file_logger, paste("Since --DESeq2_test_type was set to LRT, we are not going to use an effect size measure such as log2-fold-change in the downstream analysis.", "\n", sep = ""));
	effect_size_measure <- NULL;
}else
{
	fatal(file_logger, paste("Unable to assign a value to effect_size_measure based on --DESeq2_test_type and --stat_sig_measure.", "\n", sep = ""));
	base::sink();
	cat("Unable to assign a value to effect_size_measure based on --DESeq2_test_type and --stat_sig_measure.", "\n", sep = "");
	cat("The program: Differential_gene_expression_analysis_on_counts.R failed. Check log for more details: ", log_file, "\n", sep = "");
	q(save = "no", status = 1);
}




## =========================================================
## Custom functions
## =========================================================
## Define some color palettes -----------------------------
# Reference: https://www.r-bloggers.com/the-paul-tol-21-color-salute/
declare("RBGO4", Character(null_ok = FALSE), value = c("red", "blue", "darkgreen", "orange"), const = TRUE);
declare("tol21rainbow", Character(null_ok = FALSE), value = c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788"), const = TRUE);
declare("NPG", Character(null_ok = FALSE), value = ggsci::pal_npg(palette = "nrc", alpha = 1)(10), const = TRUE);

## Specify colors for different factors: conditions, batch, stage codes etc -----------
getNamedColorVectors <- Character() ? function(category_names = Character())
{
	declare("n", Integer(), value = length(category_names));
	declare("category_cols", Character(null_ok = TRUE));

	#if(n <= 4)
	#{
	#	category_cols <- RBGO4[1:n];
	#}else 
	if(n <= 10)
	{
		category_cols <- NPG[1:n];
	}else if( n > 10 && n <= 21 )
	{
		category_cols <- tol21rainbow[1:n]
	}else
	{
		category_cols <- gplots::colorpanel(n = n, low = "darkviolet" , mid = "darkgreen" , high = "darkorange")
	};
	
	names(category_cols) <- as.character(category_names);
	return(category_cols);
}

## Draw density plot ---------------------------------------
drawDensityPlot <- Logical() ? function(data_df = Data.frame(), x = Character(), group_by = Character(), fill_by = Character(), manual_color_map = Character(), title = Character(), x_axis_title = Character(), y_axis_title = Character(), outfile = Character())
{
	declare("p", Any());
	p <- ggplot2::ggplot(data = data_df) + theme_bw() + ggtitle(title) + xlab(x_axis_title) + ylab(y_axis_title);
	p <- p + geom_density(mapping = aes_string(x = base::as.name(x), group = base::as.name(group_by), fill = base::as.name(fill_by)), alpha = 0.4);
	p <- p + scale_fill_manual(values = manual_color_map);
	
	declare("outfile_png", Character(length = 1), value = sub(x = outfile, pattern = "\\.html", replacement = ".png"));
	png(filename = outfile_png, units = "in", height = 8, width = 10, res = 300, type = "cairo");
	print(p);
	dev.off();
	invisible(gc());

	#debug(file_logger, "Temporarily disabling generation of html file.");
	declare("fig", Any(), value = plotly::ggplotly(p));
	htmlwidgets::saveWidget(widget = fig, file = outfile, selfcontained = TRUE);

	#debug(file_logger, paste("Density plots are stored as ", outfile_png, " and ", outfile, sep = ""));
	return(TRUE);
}

## Draw bar plot -------------------------------------------
drawBarPlot <- Logical() ? function(data_df = Data.frame(), x = Character(), y = Character(), fill_by = Character(), manual_color_map = Character(), title = Character(), x_axis_title = Character(), y_axis_title = Character(), outfile = Character())
{
	declare("p", Any());
	p <- ggplot2::ggplot(data = data_df) + theme_bw() + ggtitle(title) + xlab(x_axis_title) + ylab(y_axis_title);
	p <- p + geom_bar(mapping = aes_string(x = base::as.name(x), y = base::as.name(y), fill = base::as.name(fill_by)), stat  = "identity", show.legend = TRUE);
	p <- p + scale_fill_manual(values = manual_color_map) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1.0));
	
	declare("outfile_png", Character(length = 1), value = sub(x = outfile, pattern = "\\.html", replacement = ".png"));
	png(filename = outfile_png, units = "in", height = 8, width = 10, res = 300, type = "cairo");
	print(p);
	dev.off();
	invisible(gc());

	#debug(file_logger, "Temporarily disabling generation of html file.");
	declare("fig", Any(), value = plotly::ggplotly(p));
	htmlwidgets::saveWidget(widget = fig, file = outfile, selfcontained = TRUE);

	#debug(file_logger, paste("Bar plots are stored as ", outfile_png, " and ", outfile, sep = ""));
	
	return(TRUE);
}

## Draw box plot -------------------------------------------
drawBoxPlot <- Logical() ? function(data_df = Data.frame(), x = Character(), y = Character(), fill_by = Character(), manual_color_map = Character(), title = Character(), x_axis_title = Character(), y_axis_title = Character(), outfile = Character())
{
	declare("p", Any());
	p <- ggplot2::ggplot(data = data_df, mapping = aes_string(x = base::as.name(x), y = base::as.name(y), fill = base::as.name(fill_by))) + theme_bw() + ggtitle(title) + xlab(x_axis_title) + ylab(y_axis_title);
	p <- p + geom_boxplot();
	p <- p + scale_fill_manual(values = manual_color_map) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1.0));

	declare("outfile_png", Character(length = 1), value = sub(x = outfile, pattern = "\\.html", replacement = ".png"));
	png(filename = outfile_png, units = "in", height = 8, width = 10, res = 300, type = "cairo");
	print(p);
	dev.off();

	declare("fig", Any(), value = plotly::ggplotly(p));
	htmlwidgets::saveWidget(widget = fig, file = outfile, selfcontained = TRUE);

	return(TRUE)
}

## Filter features based on their variance across samples ---
filterByVariance <- Character(null_ok = FALSE) ? function(data_mat = Matrix(), top_n = Integer())
{
	declare("var_vec", Double(null_ok = TRUE));
	declare("topn_var_vec", Character(null_ok = TRUE));

	var_vec <- rep(0, nrow(data_mat));
	var_vec <- unlist(lapply(X = rownames(data_mat), FUN = function(x) {var(data_mat[x,])}));
	names(var_vec) <- rownames(data_mat);
	invisible(gc());
	topn_var_vec <- names(head(x = sort(var_vec, decreasing = TRUE), n = top_n));
		
	return(topn_var_vec);
}

## Draw Complex Heatmap ------------------------------------
drawComplexHeatmap <- Logical() ? function(gene_list = Character(), mat = Matrix(), sample_info_df = Data.frame(), primary_factor_column = Integer(), manual_color_map = Character(), output_image_file = Character())
{
	declare("gene_list_indices", Integer(null_ok = FALSE), value = which(rownames(mat) %in% gene_list));
	declare("mat_scaled", Matrix(null_ok = TRUE));
	if(length(gene_list_indices) == 0)
	{
		warn(file_logger, "Gene list given for generating heatmap has no match with transformed matrix, hence will skip drawing Heatmap for this data.");
		return(FALSE);
	}else if(length(gene_list_indices) == 1)
	{
		warn(file_logger, "Gene list has only 1 entry, hence will skip drawing Heatmap for this data.");
		return(FALSE);
	}else
	{
		mat_scaled <- t(apply(mat[gene_list_indices,], 1, scale));
	}
		
	if(class(mat_scaled) != "NULL" && nrow(mat_scaled) >= 2)
	{
		colnames(mat_scaled) <- paste(as.character(sample_info_df[,primary_factor_column]), rownames(sample_info_df), sep = "-", collapse = NULL);
		declare("d_cols", Any(), value = dist(x = t(mat_scaled), method = "euclidean"));
		declare("hc_cols", Any(), value = hclust(d = d_cols, method = "complete"));
		declare("d_rows", Any(), value = dist(x = mat_scaled, method = "euclidean"));
		declare("hc_rows", Any(), value = hclust(d = d_rows, method = "complete"));

		declare("complexHeatmap_object", Any());
		
		#ha <- ComplexHeatmap::HeatmapAnnotation(df = data.frame(condition = sample_info_df[, primary_factor_column]), col = list(condition = get_named_color_vectors(unique(as.character(sample_info_df[, primary_factor_column])))), annotation_legend_param = list(labels_gp = gpar(fontsize = 5), ncol = 2));
		declare("ha", Any(), value = ComplexHeatmap::HeatmapAnnotation(df = data.frame(condition = sample_info_df[, primary_factor_column], stringsAsFactors = FALSE), col = list(condition = manual_color_map), annotation_legend_param = list(labels_gp = gpar(fontsize = 5), ncol = 2)));
		
		declare("row_fontsize", Double(), value = (10*100/nrow(mat_scaled)));
		row_fontsize <- ifelse(test = (row_fontsize > 10), yes = 10, no = row_fontsize);
		declare("col_fontsize", Double(), value = 1);
		
		png(file = output_image_file, units = "in", width = 12, height = 10, res = 600, type = "cairo");
		complexHeatmap_object <- ComplexHeatmap::Heatmap(matrix =  mat_scaled, name = "Z-score", color = colorpanel(n = 1000, low = "blue" , mid = "white" , high = "red"), cluster_rows = hc_rows, row_dend_side = c("left"), row_dend_reorder = FALSE, row_dend_width = unit(20, "mm"), cluster_columns = hc_cols, column_dend_side = c("top"), column_dend_reorder = FALSE, column_dend_height = unit(20, "mm"), show_row_names = TRUE, show_column_names = TRUE, column_names_max_height = unit(20, "mm"), row_names_gp = gpar(fontsize = row_fontsize), column_names_gp = gpar(fontsize = col_fontsize), top_annotation = ha, heatmap_legend_param = list(labels_gp = gpar(fontsize = 5)));
		draw(complexHeatmap_object, heatmap_legend_side = "left", annotation_legend_side = "right");
		dev.off();
	}else
	{
		warn(file_logger, "No record selected to draw heatmap, hence will skip drawing Heatmap for this data.");
		return(FALSE);
	}
	return(TRUE);
}

## Draw interactive heatmap --------------------------------
drawIHeatmap <- Logical() ? function(gene_list = Character(), mat = Matrix(), row_ann_df = Data.frame(), row_colors_list = typed::List(), row_reorder = Logical(), col_ann_df = Data.frame(), col_colors_list = typed::List(), col_reorder = Logical(), output_file = Character())
{
	declare("gene_list_indices", Integer(null_ok = FALSE), value = which(rownames(mat) %in% gene_list));
	declare("mat_scaled", Matrix(null_ok = TRUE));

	if(length(gene_list_indices) == 0)
	{
		warn(file_logger, "Gene list given for generating heatmap has no match with transformed matrix, hence will skip drawing heatmap for this data.");
		return(FALSE);
	}else if(length(gene_list_indices) == 1)
	{
		warn(file_logger, "Gene list has only 1 entry, hence will skip drawing Heatmap for this data.");
		return(FALSE);
	}else
	{
		mat_scaled <- t(apply(mat[gene_list_indices,], 1, scale));
		colnames(mat_scaled) <- colnames(mat);
	}

	if(class(mat_scaled) != "NULL" && nrow(mat_scaled) >= 2)
	{
		declare("d_cols", Any());
		declare("hc_cols", Any());
		declare("d_rows", Any());
		declare("hc_rows", Any());
        declare("iheatmap_obj", Any());
		iheatmap_obj <- iheatmapr::main_heatmap(data = mat_scaled, name = "Z-score<br>(calculated on rows)");
		
		if(col_reorder)
        {
        	d_cols <- dist(x = t(mat_scaled), method = "euclidean");
            hc_cols <- hclust(d = d_cols, method = "complete");
                
            iheatmap_obj <- iheatmap_obj %>% add_col_dendro(hc_cols, reorder = col_reorder);
        }
        if(row_reorder)
		{
			d_rows <- dist(x = mat_scaled, method = "euclidean");
			hc_rows <- hclust(d = d_rows, method = "complete");

			iheatmap_obj <- iheatmap_obj %>% add_row_dendro(hc_rows, reorder = row_reorder);
		}
			
		iheatmap_obj <- iheatmap_obj %>%	
			add_row_annotation(row_ann_df, row_colors_list) %>%
			add_col_annotation(col_ann_df, col_colors_list) %>%
			add_row_labels() %>%
			add_col_labels();
                
		iheatmapr::save_iheatmap(p = iheatmap_obj, filename = output_file, selfcontained = TRUE);
		#iheatmapr::save_iheatmap(p = iheatmap_obj, filename = "test.png");
	}else
	{
		warn(file_logger, "No record selected to draw heatmap, hence will skip drawing Heatmap for this data.");
		return(FALSE);
	}
	return(TRUE);
}

## Identify gene clusters ----------------------------------
# This function is alternative to 'degPatterns'. But we are using cor() instead of cor.test() because the latter is very slow.
identifyGeneCluster <- Any() ? function(mat = Matrix(), sample_info_df = Data.frame(), minc = Integer(), sample_grouped_by = Character())
{
	declare("mat_scaled", Matrix(null_ok = TRUE));
	declare("mat_summarized", Matrix(null_ok = TRUE));
	declare("df_summarized_rearranged", Data.frame(null_ok = TRUE));
	declare("cor_mat", Matrix());
	declare("dist_mat", Matrix());
	declare("clusters", Any());
	declare("optimum_k", Integer());
	
	mat_scaled <- t(apply(mat, 1, scale));
	colnames(mat_scaled) <- colnames(mat);

	# Summarize the expression data based on factor levels of sample_grouped_by ---
	if(!is.factor(sample_info_df[,sample_grouped_by]))
	{
		sample_info_df[, sample_grouped_by] <- factor(sample_info_df[, sample_grouped_by], levels = unique(sample_info_df[, sample_grouped_by]));
	}
	mat_summarized <- matrix(data = NA, nrow = nrow(mat_scaled), ncol = length(levels(sample_info_df[, sample_grouped_by])));
	rownames(mat_summarized) <- rownames(mat_scaled);
	colnames(mat_summarized) <- levels(sample_info_df[, sample_grouped_by]);
	for(factor_level in levels(sample_info_df[, sample_grouped_by]))
	{
		for(i in 1:nrow(mat_scaled))
		{
			mat_summarized[i,factor_level] <- mean(mat_scaled[i, rownames(sample_info_df)[which(sample_info_df[, sample_grouped_by] == factor_level)]]);
		}
	}
	
	# Calculate correlation value for all pair-wise gene expression --------------------
	cor_mat <- cor(t(mat_summarized), use = "everything", method = "kendall");

	# Create a distance matrix based on correlation values -----------------------------
	dist_mat <- 1 - cor_mat;
    
	# Cluster gene-gene distance matrix using a divisive hierarchical clustering approach -------------------------------
	clusters <- cluster::diana(x = dist_mat, diss = TRUE, metric = "euclidean", stand = FALSE, stop.at.k = FALSE);
	#names(clusters); #"order"     "height"    "dc"        "merge"     "diss"      "call" "order.lab"

	# Identify optimum number of clusters and store membership -------------------------
	# cut the tree using the divisive coefficient of the clustering
	clusters$membership <- stats::cutree(tree = as.hclust(clusters), h = clusters$dc);

	# Create a ggplot object for all gene clusters which satisfy minc -------------------
	for(i in 1:ncol(mat_summarized))
	{
		df_summarized_rearranged <- rbind(df_summarized_rearranged, data.frame(gene_id = rownames(mat_summarized), condition = colnames(mat_summarized)[i], zscore = mat_summarized[,i]));
	}
	df_summarized_rearranged$condition <- factor(as.character(df_summarized_rearranged$condition), levels = unique(as.character(df_summarized_rearranged$condition)));
	
	# rename condition column to sample_grouped_by value
	colnames(df_summarized_rearranged)[which(colnames(df_summarized_rearranged) == "condition")] <- sample_grouped_by;
	
	declare("ggplot_list", typed::List(null_ok = TRUE));
	for(i in 1:length(unique(clusters$membership)))
	{
		declare("cluster_i_genes", Character());
		declare("temp_df", Data.frame());

		cluster_i_genes <- names(which(clusters$membership == i));
		if(length(cluster_i_genes) < minc)
		{
			next;
		}
		temp_df <- df_summarized_rearranged[which(df_summarized_rearranged$gene_id %in% cluster_i_genes),];
		
		ggplot_list[[i]] <- ggplot2::ggplot(data = temp_df, mapping = aes_string(x = base::as.name(sample_grouped_by), y = "zscore"));

		ggplot_list[[i]] <- ggplot_list[[i]] + theme_bw() +
			geom_boxplot() +
			geom_point(alpha = 0.3, color = "grey", position = position_jitter(width = 0.1)) +
			geom_line(mapping = aes(group = gene_id), alpha = 0.1, color = "grey", position = position_jitter(width = 0.1));
			
		ggplot_list[[i]] <- ggplot_list[[i]] +
			labs(title = paste("Group: ", i, " - ", "genes: ", length(cluster_i_genes), sep = ""), x = sample_grouped_by, y = "Scaled expression") +
			stat_summary(mapping = aes(group = 1), fun = median, geom = "line", color = "black") +
			stat_summary(fun = median, geom = "point");
		# Ref: https://stackoverflow.com/questions/3989987/joining-means-on-a-boxplot-with-a-line-ggplot2
	}

	# remove NULL item from the list
	ggplot_list <- ggplot_list[which(!sapply(X = 1:length(ggplot_list), FUN = function(i){is.null(ggplot_list[[i]])}))];
	
	# arrange multiple ggplots over multiple pages
	clusters$plot <- ggpubr::ggarrange(plotlist = ggplot_list, common.legend = TRUE, ncol = 2, nrow = 4);
	
	
	return(clusters);
}

identifyGeneClusterMfuzz <- Any() ? function(gene_list = Character(), mat = Matrix(), sample_info_df = Data.frame(), sample_order_by = Integer())
{
	#This function is incomplete
	declare("gene_list_indices", Integer(null_ok = FALSE), value = which(rownames(mat) %in% gene_list));
	declare("exprs", Matrix(null_ok = TRUE));
	declare("pheno_data", Any());
	declare("exprs_set", Any());
	declare("sample_order", Character());
	declare("optimum_fuzzifier_m", Double());
	declate("optimum_c", Integer());

	declare("clusters", Any());

	# Prepare expression data ----------
	if(length(gene_list_indices) == 0)
	{
		warn(file_logger, "Gene list given for identifying gene clusters has no match with transformed matrix, hence will skip this step.");
	}else if(length(gene_list_indices) == 1)
	{
		warn(file_logger, "Gene list has only 1 entry, hence will skip this step.");
	}else
	{
		exprs <- mat[gene_list_indices,];
		colnames(exprs) <- colnames(mat);
	}
	
	# Get required order of samples ------------
	sample_order <- unlist(lapply(X = levels(sample_info_df[, sample_order_by]), FUN =  function(x){ which(sample_info_df[, sample_order_by] == x)}));

	# Order sample ids of expression data as per sample_order_by ---
	exprs <- exprs[,rownames(sample_info_df[sample_order,])];

	# Prepare phenotypic data ------------------
	# keeping sample order same as that of exprs
	pheno_data <- new("AnnotatedDataFrame", data = sample_info_df[colnames(exprs),]);

	# Prepare an ExpressionSet object -----------
	exprs_set <- Biobase::ExpressionSet(assayData = exprs, phenoData = pheno_data);

	# Handle missing values ---------------------
	# Remove genes with more than 25% missing values
	exprs_set <- Mfuzz::filter.NA(eset = exprs_set, thres = 0.25);
	
	# Replace missing values with average expression of the corresponding gene
	exprs_set <- Mfuzz::fill.NA(eset = exprs_set, mode  = "knn", k = 10);
	
	# Standardize expressio values --------------
	# so that avg. expression is 0 and SD is 1 for each gene
	exprs_set <- Mfuzz::standardise(eset = exprs_set);
	
	# Identify optimum fuzzifier parameter, m ---
	optimum_fuzzifier_m <- Mfuzz::mestimate(exprs_set);

	# Identify optimum number of clusters, optimum_c
	# optimum_c would be the maximum one which leads to appearance of empty clusters which are defined when none of the genes have membership value larger than 0.5
	png(file = paste(output_dir, "/", output_prefix, "mfuzz_optimum_c.png", sep = "", collapse = ""), type = "cairo", units = "in", height = 10, width = 12, res = 600);
	optimum_c <- Mfuzz::cselection(eset = exprs_set, m = optimum_fuzzifier_m, crange = seq(2, 3^length(levels(sample_info_df[, primary_factor_column])), 2), repeats = 5, visu = TRUE);
	dev.off();

	# Initially cselection() did not seem to be effective as we didnt know what should be crange, I thought of writing a code that will search for optimum_c as shown below.
	if(FALSE)
	{
		# Unfortuunately, this code block for identification of optimum_c is taking too long.
		debug(file_logger, "Identifying optimum number of clusters ...");
		empty_cluster <- FALSE;
		optimum_c <- 0;
		while(!empty_cluster)
		{
			optimum_c <- optimum_c + 4;
			debug(file_logger, paste("Processing optimum_c: ", optimum_c, sep = ""));
			clusters <- Mfuzz::mfuzz(exprs_set, c = optimum_c, m = optimum_fuzzifier_m);
			for(cluster_i in 1:ncol(clusters$membership))
			{
				if(all(clusters$membership[,cluster_i] < 0.5))
				{
					empty_cluster <- TRUE;
				}
			}
		}
	}

	debug(file_logger, paste("Optimum number of clusters: ", optimum_c, sep = ""));
	

	png(file = paste(output_dir, "/", output_prefix, "mfuzz_optimum_Dmin.png", sep = "", collapse = ""), type = "cairo", units = "in", height = 10, width = 12, res = 600);
	optimum_Dmin <- Mfuzz::Dmin(eset = exprs_set, m = optimum_fuzzifier_m, crange = seq(4, optimum_c, 4), repeats = 5, visu = TRUE);
	dev.off();

	# Soft clustering of data -------------------
	clusters <- Mfuzz::mfuzz(exprs_set, c = optimum_c, m = optimum_fuzzifier_m);
	
	#names(clusters);
	#"centers"     "size"        "cluster"     "membership"  "iter"	"withinerror" "call"
	#clusters$cluster
	#clusters$membership
	
	# Cluster stability -------------------------
	png(file = paste(output_dir, "/", output_prefix, "mfuzz_plot.png", sep = "", collapse = ""), type = "cairo", units = "in", height = 10, width = 12, res = 600);
	Mfuzz::mfuzz.plot2(exprs_set, cl = clusters, xlab = colnames(sample_info_df)[primary_factor_column], ylab = "Z-score on expression data", time.labels = pData(exprs_set)[, primary_factor_column], mfrow = c(ceiling(length(clusters$size)/4), 4), min.mem = 0.7, x11 = FALSE, centre = TRUE, centre.col = "black", centre.lwd = 2);
	dev.off();
	# Note: yellow or green colored lines correspond to genes with low membership value; red and purple colored lines correspond to genes with high membership value. 

	# Get genelist associated to each cluster ----------

}


## =============================================================
## Read expression data, and sample info. table
## =============================================================
tryCatch(
	expr = 
	{
		## Check read access to input count data and/or sample info. table -----------
		if(file.access(names = count_file, mode = 4) != 0)
		{
			error(file_logger, paste("Unable to read --count_file file: ", count_file, sep = ""));
			base::sink();
			cat("Unable to read --count_file file: ", count_file, "\n", sep = "");
			cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Read expression data, and sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}else if(file.access(names = sample_info_file, mode = 4) != 0)
		{
			error(file_logger, paste("Unable to read --sample_info_file", sep = ""));
			base::sink();
			cat("Unable to read --sample_info_file", "\n", sep = "");
			cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Read expression data, and sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}

		declare("expression_df", Data.frame(null_ok = FALSE), value = as.data.frame(data.table::fread(input = count_file, sep = "auto", header = TRUE, check.names = FALSE), stringsAsFactors = FALSE));
		invisible(gc());
		
		## Check if gene_id_column is there in the counts table ----------------
		if(!gene_id_column %in% 1:ncol(expression_df))
		{
			fatal(file_logger, paste("User provided --gene_id_column number: ", gene_id_column, " is absent in the input counts file: ", count_file, sep = ""));
			base::sink();
			cat("User provided --gene_id_column number: ", gene_id_column, " is absent in the input counts file: ", count_file, "\n", sep = "");
			cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Read expression data, and sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}
		rownames(expression_df) <- as.character(expression_df[,gene_id_column]);
		rownames(expression_df) <- gsub(rownames(expression_df), pattern = "-", replacement = "_");

		## Check if count_column is there in the counts table ----------------
		if(!count_column %in% 2:ncol(expression_df))
		{
			fatal(file_logger, paste("User provided --count_column number: ", count_column, " is absent in the input counts file: ", count_file, ". It has to be >=2.", sep = ""));
			base::sink();
			cat("User provided --count_column number: ", count_column, " is absent in the input counts file: ", count_file, ". It has to be >=2.", "\n", sep = "");
			cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Read expression data, and sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}

		declare("gene_info_df", Data.frame(null_ok = TRUE));
		gene_info_df <- as.data.frame(expression_df[,1:(count_column - 1)]);
		rownames(gene_info_df) <- rownames(expression_df);
		colnames(gene_info_df) <- colnames(expression_df)[1:(count_column - 1)];
		
		expression_mat <- as.matrix(expression_df[,count_column:ncol(expression_df)]);
		debug(file_logger, "expression_mat[1:2,1:2]");
		print(expression_mat[1:2,1:2]);
		debug(file_logger, "dim(expression_mat)");
		print(dim(expression_mat)); # 60483 genes x 265 samples

		info(file_logger, paste("Input data has features: ", nrow(expression_mat), ", samples: ", ncol(expression_mat), sep = ""));

		declare("sample_info_df", Data.frame(null_ok = FALSE), value = as.data.frame(data.table::fread(input = sample_info_file, sep = "auto", header = TRUE, check.names = FALSE), stringsAsFactors = FALSE));

		## Check if sample_id_column is there in the sample info. table ----------------
		if(!sample_id_column %in% 1:ncol(sample_info_df))
		{
			fatal(file_logger, paste("User provided --sample_id_column number: ", sample_id_column, " is absent in the sample info. table: ", sample_info_file, sep = ""));
			base::sink();
			cat("User provided --sample_id_column number: ", sample_id_column, " is absent in the sample info. table: ", sample_info_file, "\n", sep = "");
			cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Read expression data, and sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}
		
		## Check if primary_factor_column is there in the sample info. table ------------------
		if(!primary_factor_column %in% 1:ncol(sample_info_df))
		{
			fatal(file_logger, paste("User provided --primary_factor_column number: ", primary_factor_column, " is absent in the sample info. table: ", sample_info_file, sep = ""));
			base::sink();
			cat("User provided --primary_factor_column number: ", primary_factor_column, " is absent in the sample info. table: ", sample_info_file, "\n", sep = "");
			cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Read expression data, and sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}

		## Check if all conditions mentioned in the contrasts are valid (for DESeq2_test_type: Wald) ---------------
		if(DESeq2_test_type == "Wald" && !contrasts == "all")
		{
			info(file_logger, paste("Selected contrast(s): ", contrasts, sep = ""));
			contrasts <- unlist(strsplit(arguments$contrasts, split = ","));
			
			declare("conditions_in_input_contrasts", Character(null_ok = FALSE), value = unique(unlist(strsplit(contrasts, split = "_vs_"))), const = TRUE);
			if(!all(conditions_in_input_contrasts %in% unique(as.character(sample_info_df[,primary_factor_column]))))
			{
				declare("invalid_conditions", Character(null_ok = TRUE), const = FALSE);
				invalid_conditions <- conditions_in_input_contrasts[!conditions_in_input_contrasts %in% unique(as.character(sample_info_df[,primary_factor_column]))];
				fatal(file_logger, paste("Following condition(s): ", paste(invalid_conditions, collapse = ","), " not present in primary_factor_column of sample table.", sep = ""));
				base::sink();
				cat("Following condition(s): ", paste(invalid_conditions, collapse = ","), " not present in primary_factor_column of sample table.", "\n", sep = "");
				cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Read expression data, and sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
				q(save = "no", status = 1);
			}
		}else if(DESeq2_test_type == "Wald" && contrasts == "all")
		{
			# Prepare all pairwise contrasts -----------------------------------------
			declare("conditions_in_input_contrasts", value = unique(as.character(sample_info_df[,primary_factor_column])), const = TRUE);
			declare("contrasts_matrix", Matrix(null_ok = FALSE), value = utils::combn(x = conditions_in_input_contrasts, m = 2));
			contrasts <- paste(contrasts_matrix[1,], contrasts_matrix[2,], sep = "_vs_");	
		}else if(DESeq2_test_type == "LRT")
		{
			# Contrasts are not required for LRT
			info(file_logger, "--contrasts will be ignored for --DESeq2_test_type = LRT as analysis will be performed based on --full_design and --reduced_design values if provided, else will be auto generated based on --primary_factor_column and --other_known_factors_columns.");
			contrasts <- NULL;
		}

		## Check if all other_known_factors_columns exist
		if(!is.null(other_known_factors_columns) && !all(other_known_factors_columns %in% 1:ncol(sample_info_df)) )
		{
			fatal(file_logger, "Not all other_known_factors_columns exist in sample_info_file.");
			base::sink();
			cat("Not all --other_known_factors_columns exist in sample_info_file.", "\n", sep = "");
			cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Read expression data, and sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}

		## If --full_design is defined, define other_known_factors_columns
		if(!is.null(full_design))
		{
			declare("temp_factors", Character(null_ok = TRUE));
			temp_factors <- unlist(strsplit(full_design, split = "[ ~+:]"));
			temp_factors <- temp_factors[!(temp_factors %in% c(""))];
			temp_factors <- unique(temp_factors);
			
			if(!all(temp_factors %in% colnames(sample_info_df)))
			{
				error(file_logger, paste("Not all factors: ", paste(temp_factors, collapse = ","), " mentioned in the --full_design are present in --sample_info_file columns", sep = ""));
				base::sink();
				cat("Not all factors mentioned in the --full_design are present in --sample_info_file columns", "\n", sep = "");
				cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Read expression data, and sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
				q(save = "no", status = 1);
			}

			temp_factors <- temp_factors[which(temp_factors %in% colnames(sample_info_df))];
			other_known_factors_columns <- which(colnames(sample_info_df) %in% temp_factors);
			other_known_factors_columns <- setdiff(other_known_factors_columns, primary_factor_column);
			if(length(other_known_factors_columns))
			{
				debug(file_logger, paste("Since --full_design was defined by the user, values for other_known_factors_columns were inferred which is/are: ", other_known_factors_columns, "\n", "These columns will be used in interactive clustered heatmaps. If --DESeq2_test_type = LRT but --reduced_design is not set, then these columns will be used to define --reduced_design."));
			}else
			{
				other_known_factors_columns <- NULL;
			}
			rm(temp_factors);
		}
		
		## If --reduced_design is defined, check its factors
		if(!is.null(reduced_design))
		{
			declare("temp_factors", Character(null_ok = TRUE));
			temp_factors <- unlist(strsplit(reduced_design, split = "[ ~+:]"));
			temp_factors <- temp_factors[!(temp_factors %in% c("", "1"))];
			temp_factors <- unique(temp_factors);
			if(!all(temp_factors %in% colnames(sample_info_df)))
			{
				error(file_logger, "Not all factors mentioned in the --reduced_design are present in --sample_info_file columns");
				base::sink();
				cat("Not all factors mentioned in the --reduced_design are present in --sample_info_file columns", "\n", sep = "");
				cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Read expression data, and sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
				q(save = "no", status = 1);
			}
		}
		
		## Remove duplicate entries (based on sample_id_column), and assign row names ---
		sample_info_df <- sample_info_df[!duplicated(sample_info_df[,sample_id_column]),];
		rownames(sample_info_df) <- as.character(sample_info_df[,sample_id_column]);

		info(file_logger, paste("Sample info table has ", nrow(sample_info_df), " samples", sep = ""));
	},
	error = function(e)
	{
		error(file_logger, "An error occured while reading expression data and sample info. table.");
		print(e);
		base::sink();
		print(e);
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Read expression data, and sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	},
	warning = function(w)
	{
		warn(file_logger, "Caught a warning while reading expression data and sample info. table. Something unexpected must have happened.");
		print(w);
		base::sink();
		print(w);
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Read expression data, and sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	}
)
save(list=ls(), file = paste(output_dir, "/", "DGE.RData", sep = ""));

## =========================================================
## Filter data based on sample info. table
## =========================================================
tryCatch(
	expr = 
	{
		## Filter data ----------------------------------------------------------------------------------
		info(file_logger, "Filtering data based on sample info. table ...");

		declare("common_sample_ids", Character());
		common_sample_ids <- intersect(rownames(sample_info_df), colnames(expression_mat));
		info(file_logger, paste("Number of sample ids common between input expression data and sample info. table: ", length(common_sample_ids), sep = ""));
		if(length(common_sample_ids) == 0)
		{
			error(file_logger, paste("Nothing is common between files provided for --count_file and --sample_info_file based on --sample_id_column. Hence, no further processing will occur.", sep = ""));
			base::sink();
			cat("Nothing is common between files provided for --count_file and --sample_info_file based on --sample_id_column. Hence, no further processing will occur.", "\n", sep = "");
			cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Filter data based on sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}

		expression_mat <- expression_mat[, common_sample_ids];
		debug(file_logger, "dim(expression_mat)");
		print(dim(expression_mat)); # 60483 Genes x 245 samples
		info(file_logger, paste("Filtered data has features: ", dim(expression_mat)[1], ", samples: ", dim(expression_mat)[2], sep = ""));

		sample_info_df <- sample_info_df[common_sample_ids,];
		debug(file_logger, "dim(sample_info_df)");
		print(dim(sample_info_df));
		info(file_logger, paste("Sample info. table now also contains the common sample ids only. It has #rows: ", dim(sample_info_df)[1], " #cols: ", dim(sample_info_df)[2], sep = ""));


		## For DGE analysis, it is important that the columns of sample_info_file which will be part of design formula of DESeqDataSet object should be available as factors
		sample_info_df[, primary_factor_column] <- factor(as.character(sample_info_df[, primary_factor_column]), levels = unique(as.character(sample_info_df[, primary_factor_column])));
				
		# repeat for other_known_factors_columns
		if(!is.null(other_known_factors_columns))
		{
			for(i in 1:length(other_known_factors_columns))
			{
				sample_info_df[, other_known_factors_columns[i]] <- factor(as.character(sample_info_df[, other_known_factors_columns[i]]), levels = unique(as.character(sample_info_df[, other_known_factors_columns[i]])));
			}
		}

		debug(file_logger, paste("List of conditions for the samples: ", paste(sample_info_df[, primary_factor_column], collapse = ","), sep = ""));
		debug(file_logger, "Original levels:");
		print(levels(sample_info_df[, primary_factor_column]));

		## Store all sample conditions (or sample groups) in a new variable -------------------------------
		conditions <- unique(as.character(sample_info_df[, primary_factor_column]));
		debug(file_logger, "print(conditions)");
		print(conditions);

		## Label conditions with distinct colors
		condition_col <- getNamedColorVectors(conditions);

		## Prepare data to be used for sample annotation in iheatmap ----------------------------------------
		declare("col_ann_df", Data.frame());
		declare("col_colors_list", typed::List());

		col_ann_df <- as.data.frame(sample_info_df[, c(primary_factor_column, other_known_factors_columns)]);
		rownames(col_ann_df) <- rownames(sample_info_df);
		colnames(col_ann_df) <- colnames(sample_info_df)[c(primary_factor_column, other_known_factors_columns)];

		col_colors_list <- list(list1 = condition_col)
		names(col_colors_list)[1] <- colnames(sample_info_df)[primary_factor_column]

	},
	error = function(e)
	{
		error(file_logger, "An error occured while filtering data based on sample info. table.");
		print(e);
		base::sink();
		print(e);
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Filter data based on sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	},
	warning = function(w)
	{
		warn(file_logger, "Caught a warning while filtering data based on sample info. table");
		print(w);
		base::sink();
		print(w);
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Filter data based on sample info. table'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	}
)
save(list=ls(), file = paste(output_dir, "/", "DGE.RData", sep = ""));

## =========================================================
## Data filter based on abundance
## =========================================================
tryCatch(
	expr =
	{
		## Check if input count data has all counts ---------
		if(!is.integer(expression_mat))
		{
			warn(file_logger, paste("Input count data does not look like Counts (integer values). Typically that happens when the data is generated through Salmon or RSEM quantification tools. Anyway, these will be converted to the nearest integer for further processing.", sep = ""));
			debug(file_logger, "Original input count data:");
			print(head(expression_mat));
			expression_mat <- round(expression_mat);
			debug(file_logger, "After conversion:");
			print(head(expression_mat));
		}

		## Save counts data (before filter) ------------------------------------------
		write.table(as.data.frame(cbind(gene_id = rownames(expression_mat), expression_mat)), file = paste(output_dir, "/", output_prefix, "counts_before_filter.csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE);

		## Generate box plot showing counts distribution (before filter) --------
		info(file_logger, "Generating box plot for showing counts distribution (before filter) ...");
		declare("temp_df", Data.frame());
		for(i in 1:ncol(expression_mat))
		{
			temp_df <- rbind(temp_df, data.frame(gene_id = rownames(expression_mat), sample_id = colnames(expression_mat)[i], condition = as.character(sample_info_df[colnames(expression_mat)[i], primary_factor_column]), expr = log2(expression_mat[,i])));
		}
		temp_df[which(temp_df == -Inf, arr.ind = TRUE)] <- 0;	#log2(0) = -Inf
		temp_df[which(temp_df == 'NaN', arr.ind = TRUE)] <- 0;	#log2(-ve) = NaN

		drawBoxPlot(data_df = temp_df, x = "sample_id", y = "expr", fill_by = "condition", manual_color_map = condition_col, title = "Counts distribution before filter", x_axis_title = "Sample", y_axis_title = "Counts (log2)", outfile = paste(output_dir, "/", output_prefix, "Counts_distribution_before_filter_box_plot", ".html", sep = ""));
		rm(temp_df);

		## Generate bar plot showing total counts (before filter) ---------------------------------------------
		info(file_logger, "Gnerating bar plot on total counts (before filter) ...");
		write.table(as.data.frame(cbind(sample_id = colnames(expression_mat), total_counts = colSums(expression_mat))), file = paste(output_dir, "/", output_prefix, "total_counts_before_filter.csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE);

		declare("temp_df", Data.frame());
		temp_df <- data.frame(sample_id = colnames(expression_mat), total_counts = colSums(expression_mat)/1000000, condition = as.character(sample_info_df[colnames(expression_mat),primary_factor_column]));
		drawBarPlot(data_df = temp_df, x = "sample_id", y = "total_counts", fill_by = "condition", manual_color_map = condition_col, title = "Total Counts before filter", x_axis_title = "Sample", y_axis_title = "Counts (in million)", outfile = paste(output_dir, "/", output_prefix, "Total_counts_before_filter_bar_plot", ".html", sep = ""));
		rm(temp_df);

		## Filter data by abundance -------------------------
		declare("expression_cpm_mat", Matrix(), value = edgeR::cpm(y = expression_mat, normalized.lib.sizes = FALSE, log = FALSE, prior.count = pseudocounts));
		declare("expression_lcpm_mat", Matrix(), value = edgeR::cpm(y = expression_mat, normalized.lib.sizes = FALSE, log = TRUE, prior.count = pseudocounts));
		expression_lcpm_mat[which(expression_lcpm_mat == -Inf, arr.ind = TRUE)] <- 0;	#log2(0) = -Inf
		expression_lcpm_mat[which(expression_lcpm_mat == 'NaN', arr.ind = TRUE)] <- 0;	#log2(-ve) = NaN
		invisible(gc());
		write.table(as.data.frame(cbind(gene_id = rownames(expression_lcpm_mat), expression_lcpm_mat)), file = paste(output_dir, "/", output_prefix, "expression_log2cpm_before_filter.csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE);

		debug(file_logger, "dim(expression_cpm_mat)");
		print(dim(expression_cpm_mat));
		debug(file_logger, "dim(expression_lcpm_mat)");
		print(dim(expression_lcpm_mat));
		debug(file_logger, "table(rowSums(expression_mat == 0) == ncol(expression_mat))");
		print(table(rowSums(expression_mat == 0) == ncol(expression_mat)));
		
		## Generate density plot (before filter) ----------------------
		info(file_logger, "Generating density plot (before filter) ...");
		declare("temp_df", Data.frame());
		for(i in 1:ncol(expression_lcpm_mat))
		{
			temp_df <- rbind(temp_df, data.frame(gene_id = rownames(expression_lcpm_mat), sample_id = colnames(expression_lcpm_mat)[i], condition = as.character(sample_info_df[colnames(expression_lcpm_mat)[i], primary_factor_column]), expr = expression_lcpm_mat[,i]));
		}
		drawDensityPlot(data_df = temp_df, x = "expr", group_by = "sample_id", fill_by = "condition", manual_color_map = condition_col, title = "Density plot before FilterByAbundance", x_axis_title = "log2(CPM)", y_axis_title = "Density", outfile = paste(output_dir, "/", output_prefix, "densityPlot_before_filter", ".html", sep = ""));
		rm(temp_df);


		info(file_logger, "Filtering out lowly expressed features ...");

		declare("keep_expr", Integer(null_ok = TRUE));

		if(use_CPM_cutoff)
		{
			declare("min_lib_size", Double(length = 1), value = min(colSums(expression_mat)));
			info(file_logger, paste("Library size for the least sequenced sample - ", paste(colnames(expression_mat)[which(colSums(expression_mat) == min_lib_size)], collapse = ","), ": ", min_lib_size, sep = ""));
			declare("cpm_cutoff", Double(length = 1), value = (counts_threshold/min_lib_size) * 1000000);
			info(file_logger, paste("CPM cutoff threshold, as calculated based on the least sequenced sample: ", "(", counts_threshold, "/", min_lib_size, ") * 1000000 = ", cpm_cutoff, sep = ""));
			info(file_logger, paste("Genes with CPM >= ", cpm_cutoff, " in at least ", sample_coverage_threshold * 100, "% samples, will be kept for downstream analysis.", sep = ""));
			keep_expr <- which(rowSums(expression_cpm_mat > cpm_cutoff) >= as.integer(ncol(expression_cpm_mat) * sample_coverage_threshold));
		}else
		{
			info(file_logger, paste("Genes with total counts >= ", counts_threshold, " will be kept for downstream analysis.", sep = ""));
			keep_expr <- which(rowSums(expression_mat) >= counts_threshold);
		}

		info(file_logger, paste("No. of features, after filter: ", length(keep_expr), sep = ""));

		if(length(keep_expr) == 0)
		{
			warn(file_logger, "There is no feature to further analyze the data. Hence, program will terminate here. Check the sample with lowest library size if that is creating issue. Sometimes very small value of library size (of the smallest sample) can create a large CPM cutoff which filters out all other features.");
			base::sink();
			cat("There is no feature to further analyze the data. Hence, program will terminate here. Check the sample with lowest library size if that is creating issue. Sometimes very small value of library size (of the smallest sample) can create a large CPM cutoff which filters out all other features.", "\n", sep = "");
			q(save = "no", status = 1);
		}

		expression_mat <- expression_mat[keep_expr,];
		debug(file_logger, "dim(expression_mat)");
		print(dim(expression_mat)); # 14198   245
		info(file_logger, paste("Filtered data has features: ", nrow(expression_mat), ", samples: ", ncol(expression_mat), sep = ""));
		
		expression_lcpm_mat <- edgeR::cpm(y = expression_mat, normalized.lib.sizes = FALSE, log = TRUE, prior.count = pseudocounts);
		expression_lcpm_mat[which(expression_lcpm_mat == -Inf, arr.ind = TRUE)] <- 0;	# log2(0) = -Inf
		expression_lcpm_mat[which(expression_lcpm_mat == 'NaN', arr.ind = TRUE)] <- 0;	# log2(-ve) = NaN

		## Save counts and lo2cpm data after filter -------------------
		write.table(as.data.frame(cbind(gene_id = rownames(expression_mat), expression_mat)), file = paste(output_dir, "/", output_prefix, "counts_after_filter.csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE);
		write.table(as.data.frame(cbind(gene_id = rownames(expression_lcpm_mat), expression_lcpm_mat)), file = paste(output_dir, "/", output_prefix, "expression_log2cpm_after_filter.csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE);

		## Generate box plot showing counts distribution (after filter) --------
		info(file_logger, "Generating box plot for showing counts distribution (after filter) ...");
		declare("temp_df", Data.frame());
		for(i in 1:ncol(expression_mat))
		{
			temp_df <- rbind(temp_df, data.frame(gene_id = rownames(expression_mat), sample_id = colnames(expression_mat)[i], condition = as.character(sample_info_df[colnames(expression_mat)[i], primary_factor_column]), expr = log2(expression_mat[,i])));
		}
		temp_df[which(temp_df == -Inf, arr.ind = TRUE)] <- 0;	#log2(0) = -Inf
		temp_df[which(temp_df == 'NaN', arr.ind = TRUE)] <- 0;	#log2(-ve) = NaN

		drawBoxPlot(data_df = temp_df, x = "sample_id", y = "expr", fill_by = "condition", manual_color_map = condition_col, title = "Counts distribution after filter", x_axis_title = "Sample", y_axis_title = "Counts (log2)", outfile = paste(output_dir, "/", output_prefix, "Counts_distribution_after_filter_box_plot", ".html", sep = ""));
		rm(temp_df);
		
		## Generate bar plot showing total counts (after filter) ---------------
		info(file_logger, "Generating bar plot on total counts (after filter) ...");
		
		write.table(as.data.frame(cbind(sample_id = colnames(expression_mat), total_counts = colSums(expression_mat))), file = paste(output_dir, "/", output_prefix, "total_counts_after_filter.csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE);

		declare("temp_df", Data.frame());
		temp_df <- data.frame(sample_id = colnames(expression_mat), total_counts = colSums(expression_mat)/1000000, condition = sample_info_df[colnames(expression_mat),primary_factor_column]);
		
		drawBarPlot(data_df = temp_df, x = "sample_id", y = "total_counts", fill_by = "condition", manual_color_map = condition_col, title = "Total Counts after filter", x_axis_title = "Sample", y_axis_title = "Counts (in million)", outfile = paste(output_dir, "/", output_prefix, "Total_counts_after_filter_bar_plot", ".html", sep = ""));
		
		rm(temp_df);

		## Draw density plot on filtered data --------------
		info(file_logger, "Generating density plot on filtered data ...");
		
		declare("temp_df", Data.frame());
		for(i in 1:ncol(expression_lcpm_mat))
		{
			temp_df <- rbind(temp_df, data.frame(gene_id = rownames(expression_lcpm_mat), sample_id = colnames(expression_lcpm_mat)[i], condition = as.character(sample_info_df[colnames(expression_lcpm_mat)[i], primary_factor_column]), expr = expression_lcpm_mat[,i]));
		}
		
		drawDensityPlot(data_df = temp_df, x = "expr", group_by = "sample_id", fill_by = "condition", manual_color_map = condition_col, title = "Density plot after FilterByAbundance", x_axis_title = "log2(CPM)", y_axis_title = "Density", outfile = paste(output_dir, "/", output_prefix, "densityPlot_after_filter", ".html", sep = ""));
		
		rm(temp_df);

		## Update top_n_for_viz (if required) --------------
		if(is.null(top_n_for_viz))
		{
			top_n_for_viz <- nrow(expression_mat)
		}else
		{
			top_n_for_viz[which(top_n_for_viz > nrow(expression_mat))] <- nrow(expression_mat);
		}
	},
	error = function(e)
	{
		error(file_logger, "An error occured while filtering data based on abundance.");
		print(e);
		base::sink();
		print(e);
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Data filter based on abundance'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	},
	warning = function(w)
	{
		warn(file_logger, "Caught a warning while filtering data based on abundance.")
		print(w);
		base::sink();
		print(w);
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Data filter based on abundance'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	}
)
save(list=ls(), file = paste(output_dir, "/", "DGE.RData", sep = ""));

## =========================================================
## Normalization and Transformation
## =========================================================

tryCatch(
	expr = 
	{
		## Construct DESeqDataSet object on the filtered data with no design -------
		info(file_logger, "Constructing DESeqDataSet object ...");
		dds <- DESeq2::DESeqDataSetFromMatrix(countData = expression_mat, colData = sample_info_df, design = as.formula(paste("~ ", 1, sep = "")));	# ~ 1 for no design, but remember to shift to other designs for DE
		debug(file_logger, "dim(dds)");
		print(dim(dds));

		## Normalization: Estimate size factors for each sample -------------------
		info(file_logger, "Estimating size factors for each sample ..");
		dds <- DESeq2::estimateSizeFactors(object = dds, type = "ratio", locfunc = stats::median);
		debug(file_logger, "print(sizeFactors(dds))");
		print(sizeFactors(dds));
		# design formula has no effect on normalization (size factors)

		## Draw density plot for normalized data ----------------------------------
		info(file_logger, "Generating density plot on normalized data ...");
		declare("expression_lcpm_mat_normalized", Matrix(), value = log2(DESeq2::fpm(object = dds, robust = TRUE)));
		expression_lcpm_mat_normalized[which(expression_lcpm_mat_normalized == -Inf, arr.ind = TRUE)] <- 0;		#log2(0) = -Inf
		expression_lcpm_mat_normalized[which(expression_lcpm_mat_normalized == 'NaN', arr.ind = TRUE)] <- 0;	#log2(-ve) = NaN
		invisible(gc());
		
		declare("temp_df", Data.frame());
		for(i in 1:ncol(expression_lcpm_mat_normalized))
		{
			temp_df <- rbind(temp_df, data.frame(gene_id = rownames(expression_lcpm_mat_normalized), sample_id = colnames(expression_lcpm_mat_normalized)[i], condition = as.character(sample_info_df[colnames(expression_lcpm_mat_normalized)[i], primary_factor_column]), expr = expression_lcpm_mat_normalized[,i]));
		}

		drawDensityPlot(data_df = temp_df, x = "expr", group_by = "sample_id", fill_by = "condition", manual_color_map = condition_col, title = "Density plot post normalization", x_axis_title = "log2(CPM)", y_axis_title = "Density", outfile = paste(output_dir, "/", output_prefix, "densityPlot_post_normalization", ".html", sep = ""));

		rm(temp_df);
		
		write.table(cbind(gene_id = rownames(expression_lcpm_mat_normalized), expression_lcpm_mat_normalized), file = paste(output_dir, "/", output_prefix, "expression_log2cpm_normalized.csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE);

		## Data transformation ----------------------------------------------------
		if( ncol(dds) < 30 )
		{
			info(file_logger, "Performing rlog transformation since #Samples < 30");
			normalized_transformed_mat <- SummarizedExperiment::assay(DESeq2::rlog(object = dds, blind = TRUE, fitType = "parametric"));
		}else
		{
			info(file_logger, "Performing VST transformation since #Samples >= 30");
			normalized_transformed_mat <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(object = dds, blind = TRUE, fitType = "parametric"));
		}
		info(file_logger, paste("The transormed data has: ", "#features: ", nrow(normalized_transformed_mat), ", #samples: ", ncol(normalized_transformed_mat), sep = ""));
		
		write.table(cbind(gene_id = rownames(normalized_transformed_mat), normalized_transformed_mat), file = paste(output_dir, "/", output_prefix, "normalized_transformed.csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE);
	},
	error = function(e)
	{
		error(file_logger, "An error occured while normalization and transformation.");
		print(e);
		base::sink();
		print(e);
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Normalization and Transformation'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	},
	warning = function(w)
	{
		warn(file_logger, "Caught a warning while normalization and transformation. Must have encountered something unexpected.");
		print(w);
		base::sink();
		print(w);
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Normalization and Transformation'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	}
)
save(list=ls(), file = paste(output_dir, "/", "DGE.RData", sep = ""));

## =========================================================
## Sample-level QC on transformed data
## =========================================================
tryCatch(
	expr = 
	{
		## Draw Mean-Variance relationship -------------------------------------------------------------
		info(file_logger, "Generating mean-variance relationship plot on transformed data ...");
		png(file = paste(output_dir, "/", output_prefix, "meanSdPlot_transformed_data.png", sep = "", collapse = ""), type = "cairo", units = "in", height = 4, width = 6, res = 300);
		#vsn::meanSdPlot(x = assay(transformed.data), ranks = TRUE, bins = nrow(transformed.data)); # x-axis (means) plotted on the rank scale allowing better visual assessment of the sd as a function of the mean
		vsn::meanSdPlot(x = normalized_transformed_mat, ranks = TRUE, bins = 50);
		dev.off();

		## Generate various Sample-level QC plots on all top_n_for_viz gene lists ----------------------
		for(i in 1:length(top_n_for_viz))
		{
			## Get list of N most variable genes across samples ------
			## This list is for visualization purpose only.
			declare("genes_for_viz", Character(null_ok = TRUE));
			if(top_n_for_viz[i] == nrow(normalized_transformed_mat))
			{
				genes_for_viz <- rownames(normalized_transformed_mat);
			}else
			{
				genes_for_viz <- filterByVariance(data_mat = normalized_transformed_mat, top_n = top_n_for_viz[i]);
			}

			## Draw t-SNE plot on transformed data ----------------------------------------------------------

			## Draw ComplexHeatmap on transformed data -----------------------------------------------------
			info(file_logger, paste("Generating clustered heatmap on transformed data using top ", top_n_for_viz[i], " most variable genes across samples ...", sep = ""));
			drawComplexHeatmap(gene_list = genes_for_viz, mat = normalized_transformed_mat, sample_info_df = sample_info_df, primary_factor_column = primary_factor_column, manual_color_map = condition_col, output_image_file = paste(output_dir, "/", output_prefix, "hclust_complexheatmap_on_transformed_data_using_top", top_n_for_viz[i], "_genes.png", sep = "", collapse = ""));

			## Draw interactive heatmap on transformed data ------------------------------------------------
			if(top_n_for_viz[i] > 1000)
			{
				info(file_logger, "Generating interactive clustered heatmap on genes > 1000 is not supported as it might slow down the processing. Hence, interactive clustered heatmap will be generated only on 1000 most variables genes across samples.");
				drawIHeatmap(gene_list = genes_for_viz[1:1000], mat = normalized_transformed_mat, row_ann_df = NULL, row_colors_list = NULL, row_reorder = TRUE, col_ann_df = col_ann_df, col_colors_list = col_colors_list, col_reorder = TRUE, output_file = paste(output_dir, "/", output_prefix, "hclust_iheatmap_on_transformed_data_using_top", "1000", "_genes.html", sep = "", collapse = ""));
			}else
			{
				info(file_logger, paste("Generating interactive clustered heatmap on transformed data using top ", top_n_for_viz[i], " most variable genes across samples ...", sep = ""));
				drawIHeatmap(gene_list = genes_for_viz, mat = normalized_transformed_mat, row_ann_df = NULL, row_colors_list = NULL, row_reorder = TRUE, col_ann_df = col_ann_df, col_colors_list = col_colors_list, col_reorder = TRUE, output_file = paste(output_dir, "/", output_prefix, "hclust_iheatmap_on_transformed_data_using_top", top_n_for_viz[i], "_genes.html", sep = "", collapse = ""));
			}

			## Draw MDSplot on transformed data ------------------------------------------------------------
			info(file_logger, paste("Generating MDS plot on transformed data using top ", top_n_for_viz[i], " most variable genes selected based on pairwise distance between samples ...", sep = ""));
			Glimma::glimmaMDS(normalized_transformed_mat, labels = rownames(sample_info_df), groups = sample_info_df, top = top_n_for_viz[i], gene.selection = "pairwise", launch = FALSE, html = paste(output_dir, "/", output_prefix, "MDSPlot_on_transformed_data_using_top", top_n_for_viz[i], "_genes.html", sep = ""), width = 1000, height = 500);

			## Draw PCA plot for all metadata columns upto PC4 using different sets of topN variable genes --
		}
	},
	error = function(e)
	{
		error(file_logger, "An error occured while QC on transformed data.");
		print(e);
		base::sink();
		print(e);
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Sample-level QC on transformed data'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	},
	warning = function(w)
	{
		warn(file_logger, "Caught a warning while QC on transformed data. Must have encountered something unexpected.");
		print(w);
		base::sink();
		print(w);
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Sample-level QC on transformed data'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	}
)
save(list=ls(), file = paste(output_dir, "/", "DGE.RData", sep = ""));

## =========================================================
## Identification of Differentially Expressed Genes (DEGs)
##                      with
## Useful graphical representation of DEGs
## =========================================================
tryCatch(
	expr = 
	{
		## Differential Expression analysis using DESeq2 -------
		info(file_logger, "Performing DGE analysis ...")

		## Define full and reduced design ----------------------
		## full_design/reduced design will be based on primary_factor_column and other_known_factors_columns, if not defined -------
		if(is.null(full_design))
		{
			if(is.null(other_known_factors_columns))
			{
				full_design <- as.formula(paste("~ ", colnames(sample_info_df)[primary_factor_column], sep = ""));
			}else
			{
				full_design <- as.formula(paste("~", paste(colnames(sample_info_df)[c(other_known_factors_columns, primary_factor_column)], collapse = " + "), sep = ""));
			}
		}else
		{
			full_design <- as.formula(full_design);
		}

		if(is.null(reduced_design))
		{
			if(is.null(other_known_factors_columns))
			{
				reduced_design <- as.formula(paste("~ ", 1, sep = ""));	# Secondary factors should be included here. For example, if full_design is ~ condition + batch, then reduced_design will be ~ batch, otherwise ~1.
			}else
			{
				reduced_design <- as.formula(paste("~", paste(colnames(sample_info_df)[other_known_factors_columns], collapse = " + "), sep = ""));
			}
		}else
		{
			reduced_design <- as.formula(reduced_design);
		}


		if(DESeq2_test_type == "Wald")
		{
			info(file_logger, "design formula for expression data (counts): "); print(full_design);

			# This step performs size factor estimation (if not done already), dispersion estimation, negative Binomial GLM fitting and Wald statistics
			design(dds) <- full_design;
			dds <- DESeq2::DESeq(object = dds,	test = "Wald", fitType = DESeq2_dispersion_method, sfType = "ratio", betaPrior = FALSE,	full = full_design, minReplicatesForReplace = 7, modelMatrixType = "standard", useT = FALSE, minmu = 0.5, parallel = FALSE);
		}
		else if(DESeq2_test_type == "LRT")
		{
			info(file_logger, "design formula for full model: "); print(full_design);
			info(file_logger, "design formula for reduced model: "); print(reduced_design);
			
			# This step performs size factor estimation (if not done already), dispersion estimation, negative Binomial GLM fitting and LRT statistics
			design(dds) <- full_design;
			declare("minmu", Double(), value = 0.5);
			if(DESeq2_dispersion_method == "glmGamPoi"){minmu <- 1e-06};
			dds <- DESeq2::DESeq(object = dds, test = "LRT", fitType = DESeq2_dispersion_method, sfType = "ratio", betaPrior = FALSE, full = full_design, reduced = reduced_design, minReplicatesForReplace = 7, modelMatrixType = "standard", useT = FALSE, minmu = minmu, parallel = FALSE);
		}
		#dispersions(dds) found to be same for Wald and LRT with full_design ~condition and reduced_design ~1. Not sure, if it will happen under all circumsatances.
		#betaPrior: whether or not to put a zero-mean normal prior on the non-intercept coefficients; default is FALSE for version >=1.16, and and shrunken LFCs are obtained afterwards using ‘lfcShrink’
		#reduced: for test=LRT
		#useT: passed to 'nbinomWaldTest'; whether to use t-distribution as null distribution; Default FALSE = Wald statistics are assumed to follow a standard normal
		
		
		# Draw gene-wise dispersion estimates -----------------------------
		png(file = paste(output_dir, "/", output_prefix, "genewise_disp_ests.png", sep = "", collapse = ""), units = "in", width = 12, height = 10, res = 600, type = "cairo");
		plotDispEsts(dds);
		dev.off();

		## Identify DEGs --------------------------------------------------
		if(DESeq2_test_type == "Wald")
		{
			## Remove contrasts which are not applicable ------------------
			# as few constrasts might not be valid anymore due to sample-based filtering
			declare("contrasts_to_be_removed", Integer(null_ok = TRUE));
			for(i in 1:length(contrasts))
			{
				if( ! all(unlist(strsplit(contrasts[i], split = "_vs_")) %in% conditions) )
				{
					contrasts_to_be_removed <- c(contrasts_to_be_removed, i);
				}
			}
			if( ! is.null(contrasts_to_be_removed) )
			{
				contrasts <- contrasts[-contrasts_to_be_removed];
			}
			info(file_logger, "Following contrasts will be considered for DGE: ");
			print(contrasts);

			for(i in 1:length(contrasts))
			{
				declare("test", Character(length = 1), value = unlist(strsplit(contrasts[i], split = "_vs_"))[1]);
				declare("ref", Character(length = 1), value = unlist(strsplit(contrasts[i], split = "_vs_"))[2]);
												
				## Note: For LFC shrinkage using 'apeglm', following steps are must:
				#	1. design needs to be updated with 'ref' as reference level,
				#	2. run nbinomWaldTest to re-estimate MLE coefficients
				#	3. run lfcShrink specifying the coefficient of interest in resultsNames(dds)

				## Update design with 'ref' as reference level -------------
				dds[[colnames(sample_info_df)[primary_factor_column]]] <- relevel(factor(dds[[colnames(sample_info_df)[primary_factor_column]]]), ref = ref);

				## Run nbinomWaldTest to re-estimate MLE coefficients ------
				dds <- DESeq2::nbinomWaldTest(dds, betaPrior = FALSE, modelMatrixType = "standard", betaTol = 1e-08, maxit = 100, useOptim = TRUE, useT = FALSE, minmu = 0.5);

				# Output files will be primarily stored in separate foders grouped by contrasts
				declare("output_dir_contrast", Character(length = 1, null_ok = FALSE), value = paste(output_dir, "/", contrasts[i], sep = ""));
				if(!dir.exists(output_dir_contrast))
				{
					dir.create(path = output_dir_contrast, recursive = TRUE, mode = "0777");
				}

				for(lfc_cutoff in log2_fold_change_cutoff)
				{
					## Store statistically significant gene entries only ---------------------------------------------------
					for(sc in stat_cutoff)
					{
						# Output files will be stored in separate foders grouped by stat_sig_measure, stat_cutoff and log2_fold_change_cutoff
						declare("output_dir_stat", Character(length = 1, null_ok = FALSE), value = paste(output_dir_contrast, "/", "log2FC", lfc_cutoff, "_", stat_sig_measure, sc, sep = ""));
						if(!dir.exists(output_dir_stat))
						{
							dir.create(path = output_dir_stat, recursive = TRUE, mode = "0777");
						}

						info(file_logger, paste("Processing ", contrasts[i], " at log2-fold-change threshold: ", lfc_cutoff, " and statistical significance threshold: ", sc, " ...", sep = ""));
						declare("res", Any());
						declare("res_lfc", Any());
						declare("res_combined", Any());
						
						res <- DESeq2::results(object = dds, contrast = c(colnames(sample_info_df)[primary_factor_column], test, ref), lfcThreshold = lfc_cutoff, alpha = sc, parallel = FALSE);
						
						debug(file_logger, "print(res): ");
						print(res);
						debug(file_logger, "print(summary(res)): ");
						print(summary(res));

						## Run lfcShrink with specific coefficient of interest ----------------------------------------------
						# This will replace 'log2FoldChange' and 'lfcSE' in res with shrunken LFC and SE.
						declare("cf", Integer());
						cf <- which(resultsNames(dds) == paste(colnames(sample_info_df)[primary_factor_column], "_", test, "_vs_", ref, sep = ""));
						res_lfc <- DESeq2::lfcShrink(dds, coef = cf, res = res, type = "apeglm", lfcThreshold = lfc_cutoff, svalue = TRUE, returnList = FALSE, format = "DataFrame", saveCols = NULL, apeAdapt = TRUE, parallel = FALSE);

						debug(file_logger, "print(res_lfc): ");
						print(res_lfc);
						debug(file_logger, "print(summary(res_lfc)): ");
						print(summary(res_lfc));

						res_combined <- data.frame(baseMean = res$baseMean, log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj, log2FoldChange_apeglm = res_lfc$log2FoldChange, lfcSE_apeglm = res_lfc$lfcSE, svalue = res_lfc$svalue, stringsAsFactors = FALSE);

						## add gene annotation ------------------
						res_combined <- cbind(gene_info_df[rownames(res_combined),], res_combined);
						colnames(res_combined)[1:ncol(gene_info_df)] <- colnames(gene_info_df);

						debug(file_logger, "dim(res_combined): ");
						print(dim(res_combined));
						
						declare("res_sig_genes", Data.frame(), value = res_combined, const = FALSE);
						res_sig_genes <- res_sig_genes[!is.na(res_sig_genes[,stat_sig_measure]),];

						res_sig_genes <- res_sig_genes[(abs(res_sig_genes[,effect_size_measure]) > lfc_cutoff) & (res_sig_genes[,stat_sig_measure] < sc),];
												
						info(file_logger, paste("#significant genes: ", nrow(res_sig_genes), " at log2-fold-change threshold: ", lfc_cutoff, " and ", stat_sig_measure, " < ", sc, sep = ""));
						
						write.table(res_sig_genes, file = paste(output_dir_stat, "/", output_prefix, contrasts[i], "_at_log2FC", lfc_cutoff, "_", stat_sig_measure, sc, "_DEGs.csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE);

						if(nrow(res_sig_genes) > 0)
						{
							DGE_result <- rbind(DGE_result, cbind(contrast = contrasts[i], log2_fold_change_cutoff = lfc_cutoff, stat_cutoff = sc, res_sig_genes));
						}
						debug(file_logger, "dim(DGE_result): ");
						print(dim(DGE_result));

						## Generate MA plot ------------------------------------------------------------------------------------
						info(file_logger, paste("Generating MA plot with genes highlighted with |log2-fold-change| > ", lfc_cutoff, " at ", stat_sig_measure, " < ", sc, " ...", sep = ""));
						declare("html_output_file", Character(length = 1));
						html_output_file <- paste(output_dir_stat, "/", output_prefix, contrasts[i], "_MAPlot_at_log2FC", lfc_cutoff, "_", stat_sig_measure, sc, ".html", sep = "");
						
						## Not to use Glimma::glimmaMA as it automatically takes only one contrast based on groups, and uses log2_fold_change_cutoff of 0 and stat_cutoff of unknown? to determine status. Instead, use Glimma::glimmaXY.
						#Glimma::glimmaMA(dds, groups = sample_info_df[, primary_factor_column], transform.counts = "logcpm", main = "MA plot", xlab = "log2CPM", ylab = "log2FC", launch = FALSE, html = html_output_file, width = 920, height = 920);

						declare("regulation_status", Double(null_ok = TRUE), value = rep(0, nrow(res_combined)), const = FALSE);
						regulation_status[which((res_combined[,effect_size_measure] > lfc_cutoff) & (res_combined[,stat_sig_measure] < sc))] <- 1;
						regulation_status[which((res_combined[,effect_size_measure] < lfc_cutoff) & (res_combined[,stat_sig_measure] < sc))] <- -1;
						regulation_status[which(is.na(res_combined[,stat_sig_measure]))] <- NA;
						Glimma::glimmaXY(x = rowMeans(expression_lcpm_mat), y = res_combined[,effect_size_measure], xlab = "avgExpression_log2CPM", ylab = effect_size_measure, counts = expression_mat, groups = sample_info_df[, primary_factor_column], status = regulation_status, anno = res_combined, display.columns = setdiff(colnames(res_combined), colnames(expression_df)[gene_id_column]), status.cols = c("dodgerblue", "silver", "firebrick"), sample.cols = NULL, transform.counts = "logcpm", main = paste("MA plot: with genes highlighted with |log2-fold-change| > ", lfc_cutoff, " at ", stat_sig_measure, " < ", sc, sep = ""), html = html_output_file, width = 1600, height = 920);

						## Generate Volcano plot --------------------------------------------------------------------------------
						info(file_logger, "Generating Volcano plot ...");
						html_output_file <- paste(output_dir_stat, "/", output_prefix, contrasts[i], "_VolcanoPlot_at_log2FC", lfc_cutoff, "_", stat_sig_measure, sc, ".html", sep = "");
						Glimma::glimmaXY(x = res_combined[,effect_size_measure], y = -log10(res_combined[,stat_sig_measure]), xlab = effect_size_measure, ylab = paste("neg_log10_", stat_sig_measure, sep = ""), counts = expression_mat, groups = sample_info_df[, primary_factor_column], status = regulation_status, anno = res_combined, display.columns = setdiff(colnames(res_combined), colnames(expression_df)[gene_id_column]), status.cols = c("dodgerblue", "silver", "firebrick"), sample.cols = NULL, transform.counts = "logcpm", main = paste("Volcano plot: with genes highlighted with |log2-fold-change| > ", lfc_cutoff, " at ", stat_sig_measure, " < ", sc, sep = ""), html = html_output_file, width = 1400, height = 920);

						## Generate interactive heatmap --------------------------------------------------------------------------
						if(nrow(res_sig_genes) > 1)
						{
							info(file_logger, paste("Generating interactive clustered heatmap on transformed data using DEGs ...", sep = ""));
							html_output_file <- paste(output_dir_stat, "/", output_prefix, contrasts[i], "_hclust_iheatmap_on_transformed_data_using_DEGs_at_log2FC", lfc_cutoff, "_", stat_sig_measure, sc, ".html", sep = "", collapse = "");
							drawIHeatmap(gene_list = rownames(res_sig_genes), mat = normalized_transformed_mat, row_ann_df = NULL, row_colors_list = NULL, row_reorder = TRUE, col_ann_df = col_ann_df, col_colors_list = col_colors_list, col_reorder = TRUE, output_file = html_output_file);
						}else
						{
							info(file_logger, "Interactive clustered heatmap on DEGs skipped as number of genes is not more than 1");
						}
					}
				}
			}
		}else if(DESeq2_test_type == "LRT")
		{
			## Store statistically significant gene entries only ---------------------------------------------------
			for(sc in stat_cutoff)
			{
				# Output files will be stored in separate foders grouped by stat_sig_measure and stat_cutoff
				declare("output_dir_stat", Character(length = 1, null_ok = FALSE), value = paste(output_dir, "/", stat_sig_measure, sc, sep = ""));
				if(!dir.exists(output_dir_stat))
				{
					dir.create(path = output_dir_stat, recursive = TRUE, mode = "0777");
				}

				info(file_logger, paste("Processing ", "at statistical significance threshold: ", sc, " ...", sep = ""));
				declare("res", Any());
				res <- DESeq2::results(object = dds, alpha = sc, parallel = FALSE);
							
				debug(file_logger, "print(res): ");
				print(res);

				debug(file_logger, "summary(res): ");
				summary(res);

				## Remove log2FoldChange and lfcSE columns as these are not directly associated with the actual hypothesis test [Reference: 24]
				res <- res[, -which(colnames(res) %in% c("log2FoldChange", "lfcSE"))];
				#class(res); DESeqDataSet
				res <- data.frame(res);
				#class(res); data.frame
				
				## add gene annotation ------------------
				res <- cbind(gene_info_df[rownames(res),], res);
				colnames(res)[1:ncol(gene_info_df)] <- colnames(gene_info_df);

				declare("res_sig_genes", Data.frame(), value = res, const = FALSE);
				res_sig_genes <- res_sig_genes[!is.na(res_sig_genes[, stat_sig_measure]),];
				res_sig_genes <- res_sig_genes[res_sig_genes[, stat_sig_measure] < sc,];
							
				info(file_logger, paste("#significant genes: ", nrow(res_sig_genes), " at ", stat_sig_measure, " < ", sc, sep = ""));
							
				write.table(res_sig_genes, file = paste(output_dir_stat, "/", output_prefix, "acrossLevels_at_", stat_sig_measure, sc, "_DEGs.csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE);

				DGE_result <- rbind(DGE_result, cbind(stat_cutoff = sc, res_sig_genes));
								
				## Generate ComplexHeatmap -------------------------------------------------------------------------------
				

				## Generate interactive heatmap --------------------------------------------------------------------------
				if(nrow(res_sig_genes) > 1)
				{
					info(file_logger, paste("Generating interactive clustered heatmap on transformed data using DEGs ...", sep = ""));
					html_output_file <- paste(output_dir_stat, "/", output_prefix, "hclust_iheatmap_on_transformed_data_using_DEGs_at_", stat_sig_measure, sc, ".html", sep = "", collapse = "");
					drawIHeatmap(gene_list = rownames(res_sig_genes), mat = normalized_transformed_mat, row_ann_df = NULL, row_colors_list = NULL, row_reorder = TRUE, col_ann_df = col_ann_df, col_colors_list = col_colors_list, col_reorder = TRUE, output_file = html_output_file);
				}else
				{
					info(file_logger, "Interactive clustered heatmap on DEGs skipped as number of genes is not more than 1");
				}
				
				## Identify gene clusters exhibiting particular patterns across samples ---------------------------------
				declare("clusters", Any());
				if(nrow(res_sig_genes) > 15)
				{
					info(file_logger, "Identifying gene clusters ...");

					# Subset transformed normalized counts to keep only DEGs
					declare("normalized_transformed_mat_subset", Matrix());
					normalized_transformed_mat_subset <- normalized_transformed_mat[rownames(res_sig_genes),];
					
					# Identify gene clusters using hierarchical clustering approach based on pair-wise correlation
					if(FALSE)
					{
						#This section did not work as expected. Hence, will skip this. Instead will be using custom function identifyGeneCluster().
						if(is.null(other_known_factors_columns))
						{
							clusters <- DEGreport::degPatterns(ma = normalized_transformed_mat_subset, metadata = sample_info_df, minc = 15, summarize = "merge", time = colnames(sample_info_df)[primary_factor_column], col = NULL, consensusCluster = FALSE, reduce = FALSE, cutoff = 0.7, scale = TRUE, pattern = NULL, groupDifference = NULL, eachStep = FALSE, plot = FALSE, fixy = NULL);
							clusters$plot <- clusters$plot + theme(legend.position = "none");
						}else
						{
							#only one other known factor is considered for plot
							clusters <- DEGreport::degPatterns(ma = normalized_transformed_mat_subset, metadata = sample_info_df, minc = 15, summarize = "merge", time = colnames(sample_info_df)[primary_factor_column], col = colnames(sample_info_df)[other_known_factors_columns[1]], consensusCluster = FALSE, reduce = FALSE, cutoff = 0.7, scale = TRUE, pattern = NULL, groupDifference = NULL, eachStep = FALSE, plot = TRUE, fixy = NULL);
						
						}
					}

					clusters <- identifyGeneCluster(mat = normalized_transformed_mat[rownames(res_sig_genes),], sample_info_df, minc = 15, sample_grouped_by = colnames(sample_info_df)[primary_factor_column]);		

					ggpubr::ggexport(clusters$plot, filename = paste(output_dir_stat, "/", output_prefix, "gene_clusters_on_DEGs_at_", stat_sig_measure, sc, ".png", sep = ""), width = 800, height = 800);					
					
					clusters$df <- cbind(gene_info_df[names(clusters$membership),], clusters$membership);
					clusters$df <- as.data.frame(clusters$df)
					colnames(clusters$df)[1:ncol(gene_info_df)] <- colnames(gene_info_df);
					colnames(clusters$df)[ncol(clusters$df)] <- "cluster";

					debug(file_logger, "print(head(clusters$df))");
					print(head(clusters$df));

					write.table(clusters$df, file = paste(output_dir_stat, "/", output_prefix, "gene_cluster_membership_on_DEGs_at_", stat_sig_measure, sc, ".csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE);
				}else
				{
					info(file_logger, "Identification of gene clusters could not be performed due to less number of input features (<15).");
				}
			}
		}

		## Store merged DGE result -------------------------
		write.table(DGE_result, file = paste(output_dir, "/", output_prefix, "merged_differential_expression_result.csv", sep = ""), sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE);
	},
	error = function(e)
	{
		error(file_logger, "An error occured while identifying DEGs.");
		print(e);
		base::sink();
		print(e);
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Identification of DEGs'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	},
	warning = function(w)
	{
		warn(file_logger, "Caught a warning while idenifying DEGs.");
		print(w);
		base::sink();
		print(w);
		cat("The program: Differential_gene_expression_analysis_on_counts.R failed at 'Identifcation of DEGs'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1)
	}
)
save(list=ls(), file = paste(output_dir, "/", "DGE.RData", sep = ""));


info(file_logger, "Processing ends");
info(file_logger, "R session info:");

## Store R session information -----------------------------
print(sessionInfo());

## Store citation of the major packages involved -----------
declare("current_Rsession", Any(), value = sessionInfo());
print(current_Rsession);
if( length(current_Rsession$otherPkgs) > 0 )
{
	declare("other_packages", value = unlist(lapply(X = current_Rsession$otherPkgs, FUN = function(x){x$Package})));
	for(pkg in other_packages)
	{
		print(citation(package = pkg), bibtex = FALSE);
	}
}

base::sink();
cat("The program: Differential_gene_expression_analysis_on_counts.R executed successfully!", "\n", sep = "");
