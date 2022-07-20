#!/usr/bin/Rscript --vanilla

#########################################
## Program: Functional_enrichment.R
## Purpose: Perform functional enrichment analysis on a list of genes (optionally along with some numeric weightage like log2fold-change which might be originated from a differential expression analysis pipeline such as DESeq2). This will help to identify pathways/terms which are enriched by DEGs.
##
## Usage: /path/to/Rscript --vanilla /path/to/Functional_enrichment.R --genelist_file <path/to/gene list> --gene_id_column <column number for gene identifiers, e.g., 1> --log2fc_column <column number for log2-fold-change, e.g., 5> --reference_genelist_file <path/to/reference gene list> --ref_gene_id_column <column number for gene identifiers in reference gene list file> --species <name of the species, e.g., "Homo sapiens"> --data_sources <comma separated list of data sources for analysis, such as MSigDBHallmarks,MSigDBOncogenic,MSigDBImmunologic,DO,NCG,GOBP,GOCC,GOMF,Reactome> --run_SPIA <whether to run SPIA on KEGG database> --threshold <used for inferring pathway significance> --output_dir <path/to/output directory> --output_prefix <prefix for output files such as figures> --log_file <filename for log>
##
## Developer: Srikant Verma (srikant_verma@persistent.com)
## Version: 1.0: Created on May 27, 2022
## Version: 1.1: Update on June 7, 2022
## Version: 1.2: Updated on June 10, 2022
## Version: 1.3: Updated on June 12, 2022
## Version: 1.4: Updated on June 20, 2022
## Version: 1.5: Updated on June 28, 2022
#########################################

## Version: 1.0: Created on May 27, 2022
##  It performs the following major tasks:
#   1. Performs over-representation analysis (ORA) on following data sources:
#		MSigDB Hallmark gene sets
#		MSigDB oncogenic signatures
#		MSigDB immunologic signatures
#		Disease Ontology (DO) terms [For only Homo sapiens]
#		Network of Cancer Genes (NCG) database [For only Homo sapiens]
#		DisGeNET database [For only Homo sapiens]
#		GO (Biological Process)
#		GO (Molcular Function)
#		GO (Cellular Component)
#		KEGG pathways [code written but need KEGG license, hence will not run by default.]
#		Reactome pathways
#		
#   2. Performs gene set enrichment analysis (GSEA) or functional class scoring (FCS) approach based analysis on following data sources:
#		MSigDB Hallmark gene sets
#		MSigDB oncogenic signatures
#		MSigDB immunologic signatures
#		Disease Ontology (DO) terms [For only Homo sapiens]
#		Network of Cancer Genes (NCG) database [For only Homo sapiens]
#		DisGeNET database [For only Homo sapiens]
#		GO (Biological Process)
#		GO (Molecular Function)
#		GO (Cellular Component)
#		KEGG pathways [code written but need KEGG license, hence will not run by default.]
#		Reactome pathways
#		
#   3. Performs pathway topology (PT) based enrichment analysis on following data sources:
#		KEGG [For only Homo sapiens and Mouse]
#
#	4. Topology based enrichment analysis is being performed by SPIA ( which currently uses old KEGG.db (2012), but a licensee can use latest one; modification in script is required)
#
#	5. Supports major gene identifier types: entrezgene_id, ensembl_gene_id and uniprot_gn_symbol as input.
#
#	6. Since data sources - DO, NCG, DisGeNET, KEGG and Reactome work only with entrezgene_id, if user provides gene lists with other identifier such as ensembl_gene_id or uniprot_gn_symbol, it will be converted into entrezgene_id using biomaRt. SPIA uses KEGG database, hence it also works only with entrezgene_id.
#
#	7. Only following species tested so far: Homo sapiens.
#
#	8. It is advisable to use full list of genes (used in experiment - for example those genes kept after abundance based filtering) for GSEA because of several reasons including the one that says modest differences in expression can also be biologically relevant and applying an arbitrary threshold might remove those genes from the analysis [6]. However, I feel keeping such genes whose expression difference estimates are not statistically significant, might also be questionable. Therefore, for this version, the list of DEGs will be used for GSEA.
#
#	9. To store the enrichment results, we will prefer TSV over CSV because there is/are module/s which generates outputs containing comma, such as ReactomPA.
#
#	10. Often, the list of enriched GO terms is too long and contains redundant terms, which hinders effective interpretation. We are using clusterProfiler::simplify() function which employs GOSemSim R package to calculate semantic similarity among enriched terms. Highly similar terms (e.g., > 0.7) are removed to retain a representattive term (e.g., most significant term). [10]
#
#	11. Visualiation of enrichment results using lollipop plot (png + html).
#


## Version: 1.1: Updated on June 7, 2022
#	Minor changes:
#		genelist$genes will hold character values via as.character().
# 		Unique Ensembl gene ids found gene annotation table will not have NA or "".
#		Updated some logging messages.
#	
#	Enhancement:
#		Often connection to Ensembl mart does not get established on the first attempt, which leads to program failure. In this version, we have introduced a piece of code that will retry upto 10 times before giving up. This will reduce the chances of program failure due to connection issues.

## Version: 1.2: Updated on June 10, 2022
#	Enhancement:
#		Further enhancement for better connection to Ensembl Mart

## Version: 1.3: Updated on June 12, 2022
#	Enhancement:
#		Removed warning section from tryCatch when it is not required to quit on encounter.
#	Bug fixes:
#		1. variable 'reattempt' was not getting updated in 'error' section of tryCatch using the assignment operator '<-'. Now, we are using superassignment operator (<<-) instead to do the job.
#

## Version: 1.4: Updated on June 20, 2022
#	Bug fix:
#		When run_SPIA = TRUE but data_sources does not contain 'KEGG', then following error encountered:  "<simpleError in normalizePath(dirname(file), mustWork = TRUE): path[1]="temp_out/KEGG": No such file or directory>". This is now fixed. Whenever, --run_SPIA will be set TRUE, 'KEGG' sub-directory will be created under --output_dir.
#			

## Version: 1.5: Updated on June 28, 2022
#	Enhancements/minor issue fix:
#		1. Validate --gene_id_column, --log2fc_column, --ref_gene_id_column
#		2. Handle spaces (around comma for more than one input values) in following parameters: --data_sources
#		3. Send more error messages to STDOUT (which were previously directed to log files only).
#		4. Default value of --gene_id_column is now set to 1
#		5. Default value of --log2fc_column is now set to NULL
#		6. Handle NULL values to parameters (such as -log2fc_column, --output_prefix) when it is provided as command line arguments to these parameters.
#		7. DisGeNET is commented out as of now due to License issue
#		8. Removd save(list=ls(), file = paste(output_dir, "/", "FunctionalEnrichment.RData", sep = "")) from warning and error blocks of tryCatch() because it only saves warning/error message.
#		9. Updated list of supported species based on msigdbr::msigdbr_species()$species_name and ?ReactomePA::enrichPathway().
#	Bug fixes:
#		1. Fixed Ensembl gene id pattern in identifyGeneIdType() as described here: https://asia.ensembl.org/info/genome/stable_ids/index.html. Stable Ensembl IDs have the following format: ENS[species prefix][feature type prefix][a unique eleven digit number].
#		2. Fixed Uniport gene symbol pattern in identifyGeneIdType() as described here: https://www.genenames.org/about/old-guidelines/
#		3. Fixed Entrez gene id pattern in identifyGeneIdType()
#		4. Fixed an error where enricher/enrichGO/enrichDO etc. was returning NULL which was leading to following error message: "Error in UseMethod("filter") : no applicable method for 'filter' applied to an object of class "NULL"".


## Known issues:
#	1. SPIA uses old (free) version of KEGG.
#	2. clusterProfiler::enrichKEGG uses old (free) version of KEGG with parameter 'use_internal_data'. If set TRUE, it will use KEGG.db generated in 2012, while FALSE will download the latest data.
#	3. enrichplot R package provides a number of functions to plot enrichment output. However, when we try to convert the ggplot object to plotly (for html), the legends looks wrong. Therefore, we are going to use either custom scripts for generating html outputs or will stick with enrichplot library to generate png output.
#	4. GSEA results using clusterProfiler/fgsea might vary from one run to the other on the same data. This has been reported by several users.
#		References:
#			https://github.com/YuLab-SMU/clusterProfiler/issues/190
#			https://github.com/ctlab/fgsea/issues/12




## Future releases:

## Notes:
#
## What are the objectives of the different approaches for enrichment analysis:
#	ORA: Objective of ORA is to determine whether genes from pre-defined sets (e.g., those belonging to a specific GO term or KEGG pathway) are present more than would be expected (over-represented) in a subset of user input list of genes [5]. It only consider number of genes and ignores any values asscoiated with them such as fold-change.
#	GSEA: Objective of GSEA is to assess whether members of a gene set (like an MSigDB hallmark gene set) appear enriched at one end of the input ranked list of genes[4]. It considers values associated with genes in analysis.
#	Topology: Objective of pathway topology based approach is to consider the pathway topology to compute gene-level statistics.
#
## How to perform GSEA?
#The analysis is performed by:
##	1. ranking all genes in the data set
#	2. identifying the rank positions of all members of the gene set in the ranked data set
#	3. calculating an enrichment score (ES) that represents the difference between the observed rankings and that which would be expected assuming a random rank distribution.


## How to interpret the results of GSEA? [4]
#	Enrichment score (ES) reflects the degree to which a gene set is overrepresented at the top or bottom of the input ranked list of genes.
#	Positive normalized enrichment score: genes in the gene set will be mostly represented at the top of input list
#	Negative normalized enrichment score: genes in the gene set will be mostly represented at the bottom of input list
#	About leading edge: "As described in the Gene Set Enrichment Analysis PNAS paper, the leading-edge subset in a gene set are those genes that appear in the ranked list at or before the point at which the running sum reaches its maximum deviation from zero. The leading-edge subset can be interpreted as the core that accounts for the gene set’s enrichment signal." [8]

## What is the difference between over-representation and GSEA? [2]
#	Over Resentation Analysis (ORA): considers user defined genes selected after applying some threshold (e.g., DEGs after applying p.adj threshold and/or log2FC)
#	Gene Set Enrichment Analysis (GSEA): considers all the genes in the experiment (e.g., genes filtered after abundance threshold, which is input to DESeq2) because it assumes that weaker but coordinated changes in sets of functionally related genes can also have significant effects.
#
# ** It is advised to use full list of genes (used in experiment - for example those genes kept after abundance based filtering) for GSEA because of several reasons including the one that says modest differences in expression can also be biologically releveant and applying an arbitrary threshold might remove those genes from the analysis [6]. However, I feel keeping such genes whose expression difference estimates are not statistically significant, might also be questionable. 
#	
#

## Interpretation of various plots: [7]
#	Lollipop plot for ORA output: Enrichment analysis using ORA typically generates following columns in the output table: 'ID', 'Description', 'p.adjust', 'Count', 'BgRatio', 'GeneRatio'. Count denotes the number of input genes present in the term. BgRatio is ratio of all genes that are annotated in this term. BgRatio 200/4383 means there are 200 genes out of total 4383 genes (across all available terms) in the given term. GeneRatio is ratio of input genes that are annotated in a term. GeneRatio is ratio of input genes that are annotated in a term. GeneRatio 53/241 means there are 53 genes from de_genes which are part of the term while total 241 genes from de_genes would map to total 4383 genes. In the plot, the x-axis represents 'Rich Factor' which is defined as the ratio of input genes (e.g., DEGs) that are annotated in a term to all genes that are annotated in this term. For above example, rich factor = 53/200. The y-axis represents the term ID or Description, and is ordered on the basis of 'Rich Factor'. The size of the dots denote the number of input genes annotated to a term. The color denotes 'p.adjust'. [10] 
#
#	Lollipop plot for GSEA output: The x-axis represents the Normalized Enrichment Score (NES) of the terms, and y-axis represents the ID/Description of the terms. Dots are colored based on term's 'padj' value. In case of pathway, the NES indicates the shift of genes belonging to a certain pathway toward either end of the ranked list and represents pathway activation or suppression. NES is an indicator to interpret the degree of enrichment. A positive NES indicates that members of the gene set tend to appear at the top of the rank (pathway activation), and a negative NES indicates the opposite circumstance (pathway suppression). [10]
#	 
#	Lollipop plot for SPIA output: SPIA generates a data frame containing the ranked pathways and various statistics: ‘pSize’ is the number of genes on the pathway; ‘NDE’ is the number of DE genes per pathway; ‘tA’ is the observed total preturbation accumulation in the pathway; ‘pNDE’ is the probability to observe at least ‘NDE’ genes on the pathway using a hypergeometric model; ‘pPERT’ is the probability to observe a total accumulation more extreme than ‘tA’ only by chance; ‘pG’ is the p-value obtained by combining ‘pNDE’ and ‘pPERT’; ‘pGFdr’ and ‘pGFWER’ are the False Discovery Rate and respectively Bonferroni adjusted global p-values; and the ‘Status’ gives the direction in which the pathway is perturbed (activated or inhibited). ‘KEGGLINK’ gives a web link to the KEGG website that displays the pathway image with the differentially expressed genes highlighted in red. [Ref: ?spia]. In the plot, the x-axis represents 'Rich Factor' which is defined as the ratio of input genes (e.g., DEGs) that are annotated in a KEGG pathway to all genes that are annotated in this pathway. For example, if NDE is 12 and pSize is 128, then rich factor = 12/128. The y-axis represents the KEGG pathway name, and is ordered on the basis of 'Rich Factor'. The size of the dots denote the number of input genes (NDE) annotated to a term. The color denotes 'p.adjust'.
#
#	Upset plot: This plot helps to know which genes are involved in the multiple significant terms. [7]

## References:
#   1. Ten years of pathway analysis: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002375
#   2. Difference between ORA and GSEA: http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/FAQ#What_is_the_difference_between_GSEA_and_an_overlap_statistic_.28hypergeometric.29_analysis_tool.3F
#   3. Difference between ORA and GSEA: http://crazyhottommy.blogspot.com/2016/08/gene-set-enrichment-analysis-gsea.html
#   3. How to view GSEA result: https://software.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Viewing_Analysis_Results
#   4. How to interpret GSEA enrichment score: https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?_GSEAPreranked_Page
#	5. Objective of ORA: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/#:~:text=Over%2Drepresentation%20(or%20enrichment),a%20subset%20of%20your%20data.
#	6. GSEA paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1239896/
#	7. How to use clusterProfiler: https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
#	8. GSEA guide: https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html
#	9. fgsea output interpretation: https://bioconductor.org/packages/devel/bioc/manuals/fgsea/man/fgsea.pdf
#	10. clusterProfiler publication: https://www.cell.com/the-innovation/pdf/S2666-6758(21)00066-7.pdf



##########################################


cat("For help: /path/to/Rscript --vanilla /path/to/Functional_enrichment.R --help", "\n", sep = "");


## =========================================================
## Take user specified parameters
## =========================================================
library("optparse");	# command line option parser

args_list <- list(
    optparse::make_option(opt_str = c("--genelist_file"), type = "character", default = NULL, help = "List of differentially expressed genes in CSV or TSV format. Default: NULL"),
	optparse::make_option(opt_str = c("--gene_id_column"),			type = "integer",			default = 1,							help = "Column number in --genelist_file which contains gene identifiers. Default: 1"),
    optparse::make_option(opt_str = c("--log2fc_column"),			type = "integer",			default = NULL,							help = "Column number to be used for log2-fold-change information in the --genelist_file. In DESeq2, if statistical significance measure is pvalue or padj, the moderated log2fc is usually available in column number 3. If, however, the statistical significance measure is svalue, shrunken log2fc is estimated by 'apeglm' method and the estimated value is typically available in column number 8. Default: NULL"),
	optparse::make_option(opt_str = c("--reference_genelist_file"), type = "character", default = NULL, help = "List of reference gene list in CSV or TSV format. This is used in over-representation analysis (ORA) and topology based analysis (SPIA). Gene set enrichment analysis (GSEA) does not require this list. It should include all entries from --genelist_file. Typically, it the list of all genes which are part of the study, for example, it could be list of genes after abundance based filter that goes into differential expression analysis. Default: NULL"),
	optparse::make_option(opt_str = c("--ref_gene_id_column"), type = "integer", default = 1, help = "Column number in --reference_genelist_file which contains gene identifiers. Default: 1"),
	optparse::make_option(opt_str = c("--species"), type = "character", default = "Homo sapiens", help = "Name of the species. This will be used to retrieve species-specific annotation data sources for enrichment analysis. Currently only tested for following: Homo sapiens. Default: Homo sapiens"),
	optparse::make_option(opt_str = c("--data_sources"),			type = "character",			default = "all",						help = "Comma separated list of data sources for analysis, such as MSigDBHallmarks,MSigDBOncogenic,MSigDBImmunologic,DO,NCG,,GOBP,GOCC,GOMF,Reactome. Default: all"),
	optparse::make_option(opt_str = c("--run_SPIA"), type = "logical", default = TRUE, help = "Whether to run SPIA. It requires --log2fc_column. Default: TRUE."),
	optparse::make_option(opt_str = c("--threshold"), type = "double", default = 0.1, help = "A numerical threshold (between 0 and 1) on p.adjust to infer geneset/pathway significance. Default: 0.1."),
	optparse::make_option(opt_str = c("--output_dir"), type = "character", default = "./", help = "Output directory path. Default: ./"),
	optparse::make_option(opt_str = c("--output_prefix"), type = "character", default = NULL, help = "Prefix for output files. Default: NULL"),
	optparse::make_option(opt_str = c("--log_file"), type = "character", default = "Functional_enrichment.log", help = "Log file. Default: Functional_enrichment.log")
)
arguments <- optparse::parse_args(optparse::OptionParser(option_list = args_list));


## =========================================================
## Load required libraries
## =========================================================
library("typed");		    # To set types for a variable, and for arguments of a function
library("log4r");		    # Logging library
library("data.table");	    # Read input data fast
library("ggplot2");			# For drawing plots
library("ggsci");			# To support Nature Publishing Group (NPG) color palettes
library("gplots");			# To generate smoothly varying set of colors
library("plotly");			# For interactive visualization
library("htmlwidgets");		# For saving html widgets
library("msigdbr");			# Provides Molecular Signatures Database (MSigDB) gene sets
library("biomaRt");			# Provides annotation for genes
library("SPIA");			# For KEGG signaling pathway impact analysis
library("GOSemSim");		# For calculating semantic similarities among enriched terms
library("clusterProfiler");	# For Over-representation analysis (ORA) on MSigDB gene sets and GO terms for DEGs
library("DOSE");			# For Over-representation analysis (ORA) and GSEA on Disease Ontology terms, Network of Cancer Genes (NCG) gene sets and DisGeNET gene sets for DEGs
library("fgsea");			# For GSEA on MSigDB gene sets
library("ReactomePA");		# For Over-representation analysis (ORA) and GSEA on Reactome pathways
#library("forcats");		# For reordering factor levels
#library("enrichplot");		# Provides visualizations for enrichment results
#library("ggupset");		# For generating UpSet plot
#library("stringr");		# Wraps long strings


#library("ggpubr");			# For exploratory data analysis (boxplot)

# To avoid following error, Error in .External2(C_X11, paste0("png::", filename), g$width, g$height,
options(bitmapType = 'cairo');


## =========================================================
## Set imp. variables
## =========================================================
declare("genelist_file", Character(length = 1, null_ok = FALSE), value = arguments$genelist_file, const = TRUE);
declare("gene_id_column", Integer(length = 1, null_ok = FALSE), value = arguments$gene_id_column, const = TRUE);

if(!is.null(arguments$log2fc_column) && (arguments$log2fc_column == 'NULL'))
{
	arguments$log2fc_column <- NULL;
}
declare("log2fc_column", Integer(length = 1, null_ok = TRUE), value = arguments$log2fc_column, const = TRUE);
declare("reference_genelist_file", Character(length = 1, null_ok = FALSE), value = arguments$reference_genelist_file, const = TRUE);
declare("ref_gene_id_column", Integer(length = 1, null_ok = TRUE), value = arguments$ref_gene_id_column, const = TRUE);
declare("species", Character(length = 1, null_ok = FALSE), value = arguments$species, const = TRUE);
declare("data_sources", Character(null_ok = FALSE), value = gsub(unlist(strsplit(arguments$data_sources, split = ",")), pattern = "\\s", replacement = ""));
declare("run_SPIA", Logical(length = 1), value = arguments$run_SPIA);
declare("threshold", Double(length = 1), value = arguments$threshold);
declare("output_dir", Character(length = 1, null_ok = FALSE), value = arguments$output_dir, const = TRUE);

if(!is.null(arguments$output_prefix) && (arguments$output_prefix == 'NULL'))
{
	arguments$output_prefix <- NULL;
}
declare("output_prefix", Character(length = 1, null_ok = TRUE), value = arguments$output_prefix, const = FALSE);
output_prefix <- ifelse(test = is.null(output_prefix), yes = "", no = paste(output_prefix, "_", sep = ""));

declare("log_file", Character(length = 1, null_ok = FALSE), value = arguments$log_file, const = FALSE);


#declare("supported_data_sources", Character(), value = c("MSigDBHallmarks", "MSigDBOncogenic", "MSigDBImmunologic", "DO", "NCG", "DisGeNET", "GOBP", "GOCC", "GOMF", "KEGG", "Reactome"));
declare("supported_data_sources", Character(), value = c("MSigDBHallmarks", "MSigDBOncogenic", "MSigDBImmunologic", "DO", "NCG", "", "GOBP", "GOCC", "GOMF", "Reactome"));
declare("supported_species_for_MSigDB", Character(), value = msigdbr::msigdbr_species()$species_name);
declare("supported_species_for_ReactomePA", Character(), value = c("Homo sapiens", "Rattus norvegicus", "Mus musculus", "Caenorhabditis elegans", "Saccharomyces cerevisiae", "Danio rerio", "Drosophila melanogaster"));
declare("supported_species", Character(), value = unique(c(supported_species_for_MSigDB, supported_species_for_ReactomePA)));

declare("genelist_df", Data.frame());								# Gene list data - frame
declare("genelist", typed::List());									# List of all differentially expressed genesets
declare("ref_genelist", Character(null_ok = TRUE));					# Reference gene list
declare("gene_info_df", Data.frame());								# Gene annotation data frame based in bioMart
declare("mapped_entrez_gene_ids", Character());						# Mapped named vector
declare("enriched_MSigDB_hallmark_genesets_using_enricher_df", Data.frame());
declare("enriched_MSigDB_oncogenicC6_genesets_using_enricher_df", Data.frame());
declare("enriched_MSigDB_ImmunologicC7_genesets_using_enricher_df", Data.frame());
declare("enriched_DO_genesets_using_enrichDO_df", Data.frame());
declare("enriched_NCG_genesets_using_enrichNCG_df", Data.frame());
declare("enriched_DGN_genesets_using_enrichDGN_df", Data.frame());
declare("enriched_GOBP_genesets_using_enrichGO_df", Data.frame());
declare("enriched_GOMF_genesets_using_enrichGO_df", Data.frame());
declare("enriched_GOCC_genesets_using_enrichGO_df", Data.frame());
declare("enriched_kegg_pathways_using_enrichKEGG_df", Data.frame());
declare("enriched_reactome_pathways_using_enrichPathway_df", Data.frame());
declare("enriched_MSigDB_hallmark_genesets_using_fgsea_df", Data.frame());
declare("enriched_MSigDB_oncogenicC6_genesets_using_fgsea_df", Data.frame());
declare("enriched_MSigDB_immunologicC7_genesets_using_fgsea_df", Data.frame());
declare("enriched_DO_genesets_using_gseDO_df", Data.frame());
declare("enriched_NCG_genesets_using_gseNCG_df", Data.frame());
declare("enriched_DGN_genesets_using_gseDGN_df", Data.frame());
declare("enriched_GOBP_genesets_using_gseGO_df", Data.frame());
declare("enriched_GOMF_genesets_using_gseGO_df", Data.frame());
declare("enriched_GOCC_genesets_using_gseGO_df", Data.frame());
declare("enriched_kegg_pathways_using_gseKEGG_df", Data.frame());
declare("enriched_reactome_pathways_using_gsePathway_df", Data.frame());
declare("enriched_kegg_pathways_using_spia_df", Data.frame());


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


## Check validity of some parameters -----------------------
if(!(species %in% supported_species))
{
	error(file_logger, paste("User provided --species: ", species, " is not supported currently. Only following are supported: ", paste(supported_species, collapse = "/"), sep = ""));
	base::sink();
	cat("User provided --species: ", species, " is not supported currently. Only following are supported: ", paste(supported_species, collapse = "/"), "\n", sep = "")
	cat("The program: Functional_enrichment.R failed. Check log for more details: ", log_file, "\n", sep = "");
	q(save = "no", status = 1);
}

if(is.null(log2fc_column) && run_SPIA)
{
	warn(file_logger, "SPIA will work only if genes are associated with some log2-fold-change values, which is not the case here.");
	run_SPIA <- FALSE;
}

if(length(data_sources) == 1 && data_sources == 'all')
{
	data_sources = supported_data_sources;
}else
{
	if(!all(data_sources %in% supported_data_sources))
	{
		error(file_logger, paste("Only following data sources are supported for enrichment analysis", paste(supported_data_sources, collapse = ","), sep = ""));
		base::sink();
		cat("Only following data sources are supported for enrichment analysis", paste(supported_data_sources, collapse = ","), "\n", sep = "");
		cat("The program: Functional_enrichment.R failed. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	}
}

## =========================================================
## Custom functions
## =========================================================

## Extract list of genes ----------------------------------------
extractGenes <- typed::List() ? function(genelist_df = Data.frame(), gene_id_column = Integer(), log2fc_column = Integer())
{
	declare("genelist", typed::List(null_ok = TRUE));
	
	# Entire gene list will be considered a single list
	genelist$genes <- as.character(genelist_df[, gene_id_column]);

	if(!is.null(log2fc_column))
	{
		genelist$log2fc <- genelist_df[, log2fc_column];
	}

	return(genelist);
}

## Identify type of gene identifier - Ensembl or Entrez or Uniprot Gene Symbol ---------
identifyGeneIdType <- Character() ? function(ids = Character())
{
	declare("id_type", Character());

	if(all(grepl(ids, pattern = "^ENS[A-Z]*G[0-9]{11}")))
	{
		id_type <- "ensembl_gene_id";
	}else if(all(grepl(ids, pattern = "^[0-9]+$")))
	{
		id_type <- "entrezgene_id";
	}else if(all(grepl(ids, pattern = "^[A-Za-z][A-Za-z0-9-]+")))
	{
		id_type <- "uniprot_gn_symbol";
	}else
	{
		id_type <- "Unknown";
	}
	return(id_type);
}

## Connect to bioMart for gene annotation -----------------------
annotateWithBioMart <- Data.frame() ? function(biomart_database_name = Character(), species = Character(), gene_id_type = Character(), gene_ids = Character(), file_logger = Any())
{
	declare("ensembl", Any());
	declare("ensembl_version_df", Data.frame());
	declare("gene_info_df", Data.frame());

	ensembl <- tryCatch(
    	expr = {
			info(file_logger, paste("Connecting to ", "https://www.ensembl.org", " ...", sep = ""));
			ensembl <- biomaRt::useMart(biomart = biomart_database_name, host = "https://www.ensembl.org");
		},
		error = function(e)
		{
			error(file_logger, "Following error encountered: ");
			print(e);
			tryCatch(
				expr = {
					info(file_logger, paste("Connecting to ", "https://useast.ensembl.org", " ...", sep = ""));
					ensembl <- biomaRt::useMart(biomart = biomart_database_name, host = "https://useast.ensembl.org");
				},
				error = function(e)
				{
					error(file_logger, "Following error encountered: ");
					print(e);
					tryCatch(
						expr = {
							info(file_logger, paste("Connecting to ", "https://uswest.ensembl.org", " ...", sep = ""));
							ensembl <- biomaRt::useMart(biomart = biomart_database_name, host = "https://uswest.ensembl.org");
						},
						error = function(e)
						{
							error(file_logger, "Following error encountered: ");
							print(e);
							tryCatch(
								expr = {
									info(file_logger, paste("Connecting to ", "https://asia.ensembl.org", " ...", sep = ""));
									ensembl <- biomaRt::useMart(biomart = biomart_database_name, host = "https://asia.ensembl.org");
								},
								error = function(e)
								{
									error(file_logger, "Failed while connecting to bioMart for gene annotation");
									print(e);
									base::sink();
									print(e);
									cat("The program: Annotation_using_OpenTargets.R failed because it could not connect to bioMart for gene annotation.", "\n", sep = "");
									q(save = "no", status = 1);
								}
							)
						}
					)
				}
			)
		}
	)

	declare("data_retrieved", Logical(), value = FALSE);
	declare("reattempt", Integer(), value = as.integer(0));
	while(!data_retrieved)
	{
		tryCatch(
			expr = 
			{
				if(reattempt > 0)
				{
					info(file_logger, paste("Re-attempt number: ", reattempt, sep = "")); 
				}
				
				ensembl_version_df <- as.data.frame(biomaRt::listDatasets(ensembl) %>% filter(grepl(gsub(tolower(species), pattern = "([a-z]).*?\\s", replacement = "\\1"),dataset)), stringsAsFactors = FALSE);
				info(file_logger, "Ensembl gene annotation version to be used: ");
				print(ensembl_version_df);
				ensembl <- biomaRt::useDataset(dataset = ensembl_version_df$dataset, mart = ensembl);

				gene_info_df <- as.data.frame(biomaRt::getBM(attributes = c('ensembl_gene_id', 'entrezgene_id', 'uniprot_gn_symbol', 'external_gene_name', 'external_gene_source', 'description', 'gene_biotype', 'source', 'version'), filters = gene_id_type, values = gene_ids, mart = ensembl), stringsAsFactors = FALSE);

				data_retrieved <- TRUE;
			},
			error = function(e)
			{
				error(file_logger, "Following error encountered while retrieving data from Ensembl.");
				print(e);

				if(reattempt < 10)
				{
					info(file_logger, "Will retry ...");
					reattempt <<- as.integer(reattempt + 1);
				}else
				{
					fatal(file_logger, "Giving up now! Program will terminate now!");
					base::sink();
					print(e);
					cat("The program: Annotation_using_OpenTargets.R failed while retrieving data from bioMart.", "\n", sep = "");
					q(save = "no", status = 1);
				}
			}
		)   
	}

	declare("perc_overlap", Double());
	perc_overlap <- sum(gene_ids %in% gene_info_df[,gene_id_type])/length(gene_ids)*100;
	if(perc_overlap < 10)
	{
		warn(file_logger, paste("Only ", perc_overlap, "% of input genes (DEGs + reference genes) could be mapped to Ensembl!", "\n", "Check whether --species was correctly set.", sep = ""));
	}
	else if(perc_overlap < 50)
	{
		warn(file_logger, paste("Only ", perc_overlap, "% of input genes (DEGs + reference genes) could be mapped to Ensembl!", sep = ""));
	}else
	{
		info(file_logger, paste(perc_overlap, "% of input genes (DEGs + reference genes) could be mapped to Ensembl!", sep = ""));
	}

	return(gene_info_df);
}

## Draw bar plot -------------------------------------------
drawBarPlot <- Logical() ? function(data_df = Data.frame(), x = Character(), y = Character(), fill_by = Character(), title = Character(), x_axis_title = Character(), y_axis_title = Character(), outfile = Character())
{
	declare("p", Any());

	# reorder the appearance of y component. Y component should be a categorical variable which will be converted to a factor for reordering purpose based on x variable.
	data_df[, y] <- factor(data_df[,y], levels = data_df[order(data_df[,x], decreasing = TRUE), y]);

	p <- ggplot2::ggplot(data = data_df, mapping = aes_string(x = base::as.name(x), y = base::as.name(y))) +
		ggtitle(title) + xlab(x_axis_title) + ylab(y_axis_title);
	p <- p + ggplot2::geom_bar(mapping = aes_string(fill = base::as.name(fill_by)), stat  = "identity", show.legend = TRUE);
	p <- p + ggplot2::scale_fill_viridis_c();
	p <- p + theme_bw();

	declare("fig", Any(), value = plotly::ggplotly(p));
	htmlwidgets::saveWidget(widget = fig, file = outfile, selfcontained = TRUE);

	return(TRUE);
}

## Draw lollipop plot --------------------------------------
drawLollipopPlot <- Logical() ? function(data_df = Data.frame(), x = Character(), y = Character(), color_by = Character(), size_by = Character(), title = Character(), x_axis_title = Character(), y_axis_title = Character(), outfile = Character())
{
	declare("p", Any());
	declare("fig", Any());
	
	# reorder the appearance of y component. Y component should be a categorical variable which will be converted to a factor for reordering purpose based on x variable.
	data_df[, y] <- factor(data_df[,y], levels = data_df[order(data_df[,x], decreasing = TRUE), y]);

	if(FALSE)
	{
		# This section is experimental, and an alternative to ggplot2 + ggplotly
		p <- plotly::plot_ly(data = data_df,
			x = as.formula(paste("~", base::as.name(x), sep = "")), y = as.character(data_df[,y]),
			text = as.formula(paste("~", base::as.name(y), sep = "")),
			type = "scatter", mode = "markers", color = as.formula(paste("~", base::as.name(color_by), sep = "")),
			marker = list(size = as.formula(paste("~", size_by, sep = "")))
		);
		
		p <- p %>% 
			plotly::layout(
				title = title,
				xaxis = list(title = x_axis_title, showgrid = TRUE),
				yaxis = list(title = y_axis_title, showgrid = TRUE, categoryorder = "array", categoryarray = as.character(data_df[,y])),
				showlegend = TRUE
			)

		htmlwidgets::saveWidget(widget = p, file = outfile, selfcontained = TRUE);
	}

	p <- ggplot2::ggplot(data = data_df, mapping = aes_string(x = base::as.name(x), y = base::as.name(y))) +
			ggtitle(title) + xlab(x_axis_title) + ylab(y_axis_title);
	p <- p + ggplot2::geom_segment(mapping = aes_string(xend = 0, yend = base::as.name(y)));
	p <- p + ggplot2::geom_point(data = data_df, mapping = aes_string(color = base::as.name(color_by), size = base::as.name(size_by)), show.legend = TRUE);
	p <- p + ggplot2::scale_colour_viridis_c();
	p <- p + theme_bw();

	fig <- plotly::ggplotly(p);
	htmlwidgets::saveWidget(widget = fig, file = outfile, selfcontained = TRUE);

	# As currently ggplotly does not reproduce size legend, the above HTML will not have the same. Therefore, we will be generating a PNG file to overcome this. [https://github.com/plotly/plotly.R/issues/705, https://stackoverflow.com/questions/68180064/ggplotly-does-not-reproduce-size-legend] 
	png(filename = sub(pattern = "\\.html", replacement = ".png", x = outfile), width = 12, height = 10, units = "in", res = 300, type = "cairo");
	plot(p);
	dev.off();
	
	return(TRUE);
}


## Draw UpSet plot ----------------------------------------------------
drawUpSetPlot <- Logical() ? function(data_df = Data.frame(), x = Character(), y = Character(), title = Character(), outfile = Character())
{
	declare("p", Any());
	declare("fig", Any());

	p <- ggplot2::ggplot(data = data_df, mapping = aes_string(x = base::as.name(x))) +
		ggtitle(title);
	p <- p + ggplot2::geom_bar();
	p <- p + ggupset::scale_x_upset(n_intersections = Inf);

	#HTML version of the plot is not coming nicely, hence will not generate
	#fig <- plotly::ggplotly(p);
	#htmlwidgets::saveWidget(widget = fig, file = outfile, selfcontained = TRUE);

	png(filename = sub(pattern = "\\.html", replacement = ".png", x = outfile), width = 12, height = 10, units = "in", res = 300, type = "cairo");
	plot(p);
	dev.off();

	return(TRUE);

}


## Run ORA on MSigDB gene sets ----------------------------------
runORAonMSigDB <- Data.frame() ? function(de_genes = Character(), ref_genes = Character(), id_type = Character(), msigdb_gs_df = Data.frame(), gene_ann_df = Data.frame(), threshold = Double(), output_file_prefix = Character())
{
	declare("enrich_result", Any());
	declare("term2gene_df", Data.frame());
	
	if(id_type == 'ensembl_gene_id')
	{
		term2gene_df <- msigdb_gs_df[,c("gs_name", "ensembl_gene")];
	}else if(id_type == 'entrezgene_id')
	{
		term2gene_df <- msigdb_gs_df[,c("gs_name", "entrez_gene")];
	}else if(id_type == 'uniprot_gn_symbol')
	{
		term2gene_df <- msigdb_gs_df[,c("gs_name", "gene_symbol")];
	}

	# Peform ORA ---------------------------------------------------
	enrich_result <- clusterProfiler::enricher(gene = de_genes, universe = ref_genes, pvalueCutoff = threshold, pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10, maxGSSize = Inf, TERM2GENE = term2gene_df);

	# Return an empty data.frame if no gene set have size > 10 which satisfies pvalueCutoff = threshold
	if(is.null(enrich_result))
	{
		enrich_result <- as.data.frame(enrich_result);
		return(enrich_result);
	}

	enrich_result <- enrich_result %>% filter(p.adjust < threshold);

	## Add rich factor of enriched terms -----------------------------
	# GeneRatio = k/m, BgRatio = M/N, then richFactor = k/M
	enrich_result <- clusterProfiler::mutate(enrich_result, richFactor = Count/as.numeric(sub(pattern = "/\\d+", replacement = "", x = BgRatio)));
	
	## Save enrichment output as a data.frame -----------------------
	enrich_result <- as.data.frame(enrich_result);
	#colnames(enrich_result);	#ID Description GeneRatio BgRatio pvalue p.adjust qvalue geneID Count richFactor

	enrich_result <- enrich_result[order(enrich_result$richFactor, decreasing = TRUE),];
	
	## Add gene symbols, if ids are not the same already ----------
	if(nrow(enrich_result))
	{
		if(id_type != 'uniprot_gn_symbol')
		{
			enrich_result$geneName_from_bioMart <- NA;
			for(i in 1:nrow(enrich_result))
			{
				if(id_type == "ensembl_gene_id")
				{
					enrich_result$geneName_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(enrich_result$geneID[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$ensembl_gene_id == x)][1]})), collapse = "/");
				}else if(id_type == "entrezgene_id")
				{
					enrich_result$geneName_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(enrich_result$geneID[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
				}
			}
		}

		## Lollipop chart to visualize rich factors -----------------
		drawLollipopPlot(data_df = enrich_result, x = "richFactor", y = "ID", color_by = "p.adjust", size_by = "Count", title = "Lollipop plot of enriched MSigDB genesets in ORA", x_axis_title = "Rich Factor", y_axis_title = "MSigDB geneset", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));	
	}
	return(enrich_result);
}

## Run ORA on GO terms -------------------------------------------
runORAonGO <- Data.frame() ? function(de_genes = Character(), ref_genes = Character(), id_type = Character(), gene_ann_df = Data.frame(), species = Character(), ontology_type = Character(), threshold = Double(), output_file_prefix = Character())
{
	declare("enrich_result", Any());
	declare("species_initials", Character());
	declare("org_genomwide_annotation_library", Character());
	declare("key_type", Character());

	species_initials <- paste(unlist(lapply(X = strsplit(species, split = "\\s"), FUN = function(x){substr(x, start = 1, stop = 1)})), collapse = "");
	org_genomwide_annotation_library <- paste("org", species_initials, "eg", "db", sep = "."); # like org.Hs.eg.db

	if(!require(org_genomwide_annotation_library, character.only = TRUE))
	{
		BiocManager::install(org_genomwide_annotation_library, dependencies = TRUE);
		library(org_genomwide_annotation_library, character.only = TRUE);
	}

	if(id_type == "ensembl_gene_id")
	{
		key_type <- "ENSEMBL";
	}else if(id_type == "entrezgene_id")
	{
		key_type <- "ENTREZID";
	}else if(id_type == "uniprot_gn_symbol")
	{
		key_type <- "SYMBOL";
	}

	# Peform ORA ---------------------------------------------------
	enrich_result <- clusterProfiler::enrichGO(gene = de_genes, universe = ref_genes, OrgDb = org_genomwide_annotation_library, keyType = key_type, ont = ontology_type, pvalueCutoff = threshold, pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10, maxGSSize = Inf, readable = FALSE, pool = FALSE);
	
	# Return an empty data.frame if no gene set have size > 10 which satisfies pvalueCutoff = threshold
	if(is.null(enrich_result))
	{
		enrich_result <- as.data.frame(enrich_result);
		return(enrich_result);
	}
	
	enrich_result <- enrich_result %>% filter(p.adjust < threshold);

	# Remove redundant terms from enrichment result ------------------
	enrich_result <- clusterProfiler::simplify(x = enrich_result, cutoff = 0.7, by = "p.adjust", select_fun = min);
		
	## Add rich factor of enriched terms -----------------------------
	# GeneRatio = k/m, BgRatio = M/N, then richFactor = k/M
	enrich_result <- clusterProfiler::mutate(enrich_result, richFactor = Count/as.numeric(sub(pattern = "/\\d+", replacement = "", x = BgRatio)));

	enrich_result <- as.data.frame(enrich_result);
	#colnames(enrich_result);	#ID Description GeneRatio BgRatio pvalue p.adjust qvalue geneID Count richFactor

	if(nrow(enrich_result))
	{
		## Add gene symbols, if ids are not the same already ----------
		if(id_type != 'uniprot_gn_symbol')
		{
			enrich_result$geneName_from_bioMart <- NA;
			
			for(i in 1:nrow(enrich_result))
			{
				if(id_type == "ensembl_gene_id")
				{
					enrich_result$geneName_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(enrich_result$geneID[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$ensembl_gene_id == x)][1]})), collapse = "/");
				}else if(id_type == "entrezgene_id")
				{
					enrich_result$geneName_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(enrich_result$geneID[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
				}
			}
		}
		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = enrich_result, x = "Count", y = "Description", fill_by = "p.adjust", title = "Bar plot of enriched terms in ORA", x_axis_title = "Term", y_axis_title = "Size of gene set after removing genes not present in input genelist", outfile = paste(output_file_prefix, ".html", sep = ""));

		## Lollipop chart to visualize rich factors -----------------
		drawLollipopPlot(data_df = enrich_result, x = "richFactor", y = "Description", color_by = "p.adjust", size_by = "Count", title = "Lollipop plot of enriched Gene Ontology terms in ORA", x_axis_title = "Rich Factor", y_axis_title = "Gene Ontology", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));	
	}

	return(enrich_result);
}

## Run ORA on DO gene sets -----------------------------------------------------
runORAonDO <- Data.frame() ? function(de_genes = Character(), ref_genes = Character(), id_type = Character(), mapped_entrez_gene_ids = Character(), gene_ann_df = Data.frame(), threshold = Double(), output_file_prefix = Character())
{
	# It requires entrez gene ids and applicable to 'Homo sapiens' only

	declare("enrich_result", Any());

	if(id_type != "entrezgene_id")
	{
		de_genes <- mapped_entrez_gene_ids[de_genes];
		ref_genes <- mapped_entrez_gene_ids[ref_genes];
		
		## In case, id_type is other than Entrez, there might be few entries with NA. We need to remove them.
		de_genes <- de_genes[!is.na(de_genes)];
		ref_genes <- ref_genes[!is.na(ref_genes)];
	}

	# Peform ORA ---------------------------------------------------
	enrich_result <- DOSE::enrichDO(gene = de_genes, universe =ref_genes, ont = "DO", pvalueCutoff = threshold, pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10, maxGSSize = Inf, readable = FALSE);
	
	# Return an empty data.frame if no gene set have size > 10 which satisfies pvalueCutoff = threshold
	if(is.null(enrich_result))
	{
		enrich_result <- as.data.frame(enrich_result);
		return(enrich_result);
	}

	enrich_result <- enrich_result %>% filter(p.adjust < threshold);

	## Add rich factor of enriched terms -----------------------------
	# GeneRatio = k/m, BgRatio = M/N, then richFactor = k/M
	enrich_result <- clusterProfiler::mutate(enrich_result, richFactor = Count/as.numeric(sub(pattern = "/\\d+", replacement = "", x = BgRatio)));

	enrich_result <- as.data.frame(enrich_result);
	#colnames(enrich_result); #ID Description GeneRatio BgRatio pvalue p.adjust qvalue geneID Count richFactor

	if(nrow(enrich_result))
	{
		## Add gene symbols, as enrichment output cotains only entrezgene_id ----------
		
		enrich_result$geneName_from_bioMart <- NA;

		for(i in 1:nrow(enrich_result))
		{
			enrich_result$geneName_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(enrich_result$geneID[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
		}
		
		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = enrich_result, x = "Count", y = "Description", fill_by = "p.adjust", title = "Bar plot of enriched terms in ORA", x_axis_title = "Term", y_axis_title = "Size of gene set after removing genes not present in input genelist", outfile = paste(output_file_prefix, ".html", sep = ""));

		## Lollipop chart to visualize rich factors -----------------
		drawLollipopPlot(data_df = enrich_result, x = "richFactor", y = "Description", color_by = "p.adjust", size_by = "Count", title = "Lollipop plot of enriched Disease Ontology terms in ORA", x_axis_title = "Rich Factor", y_axis_title = "Disease Ontology", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));
	}

	return(enrich_result);
}

## Run ORA on NCG gene sets --------------------------------------------------------
runORAonNCG <- Data.frame() ? function(de_genes = Character(), ref_genes = Character(), id_type = Character(), mapped_entrez_gene_ids = Character(), gene_ann_df = Data.frame(), threshold = Double(), output_file_prefix = Character())
{
	# It requires entrez gene ids and applicable to 'Homo sapiens' only

	declare("enrich_result", Any());

	if(id_type != "entrezgene_id")
	{
		de_genes <- mapped_entrez_gene_ids[de_genes];
		ref_genes <- mapped_entrez_gene_ids[ref_genes];
		
		## In case, id_type is other than Entrez, there might be few entries with NA. We need to remove them.
		de_genes <- de_genes[!is.na(de_genes)];
		ref_genes <- ref_genes[!is.na(ref_genes)];
	}

	# Peform ORA ---------------------------------------------------
	enrich_result <- DOSE::enrichNCG(gene = de_genes, universe = ref_genes, pvalueCutoff = threshold, pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10, maxGSSize = Inf, readable = FALSE);

	# Return an empty data.frame if no gene set have size > 10 which satisfies pvalueCutoff = threshold
	if(is.null(enrich_result))
	{
		enrich_result <- as.data.frame(enrich_result);
		return(enrich_result);
	}

	enrich_result <- enrich_result %>% filter(p.adjust < threshold);

	## Add rich factor of enriched terms -----------------------------
	# GeneRatio = k/m, BgRatio = M/N, then richFactor = k/M
	enrich_result <- clusterProfiler::mutate(enrich_result, richFactor = Count/as.numeric(sub(pattern = "/\\d+", replacement = "", x = BgRatio)));

	enrich_result <- as.data.frame(enrich_result);
	#colnames(enrich_result); #ID Description GeneRatio BgRatio pvalue p.adjust qvalue geneID Count richFactor

	#Note: The columns 'ID' and 'Description' might have an entry '-'. This denotes unknown cancer type.

	if(nrow(enrich_result))
	{
		## Add gene symbols, as enrichment output has only entrezgene_id ----------
		enrich_result$geneName_from_bioMart <- NA;
			
		for(i in 1:nrow(enrich_result))
		{
			enrich_result$geneName_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(enrich_result$geneID[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
		}
		
		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = enrich_result, x = "Count", y = "ID", fill_by = "p.adjust", title = "Bar plot of enriched terms in ORA", x_axis_title = "Term", y_axis_title = "Size of gene set after removing genes not present in input genelist", outfile = paste(output_file_prefix, ".html", sep = ""));

		## Lollipop chart to visualize rich factors -----------------
		drawLollipopPlot(data_df = enrich_result, x = "richFactor", y = "Description", color_by = "p.adjust", size_by = "Count", title = "Lollipop plot of enriched Network of Cancer Genes in ORA", x_axis_title = "Rich Factor", y_axis_title = "NCG geneset", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));
	}

	return(enrich_result);
}

## Run ORA on DisGeNET gene sets --------------------------------------------------------
runORAonDGN <- Data.frame() ? function(de_genes = Character(), ref_genes = Character(), id_type = Character(), mapped_entrez_gene_ids = Character(), gene_ann_df = Data.frame(), threshold = Double(), output_file_prefix = Character())
{
	# It requires entrez gene ids and applicable to 'Homo sapiens' only

	declare("enrich_result", Any());

	if(id_type != "entrezgene_id")
	{
		de_genes <- mapped_entrez_gene_ids[de_genes];
		ref_genes <- mapped_entrez_gene_ids[ref_genes];
		
		## In case, id_type is other than Entrez, there might be few entries with NA. We need to remove them.
		de_genes <- de_genes[!is.na(de_genes)];
		ref_genes <- ref_genes[!is.na(ref_genes)];
	}

	# Peform ORA ---------------------------------------------------
	enrich_result <- DOSE::enrichDGN(gene = de_genes, universe = ref_genes, pvalueCutoff = threshold, pAdjustMethod = "BH", qvalueCutoff = 1, minGSSize = 10, maxGSSize = Inf, readable = FALSE);
	
	# Return an empty data.frame if no gene set have size > 10 which satisfies pvalueCutoff = threshold
	if(is.null(enrich_result))
	{
		enrich_result <- as.data.frame(enrich_result);
		return(enrich_result);
	}

	enrich_result <- enrich_result %>% filter(p.adjust < threshold);

	## Add rich factor of enriched terms -----------------------------
	# GeneRatio = k/m, BgRatio = M/N, then richFactor = k/M
	enrich_result <- clusterProfiler::mutate(enrich_result, richFactor = Count/as.numeric(sub(pattern = "/\\d+", replacement = "", x = BgRatio)));

	enrich_result <- as.data.frame(enrich_result);
	#colnames(enrich_result)	# ID	Description	GeneRatio	BgRatio	pvalue	p.adjust	qvalue	geneID	Count	richFactor

	if(nrow(enrich_result))
	{
		## Add gene symbols, as enrichment output has only entrezgene_id ----------
		enrich_result$geneName_from_bioMart <- NA;
			
		for(i in 1:nrow(enrich_result))
		{
			enrich_result$geneName_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(enrich_result$geneID[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
		}
		
		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = enrich_result, x = "Count", y = "Description", fill_by = "p.adjust", title = "Bar plot of enriched terms in ORA", x_axis_title = "Term", y_axis_title = "Size of gene set after removing genes not present in input genelist", outfile = paste(output_file_prefix, ".html", sep = ""));

		## Lollipop chart to visualize rich factors -----------------
		drawLollipopPlot(data_df = enrich_result, x = "richFactor", y = "Description", color_by = "p.adjust", size_by = "Count", title = "Lollipop plot of enriched DisGeNET genesets in ORA", x_axis_title = "Rich Factor", y_axis_title = "DisGeNET geneset", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));
	}

	return(enrich_result);
}

## Run ORA on KEGG pathways ---------------------------------------------------------------
runORAonKEGG <- Data.frame() ? function(de_genes = Character(), ref_genes = Character(), id_type = Character(), mapped_entrez_gene_ids = Character(), gene_ann_df = Data.frame(), species = Character(), threshold = Double(), output_file_prefix = Character())
{
	# It requires entrez gene ids

	declare("enrich_result", Any());
	declare("organism_code", Character());
	
	organism_code <- paste(tolower(substr(unlist(strsplit(species, split = "\\s"))[1], start = 1, stop = 1)), substr(unlist(strsplit(species, split = "\\s"))[2], start = 1, stop = 2), sep = "");

	if(id_type != "entrezgene_id")
	{
		de_genes <- mapped_entrez_gene_ids[de_genes];
		ref_genes <- mapped_entrez_gene_ids[ref_genes];
		
		## In case, id_type is other than Entrez, there might be few entries with NA. We need to remove them.
		de_genes <- de_genes[!is.na(de_genes)];
		ref_genes <- ref_genes[!is.na(ref_genes)];
	}

	# Peform ORA ---------------------------------------------------
	enrich_result <- clusterProfiler::enrichKEGG(gene = de_genes, organism = organism_code, keyType = "kegg", pvalueCutoff = threshold, pAdjustMethod = "BH", qvalueCutoff = 1, universe = ref_genes, minGSSize = 10, maxGSSize = Inf, use_internal_data = FALSE); 
	# use_internal_data = TRUE will use KEGG.db generated in 2012, while FALSE will download the latest data;
	# However, KEGG.db is not available now. Therefore, this module will not be used as of now.

	# Return an empty data.frame if no gene set have size > 10 which satisfies pvalueCutoff = threshold
	if(is.null(enrich_result))
	{
		enrich_result <- as.data.frame(enrich_result);
		return(enrich_result);
	}
	
	enrich_result <- enrich_result %>% filter(p.adjust < threshold);

	## Add rich factor of enriched terms -----------------------------
	# GeneRatio = k/m, BgRatio = M/N, then richFactor = k/M
	enrich_result <- clusterProfiler::mutate(enrich_result, richFactor = Count/as.numeric(sub(pattern = "/\\d+", replacement = "", x = BgRatio)));

	enrich_result <- as.data.frame(enrich_result);
	#colnames(enrich_result); # ID	Description	GeneRatio	BgRatio	pvalue	p.adjust	qvalue	geneID	Count	richFactor

	if(nrow(enrich_result))
	{
		## Add gene symbols, as enrichment output contains only entrezgene_id ----------
		enrich_result$geneName_from_bioMart <- NA;
			
		for(i in 1:nrow(enrich_result))
		{
			enrich_result$geneName_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(enrich_result$geneID[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
		}
		
		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = enrich_result, x = "Count", y = "Description", fill_by = "p.adjust", title = "Bar plot of enriched terms in ORA", x_axis_title = "Term", y_axis_title = "Size of gene set after removing genes not present in input genelist", outfile = paste(output_file_prefix, ".html", sep = ""));

		## Lollipop chart to visualize rich factors -----------------
		drawLollipopPlot(data_df = enrich_result, x = "richFactor", y = "Description", color_by = "p.adjust", size_by = "Count", title = "Lollipop plot of enriched KEGG pathways in ORA", x_axis_title = "Rich Factor", y_axis_title = "KEGG Pathway", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));
	}

	return(enrich_result);
}

## Run ORA on Reactome pathways ---------------------------------------------------------------
runORAonReactome <- Data.frame() ? function(de_genes = Character(),  ref_genes = Character(), id_type = Character(), mapped_entrez_gene_ids = Character(), gene_ann_df = Data.frame(), species = Character(), threshold = Double(), output_file_prefix = Character())
{
	# It requires entrez gene ids

	declare("enrich_result", Any());
	declare("species_initials", Character());
	declare("organism_code", Character());
	declare("org_genomwide_annotation_library", Character());
	
	if(species == "Homo sapiens")
	{
		organism_code <- "human";
	}else if(species == "Mus musculus")
	{
		organism_code <- "mouse";
	}else if(species == "Rattus norvegicus")
	{
		organism_code <- "rat";
	}else if(species == "Caenorhabditis elegans")
	{
		organism_code <- "celegans";
	}else if(species == "Saccharomyces cerevisiae")
	{
		organism_code <- "yeast";
	}else if(species == "Danio rerio")
	{
		organism_code <- "zebrafish"
	}else if(species == "Drosophila melanogaster")
	{
		organism_code <- "fly";
	}

	species_initials <- paste(unlist(lapply(X = strsplit(species, split = "\\s"), FUN = function(x){substr(x, start = 1, stop = 1)})), collapse = "");
	org_genomwide_annotation_library <- paste("org", species_initials, "eg", "db", sep = "."); # like org.Hs.eg.db

	if(!require(org_genomwide_annotation_library, character.only = TRUE))
	{
		BiocManager::install(org_genomwide_annotation_library, dependencies = TRUE);
		library(org_genomwide_annotation_library, character.only = TRUE);
	}

	if(id_type != "entrezgene_id")
	{
		de_genes <- mapped_entrez_gene_ids[de_genes];
		ref_genes <- mapped_entrez_gene_ids[ref_genes];
		
		## In case, id_type is other than Entrez, there might be few entries with NA. We need to remove them.
		de_genes <- de_genes[!is.na(de_genes)];
		ref_genes <- ref_genes[!is.na(ref_genes)];
	}

	# Peform ORA ---------------------------------------------------
	enrich_result <- ReactomePA::enrichPathway(gene = de_genes, organism = organism_code, pvalueCutoff = threshold, pAdjustMethod = "BH", qvalueCutoff = 1, universe = ref_genes, minGSSize = 10, maxGSSize = Inf, readable = FALSE); 

	# Return an empty data.frame if no gene set have size > 10 which satisfies pvalueCutoff = threshold
	if(is.null(enrich_result))
	{
		enrich_result <- as.data.frame(enrich_result);
		return(enrich_result);
	}

	enrich_result <- enrich_result %>% filter(p.adjust < threshold);

	## Add rich factor of enriched terms -----------------------------
	# GeneRatio = k/m, BgRatio = M/N, then richFactor = k/M
	enrich_result <- clusterProfiler::mutate(enrich_result, richFactor = Count/as.numeric(sub(pattern = "/\\d+", replacement = "", x = BgRatio)));

	enrich_result <- as.data.frame(enrich_result);
	#colnames(enrich_result);	#ID	Description"	GeneRatio	BgRatio	pvalue	p.adjust	qvalue	geneID	Count	richFactor

	if(nrow(enrich_result))
	{
		## Add gene symbols as enrichment result contains entrezgene_id only  ----------
		enrich_result$geneName_from_bioMart <- NA;
			
		for(i in 1:nrow(enrich_result))
		{
			enrich_result$geneName_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(enrich_result$geneID[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
		}
		
		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = enrich_result, x = "Count", y = "Description", fill_by = "p.adjust", title = "Bar plot of enriched terms in ORA", x_axis_title = "Term", y_axis_title = "Size of gene set after removing genes not present in input genelist", outfile = paste(output_file_prefix, ".html", sep = ""));
		

		## Lollipop chart to visualize rich factors -----------------
		drawLollipopPlot(data_df = enrich_result, x = "richFactor", y = "Description", color_by = "p.adjust", size_by = "Count", title = "Lollipop plot of enriched Reactome Pathways in ORA", x_axis_title = "Rich Factor", y_axis_title = "Reactome Pathway", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));
	}

	return(enrich_result);
}

## Run GSEA on MSigDB gene sets -------------------------------------------------------------
runGSEAonMSigDB <- Data.frame() ? function(de_genes = Double(), geneset_list = typed::List(), id_type = Character(), gene_ann_df = Data.frame(), threshold = Double(), output_file_prefix = Character())
{
	declare("gse_result", Data.frame());

	# Peform GSEA ---------------------------------------------------
	gse_result <- fgsea::fgsea(pathways = geneset_list, stats = de_genes, minSize = 10, maxSize = Inf);
	
	# Return an empty data.frame if no gene set have size > 10 
	if(is.null(gse_result))
	{
		gse_result <- as.data.frame(gse_result);
		return(gse_result);
	}

	gse_result <- gse_result[gse_result$padj < threshold,];
	
	gse_result <- as.data.frame(gse_result);
	#colnames(gse_result);	#pathway pval padj log2err ES NES size leadingEdge

	if(nrow(gse_result))
	{
		# gse_result$leadingEdge is a list. We need some data manipulation here.
		declare("leading_egde_lists", typed::List(), value = gse_result$leadingEdge);

		# Remove leadingEdge column from gse_result, and add new which will now hold characters
		gse_result$leadingEdge <- NULL;	

		for(i in 1:nrow(gse_result))
		{
			gse_result$leadingEdge[i] <- paste(leading_egde_lists[[i]], collapse = "/");
		}

		## Add gene symbols, if ids are not the same already ----------
		if(id_type != 'uniprot_gn_symbol')
		{
			gse_result$leadingEdge_geneName_from_bioMart <- NA;
			
			for(i in 1:nrow(gse_result))
			{
				if(id_type == "ensembl_gene_id")
				{
					gse_result$leadingEdge_geneName_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(gse_result$leadingEdge[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$ensembl_gene_id == x)][1]})), collapse = "/");
				}else if(id_type == "entrezgene_id")
				{
					gse_result$leadingEdge_geneName_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(gse_result$leadingEdge[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
				}
			}
		}
		
		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = gse_result, x = "NES", y = "pathway", fill_by = "padj", title = "Bar plot of enriched terms in GSEA", x_axis_title = "Normalized Enrichment Score", y_axis_title = "Pathway", outfile = paste(output_file_prefix, ".html", sep = ""));

		## Lollipop chart to visualize rich factors -----------------
		drawLollipopPlot(data_df = gse_result, x = "NES", y = "pathway", color_by = "padj", size_by = "size", title = "Lollipop plot of enriched MSigDB genesets in GSEA", x_axis_title = "Normalized Enrichment Score", y_axis_title = "MSigDB geneset", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));

		
	}

	return(gse_result);
}

## Run GSEA on GO terms -----------------------------------------------------------------
runGSEAonGO <- Data.frame() ? function(de_genes = Double(), id_type = Character(), gene_ann_df = Data.frame(), species = Character(), ontology_type = character(), threshold = Double(), output_file_prefix = Character())
{
	declare("gse_result", Any());
	declare("species_initials", Character());
	declare("org_genomwide_annotation_library", Character());

	species_initials <- paste(unlist(lapply(X = strsplit(species, split = "\\s"), FUN = function(x){substr(x, start = 1, stop = 1)})), collapse = "");
	org_genomwide_annotation_library <- paste("org", species_initials, "eg", "db", sep = "."); # like org.Hs.eg.db

	if(!require(org_genomwide_annotation_library, character.only = TRUE))
	{
		BiocManager::install(org_genomwide_annotation_library, dependencies = TRUE);
		library(org_genomwide_annotation_library, character.only = TRUE);
	}

	if(id_type == "ensembl_gene_id")
	{
		key_type <- "ENSEMBL";
	}else if(id_type == "entrezgene_id")
	{
		key_type <- "ENTREZID";
	}else if(id_type == "uniprot_gn_symbol")
	{
		key_type <- "SYMBOL";
	}

	# Peform GSEA ---------------------------------------------------
	gse_result <- clusterProfiler::gseGO(geneList = de_genes, OrgDb = get(org_genomwide_annotation_library), keyType = key_type, ont = ontology_type, exponent = 1, nPerm = 10000, minGSSize = 10, maxGSSize = Inf, pvalueCutoff = threshold, pAdjustMethod = "BH", seed = TRUE, by = "DOSE");
	gse_result <- gse_result %>% filter(p.adjust < threshold);

	# Remove redundant terms from enrichment result ------------------
	gse_result <- clusterProfiler::simplify(x = gse_result, cutoff = 0.7, by = "p.adjust", select_fun = min);

	gse_result <- as.data.frame(gse_result);
	#colnames(gse_result);	#ID Description setSize enrichmentScore NES pvalue p.adjust qvalues rank leading_edge core_enrichment

	if(nrow(gse_result))
	{
		## Add gene symbols, if ids are not the same already ----------
		if(id_type != 'uniprot_gn_symbol')
		{
			gse_result$core_enrichment_genename_from_bioMart <- NA;
			
			for(i in 1:nrow(gse_result))
			{
				if(id_type == "ensembl_gene_id")
				{
					gse_result$core_enrichment_genename_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(gse_result$core_enrichment[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$ensembl_gene_id == x)][1]})), collapse = "/");
				}else if(id_type == "entrezgene_id")
				{
					gse_result$core_enrichment_genename_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(gse_result$core_enrichment[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
				}
			}
		}

		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = gse_result, x = "setSize", y = "Description", fill_by = "p.adjust", title = "Bar plot of enriched terms in GSEA", x_axis_title = "Term", y_axis_title = "Size of gene set after removing genes not present in input genelist", outfile = paste(output_file_prefix, ".html", sep = ""));

		## Lollipop chart to visualize Normalized Enrichment Score (NES) -----------------
		drawLollipopPlot(data_df = gse_result, x = "NES", y = "Description", color_by = "p.adjust", size_by = "setSize", title = "Lollipop plot of enriched Gene Ontology terms in GSEA", x_axis_title = "Normalized Enrichment Score", y_axis_title = "Gene Ontology", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));

		## UpSet plot to visualize distribution of genes across various enriched terms --------------
		#declare("temp_gene_term_map", typed::List(null_ok = TRUE));
		#declare("temp_data_df_for_upset", Data.frame());
		
		#for(i in 1:nrow(gse_result))
		#{
		#	invisible(lapply(X = unlist(strsplit(gse_result$core_enrichment[i], split = "/")), FUN = function(x){ temp_gene_term_map[[x]][length(temp_gene_term_map[[x]])+1] <<- gse_result$Description[i] }));
		#}

		#temp_data_df_for_upset <- data.frame(Gene = names(temp_gene_term_map), Description = I(temp_gene_term_map));
		#head(temp_data_df_for_upset)
        #Gene  Description
		#MMP1     MMP1 immune s....
		#S100A9 S100A9 immune s....


		#drawUpSetPlot(data_df = gse_result, x = "Description", title = "UpSet plot on enriched terms", outfile = paste(output_file_prefix, "_upset_plot", ".html", sep = ""));
	}

	return(gse_result);
}

## Run GSEA on DO gene sets -----------------------------------------------------
runGSEAonDO <- Data.frame() ? function(de_genes = Double(), id_type = Character(), mapped_entrez_gene_ids = Character(), gene_ann_df = Data.frame(), threshold = Double(), output_file_prefix = Character())
{
	# It requires entrez gene ids and applicable to 'Homo sapiens' only

	declare("gse_result", Any());

	if(id_type != "entrezgene_id")
	{
		names(de_genes) <- mapped_entrez_gene_ids[names(de_genes)];
				
		## In case, id_type is other than Entrez, there might be few entries with NA. We need to remove them.
		de_genes <- de_genes[!is.na(names(de_genes))];

		if(length(de_genes) == 0)
		{
			gse_result <- as.data.frame(gse_result);
			return(gse_result);
		}
	}

	# Peform GSEA ---------------------------------------------------
	gse_result <- DOSE::gseDO(geneList = de_genes, exponent = 1, minGSSize = 10, maxGSSize = Inf, pvalueCutoff = threshold, pAdjustMethod = "BH", seed = FALSE, by = "fgsea");
	gse_result <- gse_result %>% filter(p.adjust < threshold);

	gse_result <- as.data.frame(gse_result);
	#colnames(gse_result);	# ID	Description	setSize	enrichmentScore	NES	pvalue	p.adjust	qvalues	rank	leading_edge	core_enrichment

	if(nrow(gse_result))
	{
		## Add gene symbols, as enrichment output cotains only entrezgene_id ----------

		gse_result$core_enrichment_genename_from_bioMart <- NA;
		
		for(i in 1:nrow(gse_result))
		{
			gse_result$core_enrichment_genename_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(gse_result$core_enrichment[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
		}

		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = gse_result, x = "setSize", y = "Description", fill_by = "p.adjust", title = "Bar plot of enriched terms in GSEA", x_axis_title = "Term", y_axis_title = "Size of gene set after removing genes not present in input genelist", outfile = paste(output_file_prefix, ".html", sep = ""));

		## Lollipop chart to visualize Normalized Enrichment Score (NES) -----------------
		drawLollipopPlot(data_df = gse_result, x = "NES", y = "Description", color_by = "p.adjust", size_by = "setSize", title = "Lollipop plot of enriched Disease Ontology terms in GSEA", x_axis_title = "Normalized Enrichment Score", y_axis_title = "Disease Ontology", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));

	}

	return(as.data.frame(gse_result));
}

## Run GSEA on NCG gene sets -----------------------------------------------------
runGSEAonNCG <- Data.frame() ? function(de_genes = Double(), id_type = Character(), mapped_entrez_gene_ids = Character(), gene_ann_df = Data.frame(), threshold = Double(), output_file_prefix = Character())
{
	# It requires entrez gene ids and applicable to 'Homo sapiens' only

	declare("gse_result", Any());

	if(id_type != "entrezgene_id")
	{
		names(de_genes) <- mapped_entrez_gene_ids[names(de_genes)];
				
		## In case, id_type is other than Entrez, there might be few entries with NA. We need to remove them.
		de_genes <- de_genes[!is.na(names(de_genes))];

		if(length(de_genes) == 0)
		{
			gse_result <- as.data.frame(gse_result);
			return(gse_result);
		}
	}

	# Peform GSEA ---------------------------------------------------
	gse_result <- DOSE::gseNCG(geneList = de_genes, exponent = 1, minGSSize = 10, maxGSSize = Inf, pvalueCutoff = threshold, pAdjustMethod = "BH", seed = FALSE, by = "fgsea");
	gse_result <- gse_result %>% filter(p.adjust < threshold);

	gse_result <- as.data.frame(gse_result);
	#colnames(gse_result);	#ID Description setSize enrichmentScore NES pvalue p.adjust qvalues rank leading_edge core_enrichment

	if(nrow(gse_result))
	{
		## Add gene symbols, as enrichment output cotains only entrezgene_id ----------

		gse_result$core_enrichment_genename_from_bioMart <- NA;
		
		for(i in 1:nrow(gse_result))
		{
			gse_result$core_enrichment_genename_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(gse_result$core_enrichment[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
		}

		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = gse_result, x = "setSize", y = "ID", fill_by = "p.adjust", title = "Bar plot of enriched terms in GSEA", x_axis_title = "Term", y_axis_title = "Size of gene set after removing genes not present in input genelist", outfile = paste(output_file_prefix, ".html", sep = ""));

		## Lollipop chart to visualize Normalized Enrichment Score (NES) -----------------
		drawLollipopPlot(data_df = gse_result, x = "NES", y = "Description", color_by = "p.adjust", size_by = "setSize", title = "Lollipop plot of enriched Network of Cancer Genes genesets in GSEA", x_axis_title = "Normalized Enrichment Score", y_axis_title = "NCG geneset", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));
	}

	return(gse_result);
}

## Run GSEA on DGN gene sets -----------------------------------------------------
runGSEAonDGN <- Data.frame() ? function(de_genes = Double(), id_type = Character(), mapped_entrez_gene_ids = Character(), gene_ann_df = Data.frame(), threshold = Double(), output_file_prefix = Character())
{
	# It requires entrez gene ids and applicable to 'Homo sapiens' only

	declare("gse_result", Any());

	if(id_type != "entrezgene_id")
	{
		names(de_genes) <- mapped_entrez_gene_ids[names(de_genes)];
				
		## In case, id_type is other than Entrez, there might be few entries with NA. We need to remove them.
		de_genes <- de_genes[!is.na(names(de_genes))];

		if(length(de_genes) == 0)
		{
			gse_result <- as.data.frame(gse_result);
			return(gse_result);
		}
	}

	# Peform GSEA ---------------------------------------------------
	gse_result <- DOSE::gseDGN(geneList = de_genes, exponent = 1, minGSSize = 10, maxGSSize = Inf, pvalueCutoff = threshold, pAdjustMethod = "BH", seed = FALSE, by = "fgsea");
	gse_result <- gse_result %>% filter(p.adjust < threshold);
	
	gse_result <- as.data.frame(gse_result);
	#colnames(gse_result)	#ID	Description	setSize	enrichmentScore	NES	pvalue	p.adjust	qvalues	rank	leading_edge	core_enrichment

	if(nrow(gse_result))
	{
		## Add gene symbols, as enrichment output cotains only entrezgene_id ----------

		gse_result$core_enrichment_genename_from_bioMart <- NA;
		
		for(i in 1:nrow(gse_result))
		{
			gse_result$core_enrichment_genename_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(gse_result$core_enrichment[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
		}

		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = gse_result, x = "setSize", y = "Description", fill_by = "p.adjust", title = "Bar plot of enriched terms in GSEA", x_axis_title = "Term", y_axis_title = "Size of gene set after removing genes not present in input genelist", outfile = paste(output_file_prefix, ".html", sep = ""));

		## Lollipop chart to visualize Normalized Enrichment Score (NES) -----------------
		drawLollipopPlot(data_df = gse_result, x = "NES", y = "Description", color_by = "p.adjust", size_by = "setSize", title = "Lollipop plot of enriched DisGeNET genesets in GSEA", x_axis_title = "Normalized Enrichment Score", y_axis_title = "DisGeNET geneset", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));
	}

	return(gse_result);
}

## Run GSEA on KEGG pathways ---------------------------------------------------------------
runGSEAonKEGG <- Data.frame() ? function(de_genes = Double(), id_type = Character(), mapped_entrez_gene_ids = Character(), gene_ann_df = Data.frame(), species = Character(), threshold = Double(), output_file_prefix = Character())
{
	# It requires entrez gene ids

	declare("gse_result", Any());
	declare("organism_code", Character());
	
	organism_code <- paste(tolower(substr(unlist(strsplit(species, split = "\\s"))[1], start = 1, stop = 1)), substr(unlist(strsplit(species, split = "\\s"))[2], start = 1, stop = 2), sep = "");

	if(id_type != "entrezgene_id")
	{
		names(de_genes) <- mapped_entrez_gene_ids[names(de_genes)];
				
		## In case, id_type is other than Entrez, there might be few entries with NA. We need to remove them.
		de_genes <- de_genes[!is.na(names(de_genes))];

		if(length(de_genes) == 0)
		{
			gse_result <- as.data.frame(gse_result);
			return(gse_result);
		}
	}

	# Peform GSEA ---------------------------------------------------
	gse_result <- clusterProfiler::gseKEGG(geneList = de_genes, organism = organism_code, keyType = "kegg", exponent = 1, minGSSize = 10, maxGSSize = Inf, pvalueCutoff = threshold, pAdjustMethod = "BH", use_internal_data = FALSE, seed = FALSE, by = "fgsea"); 
	# use_internal_data = TRUE will use KEGG.db generated in 2012, while FALSE will download the latest data;
	# However, KEGG.db is not available now. Therefore, this module will not be used as of now.
	gse_result <- gse_result %>% filter(p.adjust < threshold);

	gse_result <- as.data.frame(gse_result);
	#colnames(gse_result)	#ID	Description	setSize	enrichmentScore	NES	pvalue	p.adjust	qvalues	rank	leading_edge	core_enrichment

	if(nrow(gse_result))
	{
		## Add gene symbols, as enrichment output cotains only entrezgene_id ----------

		gse_result$core_enrichment_genename_from_bioMart <- NA;
		
		for(i in 1:nrow(gse_result))
		{
			gse_result$core_enrichment_genename_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(gse_result$core_enrichment[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
		}

		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = gse_result, x = "setSize", y = "Description", fill_by = "p.adjust", title = "Bar plot of enriched terms in GSEA", x_axis_title = "Term", y_axis_title = "Size of gene set after removing genes not present in input genelist", outfile = paste(output_file_prefix, ".html", sep = ""));

		## Lollipop chart to visualize Normalized Enrichment Score (NES) -----------------
		drawLollipopPlot(data_df = gse_result, x = "NES", y = "Description", color_by = "p.adjust", size_by = "setSize", title = "Lollipop plot of enriched KEGG pathways in GSEA", x_axis_title = "Normalized Enrichment Score", y_axis_title = "KEGG pathway", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));
	}

	return(gse_result);
}

## Run GSEA on Reactome pathways ---------------------------------------------------------------
runGSEAonReactome <- Data.frame() ? function(de_genes = Double(), id_type = Character(), mapped_entrez_gene_ids = Character(), gene_ann_df = Data.frame(), species = Character(), threshold = Double(), output_file_prefix = Character())
{
	# It requires entrez gene ids

	declare("gse_result", Any());
	declare("species_initials", Character());
	declare("organism_code", Character());
	declare("org_genomwide_annotation_library", Character());
	
	if(species == "Homo sapiens")
	{
		organism_code <- "human";
	}else if(species == "Mus musculus")
	{
		organism_code <- "mouse";
	}else if(species == "Rattus norvegicus")
	{
		organism_code <- "rat";
	}else if(species == "Caenorhabditis elegans")
	{
		organism_code <- "celegans";
	}else if(species == "Saccharomyces cerevisiae")
	{
		organism_code <- "yeast";
	}else if(species == "Danio rerio")
	{
		organism_code <- "zebrafish"
	}else if(species == "Drosophila melanogaster")
	{
		organism_code <- "fly";
	}
	
	species_initials <- paste(unlist(lapply(X = strsplit(species, split = "\\s"), FUN = function(x){substr(x, start = 1, stop = 1)})), collapse = "");
	org_genomwide_annotation_library <- paste("org", species_initials, "eg", "db", sep = "."); # like org.Hs.eg.db

	if(!require(org_genomwide_annotation_library, character.only = TRUE))
	{
		BiocManager::install(org_genomwide_annotation_library, dependencies = TRUE);
		library(org_genomwide_annotation_library, character.only = TRUE);
	}

	if(id_type != "entrezgene_id")
	{
		names(de_genes) <- mapped_entrez_gene_ids[names(de_genes)];
				
		## In case, id_type is other than Entrez, there might be few entries with NA. We need to remove them.
		de_genes <- de_genes[!is.na(names(de_genes))];

		if(length(de_genes) == 0)
		{
			gse_result <- as.data.frame(gse_result);
			return(gse_result);
		}
	}

	# Peform GSEA ---------------------------------------------------
	gse_result <- ReactomePA::gsePathway(geneList = de_genes, organism = organism_code, exponent = 1, minGSSize = 10, maxGSSize = Inf, eps = 1e-10, pvalueCutoff = threshold, pAdjustMethod = "BH", seed = FALSE, by = "fgsea"); 

	# Return an empty data.frame if no gene set have size > 10 
	if(is.null(gse_result))
	{
		gse_result <- as.data.frame(gse_result);
		return(gse_result);
	}

	gse_result <- gse_result %>% filter(p.adjust < threshold);

	gse_result <- as.data.frame(gse_result);
	#colnames(gse_result);	#ID	Description	setSize	enrichmentScore	NES	pvalue	p.adjust	qvalues	rank	leading_edge	core_enrichment

	if(nrow(gse_result))
	{
		## Add gene symbols, as enrichment output cotains only entrezgene_id ----------

		gse_result$core_enrichment_genename_from_bioMart <- NA;
		
		for(i in 1:nrow(gse_result))
		{
			gse_result$core_enrichment_genename_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(gse_result$core_enrichment[i], split = "/")), FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
		}

		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = gse_result, x = "setSize", y = "Description", fill_by = "p.adjust", title = "Bar plot of enriched terms in GSEA", x_axis_title = "Term", y_axis_title = "Size of gene set after removing genes not present in input genelist", outfile = paste(output_file_prefix, ".html", sep = ""));

		## Lollipop chart to visualize Normalized Enrichment Score (NES) -----------------
		drawLollipopPlot(data_df = gse_result, x = "NES", y = "Description", color_by = "p.adjust", size_by = "setSize", title = "Lollipop plot of enriched Reactome pathways in GSEA", x_axis_title = "Normalized Enrichment Score", y_axis_title = "Reactome pathway", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));
	}

	return(gse_result);
}

## Run SPIA ----------------------------------------------------
runSPIA <- Data.frame() ? function(de_genes = Double(), ref_genes = Character(), id_type = Character(), mapped_entrez_gene_ids = Character(), gene_ann_df = Data.frame(), species = Character(), threshold = Double(), output_file_prefix = Character())
{
	declare("organism", Character());
	declare("spia_df", Data.frame(null_ok = TRUE));

	if(species == "Homo sapiens")
	{
		organism <- "hsa";
	}else if(species == "Mus musculus")
	{
		organism <- "mmu";
	}
	
	if(id_type != "entrezgene_id")
	{
		names(de_genes) <- mapped_entrez_gene_ids[names(de_genes)];
		ref_genes <- mapped_entrez_gene_ids[ref_genes];
				
		## In case, id_type is other than Entrez, there might be few entries with NA. We need to remove them.
		de_genes <- de_genes[!is.na(names(de_genes))];
		ref_genes <- ref_genes[!is.na(ref_genes)];

		if(length(de_genes) == 0 || length(ref_genes) == 0)
		{
			return(spia_df);
		}
	}

	# Peform topology based pathway enrichment analysis -------------------------------
	spia_df <- SPIA::spia(de = de_genes, all = ref_genes, organism = organism, nB = 2000, plots = FALSE);
	#colnames(spia_df);	#Name	ID	pSize	NDE	pNDE	tA	pPERT	pG	pGFdr	pGFWER	Status	KEGGLINK

	spia_df <- spia_df[spia_df$pGFdr < threshold,];

	## Add rich factor of enriched pathways -----------------------------
	# richFactor = NDE/pSize
	spia_df$richFactor <- spia_df$NDE/spia_df$pSize;

	if(nrow(spia_df))
	{
		## Add gene symbols, as enrichment output cotains only entrezgene_id ----------

		spia_df$genename_from_bioMart <- NA;
		
		for(i in 1:nrow(spia_df))
		{
			spia_df$genename_from_bioMart[i] <- paste(unlist(lapply(X = unlist(strsplit(spia_df$KEGGLINK[i], split = "\\+"))[-1], FUN = function(x){gene_ann_df$external_gene_name[which(gene_ann_df$entrezgene_id == x)][1]})), collapse = "/");
		}

		## Generate plots on enrichment data for interpretation ----------
		#drawBarPlot(data_df = spia_df, x = "NDE", y = "Name", fill_by = "pGFdr", title = "Bar plot of enriched terms using topology based method", x_axis_title = "Term", y_axis_title = "Observed number of genes (from input list) on given pathway", outfile = paste(output_file_prefix, ".html", sep = ""));

		drawLollipopPlot(data_df = spia_df, x = "richFactor", y = "Name", color_by = "pGFdr", size_by = "NDE", title = "Lollipop plot of enriched KEGG pathways in topology based analysis", x_axis_title = "Rich Factor", y_axis_title = "KEGG pathway", outfile = paste(output_file_prefix, "_lollipop_plot", ".html", sep = ""));
	}
	return(spia_df);
}


## =============================================================
## Read input list of genes and reference list
## =============================================================
tryCatch(
	expr = 
	{
		## Check read access to input file ---------------------
		if(file.access(names = genelist_file, mode = 4) != 0)
		{
			error(file_logger, paste("Unable to read --genelist_file file: ", genelist_file, sep = ""));
			base::sink();
			cat("Unable to read --genelist_file file: ", genelist_file, "\n", sep = "");
			cat("The program: Functional_enrichment.R failed while reading input list of genes. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}

		## Read gene list --------------------------------------
		genelist_df <- as.data.frame(data.table::fread(input = genelist_file, sep = "auto", header = TRUE, check.names = FALSE), stringsAsFactors = FALSE);
		debug(file_logger, paste("Input gene list file has #rows: ", nrow(genelist_df), ", #cols: ", ncol(genelist_df), sep = ""));

		## Check validity of user specified parameters ---------
		if(!(gene_id_column %in% 1:ncol(genelist_df)))
		{
			error(file_logger, paste("User provied --gene_id_column: ", gene_id_column, " does not exist in --genelist_file ", genelist_file, sep = ""));
			base::sink();
			cat("User provied --gene_id_column: ", gene_id_column, " does not exist in --genelist_file ", genelist_file, "\n", sep = "");
			cat("The program: Functional_enrichment.R failed while reading input list of genes. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}
		if( (!is.null(log2fc_column)) && !(log2fc_column %in% 1:ncol(genelist_df)))
		{
			error(file_logger, paste("User provied --log2fc_column: ", log2fc_column, " does not exist in --genelist_file ", genelist_file, sep = ""));
			base::sink();
			cat("User provied --log2fc_column: ", log2fc_column, " does not exist in --genelist_file ", genelist_file, "\n", sep = "");
			cat("The program: Functional_enrichment.R failed while reading input list of genes. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}

		## Store genelist -------------------------------------
		genelist <- extractGenes(genelist_df, gene_id_column, log2fc_column);

		## Check read access to reference input file ---------
		if(file.access(names = reference_genelist_file, mode = 4) != 0)
		{
			error(file_logger, paste("Unable to read --reference_genelist_file: ", reference_genelist_file, sep = ""));
			base::sink();
			cat("Unable to read --reference_genelist_file: ", reference_genelist_file, "\n", sep = "");
			cat("The program: Functional_enrichment.R failed while reading refence list of genes. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}
		
		## Read reference gene list ---------------------------
		declare("reference_genelist_df", Data.frame());
		reference_genelist_df <- as.data.frame(data.table::fread(input = reference_genelist_file, sep = "auto", header = TRUE, check.names = FALSE), stringsAsFactors = FALSE);
		debug(file_logger, paste("Input reference gene list file has #rows: ", nrow(reference_genelist_df), ", #cols: ", ncol(reference_genelist_df), sep = ""));

		## Store reference gene list -------------------------
		if( !is.null(ref_gene_id_column) && !(ref_gene_id_column %in% 1:ncol(reference_genelist_df)) )
		{
			error(file_logger, paste("User provied --ref_gene_id_column: ", ref_gene_id_column, " does not exist in --reference_genelist_file ", reference_genelist_file, sep = ""));
			base::sink();
			cat("User provied --ref_gene_id_column: ", ref_gene_id_column, " does not exist in --reference_genelist_file ", reference_genelist_file, "\n", sep = "");
			cat("The program: Functional_enrichment.R failed while reading reference gene list. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}else if( !is.null(ref_gene_id_column) && (ref_gene_id_column %in% 1:ncol(reference_genelist_df)) )
		{
			ref_genelist <- as.character(reference_genelist_df[,ref_gene_id_column]);
		}
	},
	warning = function(w)
	{
		warn(file_logger, "Caught a warning while reading input list of genes. Something unexpected must have happened.");
		print(w);
		base::sink();
		print(w);
		cat("The program: Functional_enrichment.R failed while reading input list(s) of genes. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	},
	error = function(e)
	{
		error(file_logger, "An error occured while reading input list of genes.");
		print(e);
		base::sink();
		print(e);
		cat("The program: Functional_enrichment.R failed while reading input list(s) of genes. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	}
)
save(list=ls(), file = paste(output_dir, "/", "FunctionalEnrichment.RData", sep = ""));


## =====================================================================================
## Get annotations (Ensembl gene id, Entrez gene id etc.) for the all genes under study
## =====================================================================================
tryCatch(
	expr =
	{
		declare("all_genes", Character());
		declare("ensembl", Any());
		declare("ensembl_version_df", Data.frame());
		declare("gene_id_type", Character(null_ok = TRUE));
		
		all_genes <- ref_genelist;
		
		all_genes <- c(all_genes, genelist$genes);
		
		all_genes <- unique(all_genes);
		debug(file_logger, paste("Total #genes in the current study: ", length(all_genes), sep = ""));

		#gene_id_type <- identifyGeneIdType(ids = sample(x = all_genes, size = 10, replace = FALSE));
		gene_id_type <- identifyGeneIdType(ids = all_genes);
		info(file_logger, paste("Inferred gene identifier type in the current input list of genes is: ", gene_id_type, sep = ""));

		if(gene_id_type == "Unknown")
		{
			save(list=ls(), file = paste(output_dir, "/", "FunctionalEnrichment.RData", sep = ""));
			error(file_logger, "Unable to identify gene id type of input genes. Make sure that all ids are one of ensembl_gene_id (e.g., ENSG00000196611), entrezgene_id (e.g., 4312) or uniprot_gn_symbol (e.g., MMP1). Mixed types not allowed.");
			base::sink();
			cat("Unable to identify gene id type of input genes. Make sure that all ids are one of ensembl_gene_id (e.g., ENSG00000196611), entrezgene_id (e.g., 4312) or uniprot_gn_symbol (e.g., MMP1). Mixed types not allowed.", "\n", sep = "");
			cat("The program: Functional_enrichment.R failed while identifying type of gene ids. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}
		
		## Remove version information from all gene ids (in case of Ensembl, e.g., ENSG00000223972.5) --------
		if(gene_id_type == "ensembl_gene_id")
		{
			ref_genelist <-  gsub(ref_genelist, pattern = "\\..*", replacement = "");
			genelist$genes <- gsub(genelist$genes, pattern = "\\..*", replacement = "");
			all_genes <- gsub(all_genes, pattern = "\\..*", replacement = "");
		}

		## Get annotation of all genes using bioMart ---------------------------
		info(file_logger, "Preparing annotation table for all genes using bioMart ...");
		gene_info_df <- annotateWithBioMart(biomart_database_name = "ENSEMBL_MART_ENSEMBL", species = species, gene_id_type, gene_ids = all_genes, file_logger);
		if(nrow(gene_info_df) == 0)
		{
			save(list=ls(), file = paste(output_dir, "/", "FunctionalEnrichment.RData", sep = ""));
			error(file_logger, "Unable to annotate the input list of genes.");
			base::sink();
			cat("Unable to annotate the input list of genes.", "\n", sep = "");
			cat("The program: Functional_enrichment.R failed while annotating gene ids for all genes under study using bioMart. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}

		info(file_logger, "Summary on gene annotation as obtained by biomaRt:");
		cat("#Ensembl Gene Ids: ", length(unique(gene_info_df$ensembl_gene_id[! (is.na(gene_info_df$ensembl_gene_id) | (gene_info_df$ensembl_gene_id == "")) ])), "\n");
		cat("#Entrez Gene Ids: ", length(unique(gene_info_df$entrezgene_id[! (is.na(gene_info_df$entrezgene_id) | (gene_info_df$entrezgene_id == "")) ])), "\n");
		cat("#Uniprot Gene Symbols: ", length(unique(gene_info_df$uniprot_gn_symbol[! (is.na(gene_info_df$uniprot_gn_symbol) | (gene_info_df$uniprot_gn_symbol == "")) ])), "\n");
				
		## There could be scenarios where there is no 1-1 mapping between Ensembl and Entrez.
		info(file_logger, paste("Gene annotation table has ", nrow(gene_info_df), " entries.", sep = ""));
		
		write.table(x = gene_info_df, file = paste(output_dir, "/", output_prefix, "gene_info_from_BioMart.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);

		## Create a mapping for all genes to Entrez gene ids --------------------
		#if(gene_id_type != "entrezgene_id" && (any(c("DO", "NCG", "DisGeNET") %in% data_sources) || run_SPIA))
		if(gene_id_type != "entrezgene_id" && (any(c("DO", "NCG") %in% data_sources) || run_SPIA))
		{
			info(file_logger, "Some of the enrichment modules such as Disease Ontology (DO), Network of Cancer Genes (NCG), and SPIA requires entrez gene identifiers. Therefore, we will map the current identifier to entrez gene id using biomaRt.");
		}

		if(gene_id_type == "entrezgene_id")
		{
			mapped_entrez_gene_ids <- all_genes;
			names(mapped_entrez_gene_ids) <- all_genes;
		}else if(gene_id_type == "ensembl_gene_id")
		{
			mapped_entrez_gene_ids <- as.character(gene_info_df$entrezgene_id);
			names(mapped_entrez_gene_ids) <- gene_info_df$ensembl_gene_id;
			
			#Removing duplicate or NA or "" entries in 'entrezgene_id'
			mapped_entrez_gene_ids <- mapped_entrez_gene_ids[! (is.na(mapped_entrez_gene_ids) | (mapped_entrez_gene_ids == "")) ];
			mapped_entrez_gene_ids <- mapped_entrez_gene_ids[!duplicated(mapped_entrez_gene_ids)];

			#There might still be duplicate ensembl gene ids in names(mapped_entrez_gene_ids)
			mapped_entrez_gene_ids <- mapped_entrez_gene_ids[!duplicated(names(mapped_entrez_gene_ids))];

			info(file_logger, paste("The 1-1 ensembl_gene_id-entrezgene_id mapping has ", length(mapped_entrez_gene_ids), " entries", sep = ""));
		}else if(gene_id_type == "uniprot_gn_symbol")
		{
			mapped_entrez_gene_ids <- as.character(gene_info_df$entrezgene_id);
			names(mapped_entrez_gene_ids) <- gene_info_df$uniprot_gn_symbol;

			#Removing duplicate or NA or "" entries in 'entrezgene_id'
			mapped_entrez_gene_ids <- mapped_entrez_gene_ids[! (is.na(mapped_entrez_gene_ids) | (mapped_entrez_gene_ids == "")) ];
			mapped_entrez_gene_ids <- mapped_entrez_gene_ids[!duplicated(mapped_entrez_gene_ids)];

			#There might still be duplicate uniprot gene symbols in names(mapped_entrez_gene_ids)
			mapped_entrez_gene_ids <- mapped_entrez_gene_ids[!duplicated(names(mapped_entrez_gene_ids))];

			info(file_logger, paste("The 1-1 uniprot_gn_symbol-entrezgene_id mapping has ", length(mapped_entrez_gene_ids), " entries", sep = ""));
		}else
		{
			save(list=ls(), file = paste(output_dir, "/", "FunctionalEnrichment.RData", sep = ""));
			error(file_logger, "Unsupported gene identifier encountered! Currently only supported gene identifiers are: ensembl_gene_id, entrezgene_id and uniprot_gn_symbol. Therefore, program will terminate here.");
			base::sink();
			cat("Unsupported gene identifier encountered! Currently only supported gene identifiers are: ensembl_gene_id, entrezgene_id and uniprot_gn_symbol. Therefore, program will terminate here.", "\n", sep = "");
			cat("The program: Functional_enrichment.R failed while annotating gene ids. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}
	},
	warn = function(w)
	{
		warn(file_logger, "Caught a warning while preparing annotation data for the genes. Something unexpected must have happened.");
		print(w);
		base::sink();
		print(w);
		cat("The program: Functional_enrichment.R failed while annotating gene ids. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	},
	error = function(e)
	{
		error(file_logger, "An error occured while preparing annotation data for the genes.");
		print(e);
		base::sink();
		print(e);
		cat("The program: Functional_enrichment.R failed while annotating gene ids. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	}
)
save(list=ls(), file = paste(output_dir, "/", "FunctionalEnrichment.RData", sep = ""));


## =======================================================================================================
## Over-representation Analysis, Gene Set Enrichment Analysis and Topology-based Pathway Analysis
## =======================================================================================================
tryCatch(
	expr = 
	{
		## Create sub-directories to store outputs for each data source ----------------------------------
		for(data_source in data_sources)
		{
			if(!dir.exists(paste(output_dir, "/", data_source, sep = "")))
			{
				dir.create(path = paste(output_dir, "/", data_source, sep = ""), recursive = TRUE, mode = "0777");
			}
		}

		## Prepare MsigDB gene sets -------------------------
		declare("mSigDb_Hallmarks_df", Data.frame());
		declare("mSigDb_Hallmarks_list", typed::List());
		declare("mSigDb_OncogenicC6_df", Data.frame());
		declare("mSigDb_OncogenicC6_list", typed::List());
		declare("mSigDb_ImmunologicC7_df", Data.frame());
		declare("mSigDb_ImmunologicC7_list", typed::List());

		# MSigDB Hallmark gene sets
		mSigDb_Hallmarks_df <- as.data.frame(msigdbr::msigdbr(species = species, category = "H"), stringsAsFactors = FALSE);
		debug(file_logger, "dim(mSigDb_Hallmarks_df)");
		print(dim(mSigDb_Hallmarks_df));		#8209   15
		debug(file_logger, "colnames(mSigDb_Hallmarks_df)");
		print(colnames(mSigDb_Hallmarks_df));	#"gs_cat" "gs_subcat" "gs_name" "gene_symbol" "entrez_gene" "ensembl_gene" "human_gene_symbol" "human_entrez_gene" "human_ensembl_gene" "gs_id" "gs_pmid" "gs_geoid" "gs_exact_source" "gs_url" "gs_description"
		debug(file_logger, "length(unique(mSigDb_Hallmarks_df$gs_name))");
		print(length(unique(mSigDb_Hallmarks_df$gs_name)));	#50 Hallmark gene sets

		if(gene_id_type == "ensembl_gene_id")
		{
			mSigDb_Hallmarks_list <- mSigDb_Hallmarks_df %>% base::split(x = .$ensembl_gene, f = .$gs_name);
		}else if(gene_id_type == "entrezgene_id")
		{
			mSigDb_Hallmarks_list <- mSigDb_Hallmarks_df %>% base::split(x = .$entrez_gene, f = .$gs_name);
		}else if(gene_id_type == "uniprot_gn_symbol")
		{
			mSigDb_Hallmarks_list <- mSigDb_Hallmarks_df %>% base::split(x = .$gene_symbol, f = .$gs_name);
		}

		# MSigDB Oncogenic signatures
		mSigDb_OncogenicC6_df <- as.data.frame(msigdbr::msigdbr(species = species, category = "C6"), stringsAsFactors = FALSE);
		debug(file_logger, "dim(mSigDb_OncogenicC6_df)");
		print(dim(mSigDb_OncogenicC6_df));		#34974    15
		debug(file_logger, "colnames(mSigDb_OncogenicC6_df)");
		print(colnames(mSigDb_OncogenicC6_df));
		#"gs_cat" "gs_subcat" "gs_name" "gene_symbol" "entrez_gene" "ensembl_gene" "human_gene_symbol" "human_entrez_gene" "human_ensembl_gene" "gs_id" "gs_pmid""gs_geoid" "gs_exact_source" "gs_url" "gs_description"
		debug(file_logger, "length(unique(mSigDb_OncogenicC6_df$gs_name))");
		print(length(unique(mSigDb_OncogenicC6_df$gs_name)));	#189 gene sets

		if(gene_id_type == "ensembl_gene_id")
		{
			mSigDb_OncogenicC6_list <- mSigDb_OncogenicC6_df %>% base::split(x = .$ensembl_gene, f = .$gs_name);
		}else if(gene_id_type == "entrezgene_id")
		{
			mSigDb_OncogenicC6_list <- mSigDb_OncogenicC6_df %>% base::split(x = .$entrez_gene, f = .$gs_name);
		}else if(gene_id_type == "uniprot_gn_symbol")
		{
			mSigDb_OncogenicC6_list <- mSigDb_OncogenicC6_df %>% base::split(x = .$gene_symbol, f = .$gs_name);
		}

		# MSigDB Immunologic signatures
		mSigDb_ImmunologicC7_df <- as.data.frame(msigdbr::msigdbr(species = species, category = "C7"), stringsAsFactors = FALSE);
		debug(file_logger, "dim(mSigDb_ImmunologicC7_df)");
		print(dim(mSigDb_ImmunologicC7_df));	#1118964      15
		debug(file_logger, "colnames(mSigDb_ImmunologicC7_df)");
		print(colnames(mSigDb_ImmunologicC7_df));
		#"gs_cat" "gs_subcat" "gs_name"	"gene_symbol" "entrez_gene" "ensembl_gene" "human_gene_symbol" "human_entrez_gene" "human_ensembl_gene" "gs_id" "gs_pmid" "gs_geoid" "gs_exact_source" "gs_url" "gs_description"
		debug(file_logger, "length(unique(mSigDb_ImmunologicC7_df$gs_name))");
		print(length(unique(mSigDb_ImmunologicC7_df$gs_name)));	#5219 gene sets

		if(gene_id_type == "ensembl_gene_id")
		{
			mSigDb_ImmunologicC7_list <- mSigDb_ImmunologicC7_df %>% base::split(x = .$ensembl_gene, f = .$gs_name);
		}else if(gene_id_type == "entrezgene_id")
		{
			mSigDb_ImmunologicC7_list <- mSigDb_ImmunologicC7_df %>% base::split(x = .$entrez_gene, f = .$gs_name);
		}else if(gene_id_type == "uniprot_gn_symbol")
		{
			mSigDb_ImmunologicC7_list <- mSigDb_ImmunologicC7_df %>% base::split(x = .$gene_symbol, f = .$gs_name);
		}

		info(file_logger, paste("Processing genelist with #genes: ",  length(genelist$genes), " ...", sep = ""));

		if( length(genelist$genes) < 10 )
		{
			save(list=ls(), file = paste(output_dir, "/", "FunctionalEnrichment.RData", sep = ""));
			info(file_logger, "Enrichment analysis skipped for this gene set as number of genes are less than 10. This is done to avoid too much of processing time.");
			base::sink();
			cat("Enrichment analysis skipped for this gene set as number of genes are less than 10. This is done to avoid too much of processing time.", "\n", sep = "");
			cat("The program: Functional_enrichment.R did not complete its execution because the gene set has too few genes (less than 10).", "\n", sep = "");
			q(save = "no", status = 0);
		}

		## Create a named vector of log2fc ----------------------------------------------
		declare("current_genelist_with_log2fc", Double(null_ok = TRUE), value = rep(NaN, length(genelist$genes)));
		if("log2fc" %in% names(genelist))
		{
			current_genelist_with_log2fc <- genelist$log2fc;
		}
		
		names(current_genelist_with_log2fc) <- genelist$genes;

		# order genes; required for GSEA 
		current_genelist_with_log2fc = current_genelist_with_log2fc[order(current_genelist_with_log2fc, decreasing = TRUE)];

		## A. Perform Over-representation analysis (ORA) ==================================================
		# A.1 On MSigDB Hallmark gene sets ------------------------------
		if("MSigDBHallmarks" %in% data_sources)
		{
			info(file_logger, "Performing ORA on MSigDB Hallmark genesets using clusterProfiler's enricher ..");

			enriched_MSigDB_hallmark_genesets_using_enricher_df <- runORAonMSigDB(de_genes = names(current_genelist_with_log2fc), ref_genes = ref_genelist, id_type = gene_id_type, msigdb_gs_df = mSigDb_Hallmarks_df, gene_ann_df = gene_info_df, threshold = threshold, output_file_prefix = paste(output_dir, "/MSigDBHallmarks/", output_prefix, "enriched_MSigDB_hallmarks_genesets_using_ORA", sep = ""));

			info(file_logger, paste("#MSigDB hallmark signatures with p.adjust < ", threshold, ": ", nrow(enriched_MSigDB_hallmark_genesets_using_enricher_df), sep = ""));
			if(nrow(enriched_MSigDB_hallmark_genesets_using_enricher_df))
			{
				write.table(x = enriched_MSigDB_hallmark_genesets_using_enricher_df, file = paste(output_dir, "/MSigDBHallmarks/", output_prefix, "enriched_MSigDB_hallmarks_genesets_using_ORA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
			}
		}

		# A.2 On MSigDB oncogenic signatures ----------------------------
		if("MSigDBOncogenic" %in% data_sources)
		{
			info(file_logger, "Performing ORA on MSigDB oncogenic signatures using clusterProfiler's enricher ..");

			enriched_MSigDB_oncogenicC6_genesets_using_enricher_df <- runORAonMSigDB(de_genes = names(current_genelist_with_log2fc), ref_genes = ref_genelist, id_type = gene_id_type, msigdb_gs_df = mSigDb_OncogenicC6_df, gene_ann_df = gene_info_df, threshold = threshold, output_file_prefix = paste(output_dir, "/MSigDBOncogenic/", output_prefix, "enriched_MSigDB_oncogenicC6_genesets_using_ORA", sep = ""));
			
			info(file_logger, paste("#MSigDB oncogenic signatures with p.adjust < ", threshold, ": ", nrow(enriched_MSigDB_oncogenicC6_genesets_using_enricher_df), sep = ""));

			if(nrow(enriched_MSigDB_oncogenicC6_genesets_using_enricher_df))
			{
				write.table(x = enriched_MSigDB_oncogenicC6_genesets_using_enricher_df, file = paste(output_dir, "/MSigDBOncogenic/", output_prefix, "enriched_MSigDB_oncogenicC6_genesets_using_ORA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
			}
		}
			
			# A.3 On MSigDB immunologic signatures --------------------------
			if("MSigDBImmunologic" %in% data_sources)
			{
				info(file_logger, "Performing ORA on MSigDB immunologic signatures using clusterProfiler's enricher ..");

				enriched_MSigDB_ImmunologicC7_genesets_using_enricher_df <- runORAonMSigDB(de_genes = names(current_genelist_with_log2fc), ref_genes = ref_genelist, id_type = gene_id_type, msigdb_gs_df = mSigDb_ImmunologicC7_df, gene_ann_df = gene_info_df, threshold = threshold, output_file_prefix = paste(output_dir, "/MSigDBImmunologic/", output_prefix, "enriched_MSigDB_immunologicC7_genesets_using_ORA", sep = ""));

				info(file_logger, paste("#MSigDB immunologic signatures with p.adjust < ", threshold, ": ", nrow(enriched_MSigDB_ImmunologicC7_genesets_using_enricher_df), sep = ""));

				if(nrow(enriched_MSigDB_ImmunologicC7_genesets_using_enricher_df))
				{
					write.table(x = enriched_MSigDB_ImmunologicC7_genesets_using_enricher_df, file = paste(output_dir, "/MSigDBImmunologic/", output_prefix, "enriched_MSigDB_immunologicC7_genesets_using_ORA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
				}
			}

			# A.4 On GO (Biological Process) --------------------------------
			if("GOBP" %in% data_sources)
			{
				info(file_logger, "Performing ORA on GO (Biological Processs) using clusterProfiler's enrichGO ..");

				enriched_GOBP_genesets_using_enrichGO_df <- runORAonGO(de_genes = names(current_genelist_with_log2fc), ref_genes = ref_genelist, id_type = gene_id_type, gene_ann_df = gene_info_df, species = species, ontology_type = "BP", threshold = threshold, output_file_prefix = paste(output_dir, "/GOBP/", output_prefix, "enriched_GOBP_genesets_using_ORA", sep = ""));

				info(file_logger, paste("#GOBP terms with p.adjust < ", threshold, ": ", nrow(enriched_GOBP_genesets_using_enrichGO_df), sep = ""));

				if(nrow(enriched_GOBP_genesets_using_enrichGO_df))
				{
					write.table(x = enriched_GOBP_genesets_using_enrichGO_df, file = paste(output_dir, "/GOBP/", output_prefix, "enriched_GOBP_genesets_using_ORA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
				}
			}

			# A.5 On GO (Molcular Function) ---------------------------------
			if("GOMF" %in% data_sources)
			{
				info(file_logger, "Performing ORA on GO (Molecular Function) using clusterProfiler's enrichGO ..");
				
				enriched_GOMF_genesets_using_enrichGO_df <- runORAonGO(de_genes = names(current_genelist_with_log2fc), ref_genes = ref_genelist, id_type = gene_id_type, gene_ann_df = gene_info_df, species = species, ontology_type = "MF", threshold = threshold, output_file_prefix = paste(output_dir, "/GOMF/", output_prefix, "enriched_GOMF_genesets_using_ORA", sep = ""));

				info(file_logger, paste("#GOMF terms with p.adjust < ", threshold, ": ", nrow(enriched_GOMF_genesets_using_enrichGO_df), sep = ""));

				if(nrow(enriched_GOMF_genesets_using_enrichGO_df))
				{
					write.table(x = enriched_GOMF_genesets_using_enrichGO_df, file = paste(output_dir, "/GOMF/", output_prefix, "enriched_GOMF_genesets_using_ORA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
				}
			}

			# A.6 On GO (Cellular Component) --------------------------------
			if("GOCC" %in% data_sources)
			{
				info(file_logger, "Performing ORA on GO (Cellular Component) using clusterProfiler's enrichGO ..");

				enriched_GOCC_genesets_using_enrichGO_df <- runORAonGO(de_genes = names(current_genelist_with_log2fc), ref_genes = ref_genelist, id_type = gene_id_type, gene_ann_df = gene_info_df, species = species, ontology_type = "CC", threshold = threshold, output_file_prefix = paste(output_dir, "/GOCC/", output_prefix, "enriched_GOCC_genesets_using_ORA", sep = ""));
				
				info(file_logger, paste("#GOCC terms with p.adjust < ", threshold, ": ", nrow(enriched_GOCC_genesets_using_enrichGO_df), sep = ""));

				if(nrow(enriched_GOCC_genesets_using_enrichGO_df))
				{
					write.table(x = enriched_GOCC_genesets_using_enrichGO_df, file = paste(output_dir, "/GOCC/", output_prefix, "enriched_GOCC_genesets_using_ORA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
				}
			}

			if(species == "Homo sapiens")
			{
				## DO, NCG and DGN enrichment is applicable only for 'Homo sapiens' and can work only with Entrez gene ids.

			# A.7 On Disease Ontology (DO) terms ----------------------------
			if("DO" %in% data_sources)
			{
				info(file_logger, "Performing ORA on DO genesets using DOSE's enrichDO ..");
				
				enriched_DO_genesets_using_enrichDO_df <- runORAonDO(de_genes = names(current_genelist_with_log2fc), ref_genes = ref_genelist, id_type = gene_id_type, mapped_entrez_gene_ids, gene_ann_df = gene_info_df, threshold = threshold, output_file_prefix = paste(output_dir, "/DO/", output_prefix, "enriched_DO_genesets_using_ORA", sep = ""));

				info(file_logger, paste("#DO genesets with p.adjust < ", threshold, ": ", nrow(enriched_DO_genesets_using_enrichDO_df), sep = ""));

				if(nrow(enriched_DO_genesets_using_enrichDO_df))
				{
					write.table(x = enriched_DO_genesets_using_enrichDO_df, file = paste(output_dir, "/DO/", output_prefix, "enriched_DO_genesets_using_ORA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
				}
			}

			# A.8 On Network of Cancer Genes (NCG) database -----------------
			if("NCG" %in% data_sources)
			{
				info(file_logger, "Performing ORA on NCG genesets using DOSE's enrichNCG ..");

				enriched_NCG_genesets_using_enrichNCG_df <- runORAonNCG(de_genes = names(current_genelist_with_log2fc), ref_genes = ref_genelist, id_type = gene_id_type, mapped_entrez_gene_ids = mapped_entrez_gene_ids, gene_ann_df = gene_info_df, threshold = threshold, output_file_prefix = paste(output_dir, "/NCG/", output_prefix, "enriched_NCG_genesets_using_ORA", sep = ""));
				
				info(file_logger, paste("#NCG genesets with p.adjust < ", threshold, ": ", nrow(enriched_NCG_genesets_using_enrichNCG_df), sep = ""));

				if(nrow(enriched_NCG_genesets_using_enrichNCG_df))
				{
					write.table(x = enriched_NCG_genesets_using_enrichNCG_df, file = paste(output_dir, "/NCG/", output_prefix, "enriched_NCG_genesets_using_ORA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
				}
			}

			# A.9 On DisGeNET database --------------------------------------
			if("DisGeNET" %in% data_sources)
			{
				info(file_logger, "Performing ORA on DisGeNET genesets using DOSE's enrichDGN ..");
				
				enriched_DGN_genesets_using_enrichDGN_df <- runORAonDGN(de_genes = names(current_genelist_with_log2fc), ref_genes = ref_genelist, id_type = gene_id_type, mapped_entrez_gene_ids = mapped_entrez_gene_ids, gene_ann_df = gene_info_df, threshold = threshold, output_file_prefix = paste(output_dir, "/DisGeNET/", output_prefix, "enriched_DisGeNET_genesets_using_ORA", sep = ""));

				info(file_logger, paste("#DisGeNET genesets with p.adjust < ", threshold, ": ", nrow(enriched_DGN_genesets_using_enrichDGN_df), sep = ""));

				if(nrow(enriched_DGN_genesets_using_enrichDGN_df))
				{
					write.table(x = enriched_DGN_genesets_using_enrichDGN_df, file = paste(output_dir, "/DisGeNET/", output_prefix, "enriched_DisGeNET_genesets_using_ORA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
				}
			}
			}else
			{
				#debug(file_logger, "ORA on Disease Ontology (DO), Network of Cancer Genes (NCG) and DisGeNET are only valid for Homo sapiens, and hence will be ignored in the current execution.");
				debug(file_logger, "ORA on Disease Ontology (DO) and Network of Cancer Genes (NCG) are only valid for Homo sapiens, and hence will be ignored in the current execution.");
			}


			# A.10 On KEGG pathways -----------------------------------------
			if("KEGG" %in% data_sources)
			{
				info(file_logger, "Performing ORA on KEGG pathways using clusterProfiler's enrichKEGG ..");
				
				enriched_kegg_pathways_using_enrichKEGG_df <- runORAonKEGG(de_genes = names(current_genelist_with_log2fc), ref_genes = ref_genelist, id_type = gene_id_type, mapped_entrez_gene_ids = mapped_entrez_gene_ids, gene_ann_df = gene_info_df, species = species, threshold = threshold, output_file_prefix = paste(output_dir, "/KEGG/", output_prefix, "enriched_KEGG_pathways_using_ORA", sep = ""));
				
				info(file_logger, paste("#KEGG pathways with p.adjust < ", threshold, ": ", nrow(enriched_kegg_pathways_using_enrichKEGG_df), sep = ""));
				
				if(nrow(enriched_kegg_pathways_using_enrichKEGG_df))
				{
					write.table(x = enriched_kegg_pathways_using_enrichKEGG_df, file = paste(output_dir, "/KEGG/", output_prefix, "enriched_KEGG_pathways_using_ORA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
				}
			}

			if(species %in% supported_species_for_ReactomePA)
			{
			# A.11 On Reactome pathways -------------------------------------
			if("Reactome" %in% data_sources)
			{
				info(file_logger, "Performing ORA on Reactome pathways using ReactomePA ..");

				enriched_reactome_pathways_using_enrichPathway_df <- runORAonReactome(de_genes = names(current_genelist_with_log2fc), ref_genes = ref_genelist, id_type = gene_id_type, mapped_entrez_gene_ids = mapped_entrez_gene_ids, gene_ann_df = gene_info_df, species = species, threshold = threshold, output_file_prefix = paste(output_dir, "/Reactome/", output_prefix, "enriched_Reactome_pathways_using_ORA", sep = ""));

				info(file_logger, paste("#Reactome pathways with p.adjust < ", threshold, ": ", nrow(enriched_reactome_pathways_using_enrichPathway_df), sep = ""));

				if(nrow(enriched_reactome_pathways_using_enrichPathway_df))
				{
					write.table(x = enriched_reactome_pathways_using_enrichPathway_df, file = paste(output_dir, "/Reactome/", output_prefix, "enriched_Reactome_pathways_using_ORA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
				}
			}
			}else
			{
				debug(file_logger, paste("ORA on Reactome pathways is only supported for species: ", paste(supported_species_for_ReactomePA, collapse = ","), ". Hence will be ignored for current execution.", sep = ""));
			}

			save(list=ls(), file = paste(output_dir, "/", "FunctionalEnrichment.RData", sep = ""));

			## B. Perfrom Gene set enrichment analsysis (GSEA) ================================================
			if("log2fc" %in% names(genelist))
			{
				# B.1 On MSigDB Hallmark gene sets ------------------------------
				if("MSigDBHallmarks" %in% data_sources)
				{
					info(file_logger, "Performing GSEA on MSigDB Hallmark genesets using fgsea ..");
					
					enriched_MSigDB_hallmark_genesets_using_fgsea_df <- runGSEAonMSigDB(de_genes = current_genelist_with_log2fc, geneset_list = mSigDb_Hallmarks_list, id_type = gene_id_type, gene_ann_df = gene_info_df, threshold = threshold, output_file_prefix = paste(output_dir, "/MSigDBHallmarks/", output_prefix, "enriched_MSigDB_hallmarks_genesets_using_GSEA", sep = ""));

					info(file_logger, paste("#MSigDB hallmark signatures with p.adjust < ", threshold, ": ", nrow(enriched_MSigDB_hallmark_genesets_using_fgsea_df), sep = ""));

					if(nrow(enriched_MSigDB_hallmark_genesets_using_fgsea_df))
					{
						write.table(x = enriched_MSigDB_hallmark_genesets_using_fgsea_df, file = paste(output_dir, "/MSigDBHallmarks/", output_prefix, "enriched_MSigDB_hallmarks_genesets_using_GSEA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
					}
				}

				# B.2 On MSigDB oncogenic signatures ----------------------------
				if("MSigDBOncogenic" %in% data_sources)
				{
					info(file_logger, "Performing GSEA on MSigDB oncogenic signatures using fgsea ..");

					enriched_MSigDB_oncogenicC6_genesets_using_fgsea_df <- runGSEAonMSigDB(de_genes = current_genelist_with_log2fc, geneset_list = mSigDb_OncogenicC6_list, id_type = gene_id_type, gene_ann_df = gene_info_df, threshold = threshold, output_file_prefix = paste(output_dir, "/MSigDBOncogenic/", output_prefix, "enriched_MSigDB_oncogenicC6_genesets_using_GSEA", sep = ""));

					info(file_logger, paste("#MSigDB oncogenic signatures with p.adjust < ", threshold, ": ", nrow(enriched_MSigDB_oncogenicC6_genesets_using_fgsea_df), sep = ""));

					if(nrow(enriched_MSigDB_oncogenicC6_genesets_using_fgsea_df))
					{
						write.table(x = enriched_MSigDB_oncogenicC6_genesets_using_fgsea_df, file = paste(output_dir, "/MSigDBOncogenic/", output_prefix, "enriched_MSigDB_oncogenicC6_genesets_using_GSEA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
					}
				}

				# B.3 On MSigDB immunologic signatures --------------------------
				if("MSigDBImmunologic" %in% data_sources)
				{
					info(file_logger, "Performing GSEA on MSigDB immunologic signatures using fgsea ..");
					
					enriched_MSigDB_immunologicC7_genesets_using_fgsea_df <- runGSEAonMSigDB(de_genes = current_genelist_with_log2fc, geneset_list = mSigDb_ImmunologicC7_list, id_type = gene_id_type, gene_ann_df = gene_info_df, threshold = threshold, output_file_prefix = paste(output_dir, "/MSigDBImmunologic/", output_prefix, "enriched_MSigDB_immunologicC7_genesets_using_GSEA", sep = ""));

					info(file_logger, paste("#MSigDB immunologic signatures with p.adjust < ", threshold, ": ", nrow(enriched_MSigDB_immunologicC7_genesets_using_fgsea_df), sep = ""));

					if(nrow(enriched_MSigDB_immunologicC7_genesets_using_fgsea_df))
					{
						write.table(x = enriched_MSigDB_immunologicC7_genesets_using_fgsea_df, file = paste(output_dir, "/MSigDBImmunologic/", output_prefix, "enriched_MSigDB_immunologicC7_genesets_using_GSEA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
					}
				}
				
				# B.4 On GO (Biological Process) --------------------------------
				if("GOBP" %in% data_sources)
				{
					info(file_logger, "Performing GSEA on GO (Biological Processs) using clusterProfiler's gseGO ..");
					
					enriched_GOBP_genesets_using_gseGO_df <- runGSEAonGO(de_genes = current_genelist_with_log2fc, id_type = gene_id_type, gene_ann_df = gene_info_df, species = species, ontology_type = "BP", threshold = threshold, output_file_prefix = paste(output_dir, "/GOBP/", output_prefix, "enriched_GOBP_genesets_using_GSEA", sep = ""));
					
					info(file_logger, paste("#GOBP terms with p.adjust < ", threshold, ": ", nrow(enriched_GOBP_genesets_using_gseGO_df), sep = ""));

					if(nrow(enriched_GOBP_genesets_using_gseGO_df))
					{
						write.table(x = enriched_GOBP_genesets_using_gseGO_df, file = paste(output_dir, "/GOBP/", output_prefix, "enriched_GOBP_genesets_using_GSEA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
					}
				}

				# B.5 On GO (Molcular Function) ---------------------------------
				if("GOMF" %in% data_sources)
				{
					info(file_logger, "Performing GSEA on GO (Moleculr Function) using clusterProfiler's gseGO ..");
					
					enriched_GOMF_genesets_using_gseGO_df <- runGSEAonGO(de_genes = current_genelist_with_log2fc, id_type = gene_id_type, gene_ann_df = gene_info_df, species = species, ontology_type = "MF", threshold = threshold, output_file_prefix = paste(output_dir, "/GOMF/", output_prefix, "enriched_GOMF_genesets_using_GSEA", sep = ""));
					
					info(file_logger, paste("#GOMF terms with p.adjust < ", threshold, ": ", nrow(enriched_GOMF_genesets_using_gseGO_df), sep = ""));

					if(nrow(enriched_GOMF_genesets_using_gseGO_df))
					{
						write.table(x = enriched_GOMF_genesets_using_gseGO_df, file = paste(output_dir, "/GOMF/", output_prefix, "enriched_GOMF_genesets_using_GSEA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
					}
				}

				# B.6 On GO (Cellular Component) --------------------------------
				if("GOCC" %in% data_sources)
				{
					info(file_logger, "Performing GSEA on GO (Cellular Component) using clusterProfiler's gseGO ..");
					
					enriched_GOCC_genesets_using_gseGO_df <- runGSEAonGO(de_genes = current_genelist_with_log2fc, id_type = gene_id_type, gene_ann_df = gene_info_df, species = species, ontology_type = "CC", threshold = threshold, output_file_prefix = paste(output_dir, "/GOCC/", output_prefix, "enriched_GOCC_genesets_using_GSEA", sep = ""));
					
					info(file_logger, paste("#GOCC terms with p.adjust < ", threshold, ": ", nrow(enriched_GOCC_genesets_using_gseGO_df), sep = ""));

					if(nrow(enriched_GOCC_genesets_using_gseGO_df))
					{
						write.table(x = enriched_GOCC_genesets_using_gseGO_df, file = paste(output_dir, "/GOCC/", output_prefix, "enriched_GOCC_genesets_using_GSEA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
					}
				}

				if(species == "Homo sapiens")
				{
				# B.7 On Disease Ontology (DO) terms ----------------------------
				if("DO" %in% data_sources)
				{
					info(file_logger, "Performing GSEA on Disease Ontology (DO) genesets using DOSE's gseDO ..");

					enriched_DO_genesets_using_gseDO_df <- runGSEAonDO(de_genes = current_genelist_with_log2fc, id_type = gene_id_type, gene_ann_df = gene_info_df, mapped_entrez_gene_ids, threshold = threshold, output_file_prefix = paste(output_dir, "/DO/", output_prefix, "enriched_DO_genesets_using_GSEA", sep = ""));

					info(file_logger, paste("#DO genesets with p.adjust < ", threshold, ": ", nrow(enriched_DO_genesets_using_gseDO_df), sep = ""));

					if(nrow(enriched_DO_genesets_using_gseDO_df))
					{
						write.table(x = enriched_DO_genesets_using_gseDO_df, file = paste(output_dir, "/DO/", output_prefix, "enriched_DO_genesets_using_GSEA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
					}
				}

				# B.8 On Network of Cancer Genes (NCG) database -----------------
				if("NCG" %in% data_sources)
				{
					info(file_logger, "Performing GSEA on Network of Cancer Genes (NCG) genesets using DOSE's gseNCG ..");
					
					enriched_NCG_genesets_using_gseNCG_df <- runGSEAonNCG(de_genes = current_genelist_with_log2fc, id_type = gene_id_type, mapped_entrez_gene_ids, gene_ann_df = gene_info_df, threshold = threshold, output_file_prefix = paste(output_dir, "/NCG/", output_prefix, "enriched_NCG_genesets_using_GSEA", sep = ""));

					info(file_logger, paste("#NCG genesets with p.adjust < ", threshold, ": ", nrow(enriched_NCG_genesets_using_gseNCG_df), sep = ""));

					if(nrow(enriched_NCG_genesets_using_gseNCG_df))
					{
						write.table(x = enriched_NCG_genesets_using_gseNCG_df, file = paste(output_dir, "/NCG/", output_prefix, "enriched_NCG_genesets_using_GSEA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
					}
				}

				# B.9 On DisGeNET database --------------------------------------
				if("DisGeNET" %in% data_sources)
				{
					info(file_logger, "Performing GSEA on DisGeNET genesets using DOSE's gseDGN ..");
					
					enriched_DGN_genesets_using_gseDGN_df <- runGSEAonDGN(de_genes = current_genelist_with_log2fc, id_type = gene_id_type, mapped_entrez_gene_ids, gene_ann_df = gene_info_df, threshold = threshold, output_file_prefix = paste(output_dir, "/DisGeNET/", output_prefix, "enriched_DisGeNET_genesets_using_GSEA", sep = ""));
					
					info(file_logger, paste("#DisGeNET genesets with p.adjust < ", threshold, ": ", nrow(enriched_DGN_genesets_using_gseDGN_df), sep = ""));

					if(nrow(enriched_DGN_genesets_using_gseDGN_df))
					{
						write.table(x = enriched_DGN_genesets_using_gseDGN_df, file = paste(output_dir, "/DisGeNET/", output_prefix, "enriched_DisGeNET_genesets_using_GSEA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
					}
				}
				}else
				{
					debug(file_logger, "GSEA on Disease Ontology (DO), Network of Cancer Genes (NCG) and DisGeNET is only valid for Homo sapiens, and hence will be ignored in the current execution.");
				}

				# B.10 On KEGG pathways -----------------------------------------
				if("KEGG" %in% data_sources)
				{
					info(file_logger, "Performing GSEA on KEGG pathways using clusterProfiler's gseKEGG ..");
					
					enriched_kegg_pathways_using_gseKEGG_df <- runGSEAonKEGG(de_genes = current_genelist_with_log2fc, id_type = gene_id_type, mapped_entrez_gene_ids, gene_ann_df = gene_info_df, species = species, threshold = threshold, output_file_prefix = paste(output_dir, "/KEGG/", output_prefix, "enriched_KEGG_pathways_using_GSEA", sep = ""));
					
					info(file_logger, paste("#KEGG pathways with p.adjust < ", threshold, ": ", nrow(enriched_kegg_pathways_using_gseKEGG_df), sep = ""));
					
					if(nrow(enriched_kegg_pathways_using_gseKEGG_df))
					{
						write.table(x = enriched_kegg_pathways_using_gseKEGG_df, file = paste(output_dir, "/KEGG/", output_prefix, "enriched_KEGG_pathways_using_GSEA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
					}
				}

				if(species %in% supported_species_for_ReactomePA)
				{
				# B.11 On Reactome pathways -----------------------------------
				if("Reactome" %in% data_sources)
				{
					info(file_logger, "Performing GSEA on Reactome pathways using ReactomePA's gsePathway ..");

					enriched_reactome_pathways_using_gsePathway_df <- runGSEAonReactome(de_genes = current_genelist_with_log2fc, id_type = gene_id_type, mapped_entrez_gene_ids = mapped_entrez_gene_ids, gene_ann_df = gene_info_df, species = species, threshold = threshold, output_file_prefix = paste(output_dir, "/Reactome/", output_prefix, "enriched_Reactome_pathways_using_GSEA", sep = ""));

					info(file_logger, paste("#Reactome pathways with p.adjust < ", threshold, ": ", nrow(enriched_reactome_pathways_using_gsePathway_df), sep = ""));

					if(nrow(enriched_reactome_pathways_using_gsePathway_df))
					{
						write.table(x = enriched_reactome_pathways_using_gsePathway_df, file = paste(output_dir, "/Reactome/", output_prefix, "enriched_Reactome_pathways_using_GSEA.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
					}
				}
				}else
				{
					info(file_logger, paste("ORA on Reactome pathways is only supported for species: ", paste(supported_species_for_ReactomePA, collapse = ","), ". Hence will be ignored for current execution.", sep = ""));
				}

			}else
			{
				info(file_logger, "Skipped GSEA analysis because input list of genes do not have associated log2-fold-change data.");
			}
			save(list=ls(), file = paste(output_dir, "/", "FunctionalEnrichment.RData", sep = ""));


			# C. Perform Topology-based pathway enrichment analysis ===========================================
			if(species %in% c("Homo sapiens", "Mus musculus"))
			{
			if(run_SPIA && !is.null(ref_genelist) && ("log2fc" %in% names(genelist)))
			{
				info(file_logger, "Currently we are using an outdated KEGG data (downloaded on 09/07/2012 and bundled with R package SPIA) for pathway enrichment using SPIA for only Homo sapiens and Mus musculus. There is a provision to create appropriate Rdata on latest KEGG data using SPIA::makeSPIAdata function, but currently that is out of scope.");
				
				# C.1 On KEGG pathways -------------------------------------------
				if(!dir.exists(paste(output_dir, "/", "KEGG", sep = "")))
				{
					dir.create(path = paste(output_dir, "/", "KEGG", sep = ""), recursive = TRUE, mode = "0777");
				}

				info(file_logger, "Performing topology-based KEGG pathway enrichment analysis using SPIA ...");
				enriched_kegg_pathways_using_spia_df <- runSPIA(de_genes = current_genelist_with_log2fc, ref_genes = ref_genelist, id_type = gene_id_type, mapped_entrez_gene_ids, species = species, gene_ann_df = gene_info_df, threshold = threshold, output_file_prefix = paste(output_dir, "/KEGG/", output_prefix, "enriched_KEGG_pathways_using_topology", sep = ""));
				
				info(file_logger, paste("\n", "#KEGG pathways with pGFdr < ", threshold, ": ", nrow(enriched_kegg_pathways_using_spia_df), sep = ""));
				
				if(nrow(enriched_kegg_pathways_using_spia_df))
				{
					write.table(x = enriched_kegg_pathways_using_spia_df, file = paste(output_dir, "/KEGG/", output_prefix, "enriched_KEGG_pathways_using_topology.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
				}
				
			}else if(run_SPIA && (is.null(ref_genelist) || !("log2fc" %in% names(genelist))))
			{
				info(file_logger, "Skipped topology-based KEGG pathway enrichment analysis using SPIA because either reference gene list is not available or input list of genes do not have associated log2-fold-change data.");
			}
			}else
			{
				info(file_logger, "Currently there is support for only Homo sapiens and Mus musculus for SPIA analysis. Therefore, this analysis will be skipped for current execution.");
			}

			save(list=ls(), file = paste(output_dir, "/", "FunctionalEnrichment.RData", sep = ""));
	},
	warn = function(w)
	{
		warn(file_logger, "Caught a warning while performing enrichment analysis. Something unexpected must have happened.");
		print(w);
		base::sink();
		print(w);
		cat("The program: Functional_enrichment.R failed at 'Enrichment analysis'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	},
	error = function(e)
	{
		error(file_logger, "An error occured while performing enrichment analysis.");
		print(e);
		base::sink();
		print(e);
		cat("The program: Functional_enrichment.R failed at 'Enrichment analysis'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	}
)
save(list=ls(), file = paste(output_dir, "/", "FunctionalEnrichment.RData", sep = ""));


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
cat("The program: Functional_enrichment.R executed successfully!", "\n", sep = "");

