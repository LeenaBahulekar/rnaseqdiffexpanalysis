#!/usr/bin/Rscript --vanilla

#####################################################################################
###----------------------------------------------------------------------------------
###   Program: Calculate Topology and DrugTarget proximity of the query genes in PPI network
###   Developed by: Dr Arindam Deb         
###----------------------------------------------------------------------------------
###==================================================================================
### User needs to pass the following parameters as arguments -
###       Input 1. Query gene data (must have a column named as 'gene_name'), (.csv file).
###       Input 2. Protein protein interaction data (two columns - Protein_A,Protein_B), (.csv file).
###       Input 3. Drug target data (two columns - DrugTarget,DrugName. (tab separated .txt file).
###                Each row of 'DrugName' consists list of Drug names separated by '|').
###       Input 4. Column name that contains the HGNC gene symbols.
###       Input 5. Optional - Column name that contains the log2 Fold change values.
###                It is used to distinguishably colour (Up/Down) the nodes in the network.
###       Input 6. Output directory Path and Directory name (e.g., - D:\...\...\...\Output)
### ** Note: All the gene names (query_gene/DrugTarget) should be in HGNC symbol.
### ** Note: The proteins in PPI are transformed into their respective HGNC gene symbols (for mapping purpose).
###
###
### Program generates TWO output data files anf 2 Interactive Image files -
###       Output 1. Gene_Topology_attribute_in_PPI (.csv file)
###       Output 2. QueryGene_DrugTarget_ShortestPath_matrix (.csv file)
###       Image  1. SubNet_PPI_QG_DT.html
###       Image  2. Heatmap_QG_DT_ShortestPath.html
###
### 
### ** Note: Output1 is used as attribute file for network building.
###
###
###==================================================================================
###----------------------------------------------------------------------------------
#####################################################################################

# setwd("D:\\MultiOmics\\RNAseq\\Gene_Drug_Association\\")

library("optparse");	      # command line option parser
library("data.table");	    # Read input data fast
library("log4r");           # Logging library
library("biomaRt");         # Provides annotation for genes
library("igraph")
library("dplyr")
library("tidyr")
library("CINNA")
library("visNetwork")
library("heatmaply")

###--------------------------------------------------------------------------###
###                                   INPUT
###--------------------------------------------------------------------------###
args_list <- list(
  optparse::make_option(opt_str = c("--genelist_file"), type = "character", default = NULL, help = "List of query genes in CSV or TSV format. Default: NULL"),
  optparse::make_option(opt_str = c("--ppi_file"), type = "character", default = NULL, help = "PPI file in CSV or TSV format. Default: NULL"),
  optparse::make_option(opt_str = c("--drug_target_file"), type = "character", default = NULL, help = "Drug-Target mapping file in TSV format. Default: NULL"),
  optparse::make_option(opt_str = c("--gene_id_column"), type = "integer", default = 1, help = "Column number in --genelist_file which contains gene identifiers. Default: 1"),
  optparse::make_option(opt_str = c("--log2fc_column"), type = "character", default = NULL, help = "Column name in --genelist_file which contains log2-fold-change values. This is typically available if the gene list is an output of some differential expression analysis. Default: NULL"),
  optparse::make_option(opt_str = c("--output_dir"), type = "character", default = "./", help = "Output directory path. Default: ./"),
  optparse::make_option(opt_str = c("--log_file"), type = "character", default = "Drug_Gene_association.log", help = "Log file. Default: Drug_Gene_association.log")
)
arguments <- optparse::parse_args(optparse::OptionParser(option_list = args_list));


###--------------------------------------------------------------------------###
##                             Set imp. variables
###--------------------------------------------------------------------------###
genelist_file <- arguments$genelist_file;
ppi_file <- arguments$ppi_file;
drug_target_file <- arguments$drug_target_file;
gene_id_column <- arguments$gene_id_column;

if(!is.null(arguments$log2fc_column) && (arguments$log2fc_column == 'NULL'))
{
	arguments$log2fc_column <- NULL;
}
log2fc_column <- arguments$log2fc_column;
output_dir <- arguments$output_dir;
log_file <- arguments$log_file;

if(!dir.exists(output_dir))
{
	dir.create(output_dir,recursive = TRUE, mode = "0777")
}

log_dir <- paste(output_dir, "/", "log", sep = "");
if(!dir.exists(log_dir))
{
	dir.create(path = log_dir, recursive = TRUE, mode = "0777");
}

log_file <- paste(log_dir, "/", log_file, sep = "");
file_logger <- log4r::logger(threshold = "DEBUG", appenders = log4r::file_appender(file = log_file, append = TRUE));
info(file_logger, "Processing starts");
# Additionally, for capturing stdouts 
base::sink(file = log_file, append = TRUE);

info(file_logger, paste("Input parameters are:\n", paste(head(x = names(arguments), n = -1), head(x = as.character(arguments), n = -1), sep = ":", collapse = ", "), sep = ""));


# =========================================================
## Custom functions
## ========================================================

## Connect to bioMart for gene annotation -----------------------
annotateWithBioMart <- function(biomart_database_name, species, gene_id_type, gene_ids, file_logger)
{
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
									base::sink(); cat("The program: gene_drug_association.R failed because it could not connect to bioMart for gene annotation.", "\n", sep = "");
									q(save = "no", status = 1);
								}
							)
						}
					)
				}
			)
		}
	)

	data_retrieved <-  FALSE;
	reattempt <- as.integer(0)
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
					base::sink(); cat("The program: gene_drug_association.R failed while retrieving data from bioMart.", "\n", sep = "");
					q(save = "no", status = 1);
				}
			}
		)   
	}

	return(gene_info_df);
}



## Identify type of gene identifier - Ensembl or Entrez or Uniprot Gene Symbol ---------
identifyGeneIdType <- function(ids)
{
	id_type <- NULL;

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


###--------------------------------------------------------------------------###
### Read input list of genes
###--------------------------------------------------------------------------###

tryCatch(
  expr = 
  {
    # query_gene <- read.table(".\\Input\\30UpDown_Carcinoma_vs_Adenoma_logFC1_adjpval0.05.csv", sep = ",", header = T)
    query_gene_df <- as.data.frame(data.table::fread(input = genelist_file, sep = "auto", header = TRUE, check.names = FALSE), stringsAsFactors = FALSE);
    ppi_df <- as.data.frame(data.table::fread(input = ppi_file, sep = "auto", header = TRUE, check.names = FALSE), stringsAsFactors = FALSE);
    drug_target_df <- as.data.frame(data.table::fread(input = drug_target_file, sep = "\t", header = TRUE, check.names = FALSE), stringsAsFactors = FALSE);
  },
  warning = function(w)
  {
    warn(file_logger, "Following warning encountered while reading input files.");
    print(w);
    base::sink(); cat("The program: gene_drug_association.R failed while reading input files. Check log for details: ", log_file, "\n", sep = "");
    q(save = "no", status = 1);
  },
  error = function(e)
  {
    error(file_logger, "Following error occured while reading input files: ");
    print(e);
    base::sink(); cat("The program: gene_drug_association.R failed while reading input files. Check log for details: ", log_file, "\n", sep = "");
    q(save = "no", status = 1);
  }
)

save(list = ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));


if(!gene_id_column %in% 1:ncol(query_gene_df))
{
  save(list = ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));
  error(file_logger, "Specified column for gene list is NOT found in --genelist_file");
  base::sink(); cat("The program: gene_drug_association.R failed because specified column for gene list via --gene_id_column is NOT found in --genelist_file ", "\n", sep = "");
  q(save = "no", status = 1);
}

## Get list of query genes ------------------------------------------
query_genelist <- as.character(query_gene_df[,gene_id_column]);


###--------------------------------------------------------------------------###
### Check gene identifier type, and perform required conversion, if any
###--------------------------------------------------------------------------###
info(file_logger, "Checking gene identifier types in all input files, and performing suitable conversion (if required) ...");

tryCatch(
  expr = 
  {
    ## Determine type of gene identifiers in input gene list, PPI and Drug-Target -----------------------
    gene_id_type <- identifyGeneIdType(ids = query_genelist);
    info(file_logger, paste("Inferred gene identifier type in the input list of genes is: ", gene_id_type, "\n", sep = ""));

    ppi_gene_id_type <- identifyGeneIdType(ids = unique(as.character(ppi_df[,'Protein_A'])));
    info(file_logger, paste("Inferred gene identifier type in the PPI data is: ", ppi_gene_id_type, "\n", sep = ""));

    drug_target_gene_id_type <- identifyGeneIdType(ids = unique(as.character(drug_target_df$DrugTarget)));
    info(file_logger, paste("Inferred gene identifier type in the Drug-Target data is: ", drug_target_gene_id_type, "\n", sep = ""));

    if(any(c(gene_id_type, ppi_gene_id_type, drug_target_gene_id_type) == "Unknown"))
		{
      declare("temp_id_types", Character(), value = c(gene_id_type, ppi_gene_id_type, drug_target_gene_id_type));
      names(temp_id_types) <- c("genes_of_interest", "ppi_data", "drug_target_data");

      declare("unknown_ids", Character());
      unknown_ids <- names(temp_id_types)[which(c(temp_id_types) == "Unknown")];
      unknown_ids <- paste(unknown_ids, collapse = ",");

			save(list=ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));
			error(file_logger, paste("Unable to identify gene id type of ", unknown_ids, ". Make sure that all ids are one of ensembl_gene_id (e.g., ENSG00000196611), entrezgene_id (e.g., 4312) or uniprot_gn_symbol (e.g., MMP1).", sep = ""));
			base::sink();
			cat("Unable to identify gene id type of ", unknown_ids, ". Make sure that all ids are one of ensembl_gene_id (e.g., ENSG00000196611), entrezgene_id (e.g., 4312) or uniprot_gn_symbol (e.g., MMP1).", "\n", sep = "");
			cat("The program: gene_drug_association.R failed while identifying type of gene ids. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}

    if(ppi_gene_id_type != drug_target_gene_id_type)
    {
      save(list = ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));
      error(file_logger, "--ppi_file and --drug_target_file have different types of gene identifier.");
      base::sink(); cat("The program: gene_drug_association.R failed because --ppi_file and --drug_target_file have different types of gene identifier.", "\n", sep = "");
      q(save = "no", status = 1);
    }else
    {
      # Define a variable to hold the common gene type 
      ref_gene_id_type <- ppi_gene_id_type;
    }

    ## Convert identifiers of input genes of interest, if required -------------------------------
    if( (gene_id_type != ref_gene_id_type) )
    {
      info(file_logger, paste("Gene identifier type in --genelist_file is different from --ppi_file and --drug_target_file", "\n", "Therefore, ids in --genelist_file will be converted from ", gene_id_type, " to ", ref_gene_id_type, "\n", sep = ""));

      if(ref_gene_id_type == "ensembl_gene_id")
      {
        ## Remove version information from all gene ids (in case of Ensembl, e.g., ENSG00000223972.5) --------
        if(ppi_gene_id_type == "ensembl_gene_id")
        {
          ppi_df[,'Protein_A'] <- gsub(as.character(ppi_df[,'Protein_A']), pattern = "\\..*", replacement = "");
          ppi_df[,'Protein_B'] <- gsub(as.character(ppi_df[,'Protein_B']), pattern = "\\..*", replacement = "");
        }

        if(drug_target_gene_id_type == "ensembl_gene_id")
        {
          drug_target_df$DrugTarget <- gsub(as.character(drug_target_df$DrugTarget), pattern = "\\..*", replacement = "");
        }
      }

      ## Get annotation of all genes using bioMart ----------------------------------------------------------
      #gene_info_df <- connectToBioMart(biomart_database_name = "ENSEMBL_MART_ENSEMBL", species = "Homo sapiens", file_logger);
      gene_info_df <- annotateWithBioMart(biomart_database_name = "ENSEMBL_MART_ENSEMBL", species = "Homo sapiens", gene_id_type, gene_ids = query_genelist, file_logger);

      if(nrow(gene_info_df) == 0)
      {
        save(list = ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));
        error(file_logger, "Unable to annotate any gene in the study using bioMart");
        base::sink(); cat(file_logger, "The program: gene_drug_association.R failed while annotating gene ids for all genes under study using bioMart.", "\n", sep = "");
        q(save = "no", status = 1);
      }

      info(file_logger, "Summary on gene annotation as obtained by biomaRt:");
      info(file_logger, paste("#Ensembl Gene Ids: ", length(unique(gene_info_df$ensembl_gene_id[! (is.na(gene_info_df$ensembl_gene_id) | (gene_info_df$ensembl_gene_id == "")) ])), seo = ""));
      info(file_logger, paste("#Entrez Gene Ids: ", length(unique(gene_info_df$entrezgene_id[! (is.na(gene_info_df$entrezgene_id) | (gene_info_df$entrezgene_id == "")) ])), sep = ""));
      info(file_logger, paste("#Uniprot Gene Symbols: ", length(unique(gene_info_df$uniprot_gn_symbol[! (is.na(gene_info_df$uniprot_gn_symbol) | (gene_info_df$uniprot_gn_symbol == "")) ])), sep = ""));
              
      info(file_logger, paste("Gene annotation table has ", nrow(gene_info_df), " entries.", sep = ""));
          
      write.table(x = gene_info_df, file = paste(output_dir, "/", "gene_info_from_BioMart.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);

      ## Creating a mapping between gene_id_type and ref_gene_id_type ----------------
      if(ref_gene_id_type == "uniprot_gn_symbol")
      {
        mapped_ids <- as.character(gene_info_df$uniprot_gn_symbol);
        names(mapped_ids) <- as.character(gene_info_df[, gene_id_type]);
      }else if(ref_gene_id_type == "ensembl_gene_id")
      {
        mapped_ids <- as.character(gene_info_df$ensembl_gene_id);
        names(mapped_ids) <- as.character(gene_info_df[, gene_id_type]);
      }else if(ref_gene_id_type == "entrezgene_id")
      {
        mapped_ids <- as.character(gene_info_df$entrezgene_id);
        names(mapped_ids) <- as.character(gene_info_df[, gene_id_type]);
      }else
      {
        save(list = ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));
        error(file_logger, "Unsupported gene identifier encountered! Currently only supported gene identifiers are: ensembl_gene_id, entrezgene_id and uniprot_gn_symbol.");
        base::sink(); cat(file_logger, "The program: gene_drug_association.R failed because an unsupported gene identifier encountered.", "\n", sep = "");
        q(save = "no", status = 1);
      }

      ## Removing NA or empty entries and duplicates in mapped_ids -------
      mapped_ids <- mapped_ids[! (is.na(mapped_ids) | (mapped_ids == "")) ];
      mapped_ids <- mapped_ids[!duplicated(mapped_ids)];
      mapped_ids <- mapped_ids[!duplicated(names(mapped_ids))];
      #print(mapped_ids);

      info(file_logger, paste("The 1-1 ", gene_id_type, " to ", ref_gene_id_type," mapping has ", length(mapped_ids), " entries", sep = ""));

      ## Update input genelist --------------------
      query_genelist <- mapped_ids[query_genelist];

      ## Remove NA
      query_genelist <- query_genelist[!is.na(names(query_genelist))];

      ## Remove genes entries ""
      query_genelist <- query_genelist[query_genelist != ""];

      info(file_logger, paste("Updated query_genelist has number of genes: ", length(query_genelist), sep = ""));
    }
  },
  warning =  function(w)
  {
    warn(file_logger, "Following warning encountered while checking gene identifier types and performing suitable conversion ...");
    print(w);
    base::sink(); cat("The program: gene_drug_association.R failed as it is unable to identify gene identifier type or to perform the required conversion.", "\n", sep = "");
    q(save = "no", status = 1);
  },
  error = function(e)
  {
    error(file_logger, "Following error encountered while checking gene identifier types and performing suitable conversion ..,");
    print(e);
    base::sink(); cat("The program: gene_drug_association.R failed as it is unable to identify gene identifier type or to perform the required conversion.", "\n", sep = "");
    q(save = "no", status = 1);
  }
)


###----------------------------------------------------------------------------------###
###       Overlap between query genes and PPI, and between PPI and drug targets
###----------------------------------------------------------------------------------###

info(file_logger, "Determining overlap between query genes and PPI, and between PPI and drug targets ...");

tryCatch(
  expr = 
  {
    ## Check overlap between query genes and PPI, and between PPI and drug targets --------------------
    net_edge1 <- ppi_df[,c('Protein_A','Protein_B')]                                         ### include columns those are required
    net_edge1 <- as.matrix(net_edge1)                                                        ### transform data into matrix
    mynetwork1 <- graph_from_data_frame(net_edge1, directed=FALSE, vertices = NULL)          ### Data-frame for Network (Un-Directed)

    mapped_genes_in_PPI <- intersect(V(mynetwork1)$name, query_genelist)                     ### query genes mapped in PPI
    mapped_drugtarget_in_PPI <- intersect(V(mynetwork1)$name, as.character(drug_target_df$DrugTarget))     ### drugtargets mapped in PPI

    if(length(mapped_genes_in_PPI) == 0)
    {
      save(list = ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));
      warn(file_logger, "Warning! No query gene mapped to PPI. Therefore, Gene-Drug association analysis will not proceed further.");
      base::sink(); cat("No query gene mapped to PPI. Therefore, the program: gene_drug_association.R will not proceed further.", "\n", sep = "");
      info(file_logger, "Processing ends");
      q(save = "no", status = 0);
    }

    if(length(mapped_drugtarget_in_PPI) == 0)
    {
      save(list = ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));
      warn(file_logger, "Warning! No drug target mapped to PPI. Therefore, Gene-Drug association analysis will not proceed further.");
      base::sink(); cat("No drug target mapped to PPI. Therefore, the program: gene_drug_association.R will not proceed further.", "\n", sep = "");
      info(file_logger, "Processing ends");
      q(save = "no", status = 0);
    }
  },
  warning = function(w)
  {
    warn(file_logger, "Following warning encountered while determining overlap between query genes and PPI, and between PPI and drug targets:");
    print(w);
    base::sink(); cat("The program: gene_drug_association.R failed while determining overlap between query genes and PPI, and between PPI and drug targets. Check log: ", log_file, " for more details.", "\n", sep = "");
    q(save = "no", status = 1);
  },
  error = function(e)
  {
    error(file_logger, "Following error encountered while determining overlap between query genes and PPI, and between PPI and drug targets:");
    print(e);
    base::sink(); cat("The program: gene_drug_association.R failed while determining overlap between query genes and PPI, and between PPI and drug targets. Check log: ", log_file, " for more details.", "\n", sep = "");
    q(save = "no", status = 1);
  }
)

###--------------------------------------------------------------------------###
###                                   OUTPUT
###--------------------------------------------------------------------------###

OutputFilename1 = paste(output_dir, "/", "Gene_Topology_attribute_in_PPI.csv", sep = "")
OutputFilename2 = paste(output_dir, "/", "QueryGene_DrugTarget_ShortestPath_matrix.csv", sep = "")
OutputImage1 = paste(output_dir, "/", "SubNet_PPI_QG_DT.html", sep = "")
OutputImage2 = paste(output_dir, "/", "Heatmap_QG_DT_ShortestPath.html", sep = "")

# OutputFilename1 <- "./Output/Gene_Topology_attribute_in_PPI.csv"
# OutputFilename2 <- "./Output/QueryGene_DrugTarget_ShortestPath_matrix.csv"


###--------------------------------------------------------------------------###
###       Overlap between Drug-targets and query genes(present in PPI)
###--------------------------------------------------------------------------###
info(file_logger, "Determining overlap between Drug-targets and query genes(present in PPI) ...");

tryCatch(
  expr =
  {
    com_gene_qgdt <- intersect(mapped_drugtarget_in_PPI,mapped_genes_in_PPI)
    non_overlap_dt <- mapped_drugtarget_in_PPI[!mapped_drugtarget_in_PPI %in% com_gene_qgdt]
    non_overlap_qg <- mapped_genes_in_PPI[!mapped_genes_in_PPI %in% com_gene_qgdt]

    uniq_dt_qg <- union(mapped_drugtarget_in_PPI,mapped_genes_in_PPI)
    other_proteins_in_PPI <- V(mynetwork1)$name[!V(mynetwork1)$name %in% uniq_dt_qg]   ### nodes other than DT, QG, DTQG

    dt_arr <- c()
    qg_arr <- c()
    qgdt_arr <- c()
    others_arr <- c()

    for (dt in 1:length(non_overlap_dt))
    {
      tmp1 <- paste(non_overlap_dt[dt], "DT", collapse = "\t")
      dt_arr <- append(dt_arr, tmp1)
    }
      
    #### Check if Query gene data comprises any column log2FC

    if(!is.null(log2fc_column) && log2fc_column %in% colnames(query_gene_df))
    {
      for (qg in 1:length(non_overlap_qg))
      {
        for (lfc in 1:length(query_genelist))
        {
          if (non_overlap_qg[qg] == query_genelist[lfc] && query_gene_df[,log2fc_column][lfc] > 0)
          {
            tmp2 <- paste(non_overlap_qg[qg], "QG_Up", collapse = "\t")
            qg_arr <- append(qg_arr, tmp2)
          }
          if (non_overlap_qg[qg] == query_genelist[lfc] && query_gene_df[,log2fc_column][lfc] < 0)
          {
            tmp2 <- paste(non_overlap_qg[qg], "QG_Down", collapse = "\t")
            qg_arr <- append(qg_arr, tmp2)
          }
        }
      }
      
      if (length(com_gene_qgdt)>0)
      {
        for (qgdt in 1:length(com_gene_qgdt))
        {
          for (lfc2 in 1:length(query_genelist))
          {
            if (com_gene_qgdt[qgdt] == query_genelist[lfc2] && query_gene_df[,log2fc_column][lfc2] > 0)
            {
              tmp3 <- paste(com_gene_qgdt[qgdt], "QGDT_Up", collapse = "\t")
              qgdt_arr <- append(qgdt_arr, tmp3)
            }
            if (com_gene_qgdt[qgdt] == query_genelist[lfc2] && query_gene_df[,log2fc_column][lfc2] < 0)
            {
              tmp3 <- paste(com_gene_qgdt[qgdt], "QGDT_Down", collapse = "\t")
              qgdt_arr <- append(qgdt_arr, tmp3)
            }
          }
        }
      }
    } else
    {
      for (qg in 1:length(non_overlap_qg))
      {
        tmp3 <- paste(non_overlap_qg[qg], "QG", collapse = "\t")
        qg_arr <- append(qg_arr, tmp3)
      }
      for (qgdt in 1:length(com_gene_qgdt))
      {
        tmp1 <- paste(com_gene_qgdt[qgdt], "QGDT", collapse = "\t")
        qgdt_arr <- append(qgdt_arr, tmp1)
      }
    }

    for (op in 1:length(other_proteins_in_PPI))
    {
      tmp4 <- paste(other_proteins_in_PPI[op], "Others", collapse = "\t")
      others_arr <- append(others_arr, tmp4)
    }

    combined_arr <- c(qgdt_arr,dt_arr,qg_arr,others_arr)
    df <- data.frame(combined_arr)
    com_arr_col <- df %>% extract(combined_arr, c("id", "gene_type"), "(.+) (.+)")
  },
  warning = function(w)
  {
    warn(file_logger, "Following warning encountered while determining overlap between Drug-targets and query genes(present in PPI).");
    print(w);
    base::sink(); cat("The program: gene_drug_association.R failed while determining overlap between Drug-targets and query genes(present in PPI). Check log: ", log_file, " for more details.", "\n", sep = "");
    q(save = "no", status = 1);
  },
  error = function(e)
  {
    error(file_logger, "Failed while determining overlap between Drug-targets and query genes(present in PPI) due to following reason: ");
    print(e);
    base::sink(); cat("The program: gene_drug_association.R failed while determining overlap between Drug-targets and query genes(present in PPI). Check log: ", log_file, " for more details.", "\n", sep = "");
    q(save = "no", status = 1);
  }
)

save(list = ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));


###--------------------------------------------------------------------------###
###                           Topology Calculation
###--------------------------------------------------------------------------###
info(file_logger, "Calculating topology ...");

tryCatch(
  expr = 
  {
    degcentrality <- degree(mynetwork1)
    eigencentrality <- evcent(mynetwork1)[[1]]
    closecentrality <- harmonic_centrality(mynetwork1, vids = V(mynetwork1), mode = "all", weights = NULL)

    net_Deg <- c()
    net_Eigen <- c()
    net_Close <- c()
    node_attr <- c()

    for (a in 1:length(V(mynetwork1)))
    {
      net_Deg <- append(net_Deg, degcentrality[which (names(degcentrality) %in% V(mynetwork1)$name[a])])
      net_Eigen <- append(net_Eigen, eigencentrality[which (names(eigencentrality) %in% V(mynetwork1)$name[a])])
      net_Close <- append(net_Close, closecentrality[which (names(closecentrality) %in% V(mynetwork1)$name[a])])
      node_attr <- append(node_attr, com_arr_col$gene_type[which (com_arr_col$id %in% V(mynetwork1)$name[a])])
    }

    df_NodeFeature <- cbind(V(mynetwork1)$name, net_Deg, net_Eigen, net_Close, node_attr)
    write.table(df_NodeFeature, file=OutputFilename1,sep=",", quote=F, row.names=F, col.names=c('id','Deg_Cen','Eigen_Cen', 'Close_Cen','gene_type'))
  },
  warning = function(w)
  {
    warn(file_logger, "Following warning encountered while calculating topology.");
    print(w);
    base::sink(); cat("The program: gene_drug_association.R failed while calculating topology. Check log: ", log_file, " for more details.", "\n", sep = "");
    q(save = "no", status = 1);
  },
  error = function(e)
  {
    error(file_logger, "Failed while calculating topology due to following reason: ");
    print(e);
    base::sink(); cat("The program: gene_drug_association.R failed while calculating topology. Check log: ", log_file, " for more details.", "\n", sep = "");
    q(save = "no", status = 1);
  }
)

save(list = ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));


###--------------------------------------------------------------------------###
###      Proximity between Drug-targets and query genes(present in PPI)
###--------------------------------------------------------------------------###
info(file_logger, "Calculating proximity between Drug-targets and query genes(present in PPI) ...");

tryCatch(
  expr = 
  {
    #dm1 <- shortest.paths(mynetwork1, v=mapped_genes_in_PPI, to=mapped_drugtarget_in_PPI)
    dm1 <- igraph::distances(graph = mynetwork1, v = mapped_genes_in_PPI, to = mapped_drugtarget_in_PPI, mode = "all", algorithm = "automatic");
        
    for (b in 1:length(colnames(dm1)))
    {
      for (c in 1:length(drug_target_df$DrugTarget))
      {
        if(drug_target_df$DrugTarget[c] == colnames(dm1)[b])
        {
          p <- paste(drug_target_df$DrugTarget[c], drug_target_df$DrugName[c], sep=".")
          colnames(dm1)[b] <- c(p)
        }
      }
    }
    write.table(cbind(gene_id = rownames(dm1), dm1), file=OutputFilename2, sep=",", quote=F, row.names=FALSE, col.names = TRUE)
  }, 
  warning = function(w)
  {
    warn(file_logger, "Following warning encountered while calculating proximity between Drug-targets and query genes(present in PPI).");
    print(w);
    base::sink(); cat("The program: gene_drug_association.R failed while calculating proximity between Drug-targets and query genes(present in PPI). Check log: ", log_file, " for more details.", "\n", sep = "");
    q(save = "no", status = 1);
  },
  error = function(e)
  {
    error(file_logger, "Failed while calculating proximity between Drug-targets and query genes(present in PPI) due to following reason: ");
    print(e);
    base::sink(); cat("The program: gene_drug_association.R failed while calculating proximity between Drug-targets and query genes(present in PPI). Check log: ", log_file, " for more details.", "\n", sep = "");
    q(save = "no", status = 1);
  }
)

save(list = ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));


###--------------------------------------------------------------------------###
###                         Network Visualization
###--------------------------------------------------------------------------###

info(file_logger, "Constructing Sub-Graph from genes and their neighbors for visualization ...");

tryCatch(
  expr =
  {
    nbr_qg <- ego(mynetwork1, order = 1, mapped_genes_in_PPI)         ### identify direct (order=1) neighbors of QG
    nbr_dt <- ego(mynetwork1, order = 1, mapped_drugtarget_in_PPI)    ### identify direct (order=1) neighbors of DT

    nb_qg_arr <- c()
    nb_dt_arr <- c()
    for (nbqg in 1:length(nbr_qg))
    {
      nb_qg_arr <- c(nb_qg_arr, nbr_qg[[nbqg]]$name)
    }
    for (nbdt in 1:length(nbr_dt))
    {
      nb_dt_arr <- c(nb_dt_arr, nbr_dt[[nbdt]]$name)
    }
    comm_nbr_arr <- intersect(nb_qg_arr,nb_dt_arr)
    full_node_arr <- c(comm_nbr_arr,mapped_genes_in_PPI,mapped_drugtarget_in_PPI)
    final_nodelist <- unique(full_node_arr)                           ### Final nodelist comprising QG, DT and their common neighbors


    subnet <- induced.subgraph(mynetwork1,final_nodelist)             ### making sub-graph
    subnet_edges <- as.data.frame(get.edgelist(subnet))               ### Transform igraph object to data frame

    gene_attr_table <- read.table(OutputFilename1, sep = ",", header = T)

    subnet_gene_attr <- data.frame(id=character(),Deg_Cen=integer(),Eigen_Cen=double(),Close_Cen=double(), gene_type = character(), stringsAsFactors=FALSE)
    for (x in 1:length(final_nodelist))
    {
      for(y in 1:length(gene_attr_table$id))
      {
        if (final_nodelist[x] == gene_attr_table$id[y])
        {
          subnet_gene_attr[x, ] <- gene_attr_table[y, ]
        }
      }
    }

    # nodes <- read.csv(".\\Output\\Gene_Topology_attribute_in_PPI.csv", header=T, as.is=T)
    nodes <- subnet_gene_attr
    edges <- data.frame(from = subnet_edges$V1, to = subnet_edges$V2)

    nodes$group <- nodes$gene_type          # group nodes
    nodes$size <- (nodes$Close_Cen/100)     # Node size
    nodes$title <- nodes$id                 # Text on click/hover
    nodes$label <- nodes$id                 # node label
    borderWidth <- 1                        # Node border width
    nodes$shadow <- TRUE                    # Nodes will drop shadow
    edges$hoverWidth <- 6
    edges$selectionWidth <- 6
    # edges$color <- "#DBCBF2"
    # edges$highlight <- "#3355FF"
    # edges$hover <- "#3355FF"

    if(!is.null(log2fc_column) && log2fc_column %in% colnames(query_gene_df))
    {
      final_net <- visNetwork(nodes, edges, height= "700px", width = "100%") %>%
        visInteraction(hover = TRUE, multiselect = TRUE) %>%
        visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T), selectedBy = "gene_type", nodesIdSelection = T, manipulation = T, collapse = T) %>%
        visGroups(groupname = "DT", color = list(background = "#00E673", border = "#00B33C", highlight = "yellow"), shape = "dot") %>%
        visGroups(groupname = "QG_Up", color = list(background = "#FF9580", border = "#FF6A4D", highlight = "yellow"), shape = "dot") %>%
        visGroups(groupname = "QG_Down", color = list(background = "#4DC3FF", border = "#0088CC", highlight = "yellow"), shape = "dot") %>%
        visGroups(groupname = "QGDT_Up", color = list(background = "#FF9580", border = "#FF6A4D", highlight = "yellow"), shape = "triangle") %>%
        visGroups(groupname = "QGDT_Down", color = list(background = "#4DC3FF", border = "#0088CC", highlight = "yellow"), shape = "triangle") %>%
        visGroups(groupname = "Others", color = list(background = "#FFFFE6", border = "#DBC658", highlight = "yellow"), shape = "dot") %>%
        visPhysics(enabled = T, solver = "forceAtlas2Based", forceAtlas2Based = list(gravitationalConstant = -20, avoidOverlap = 0)) %>%
        visLegend() %>%
        visLayout(randomSeed = 111)
    } else
    {
      final_net <- visNetwork(nodes, edges, height= "700px", width = "100%") %>%
        visInteraction(hover = TRUE, multiselect = TRUE) %>%
        visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T), selectedBy = "gene_type", nodesIdSelection = T, manipulation = T, collapse = T) %>%
        visGroups(groupname = "DT", color = list(background = "#00E673", border = "#00B33C", highlight = "yellow"), shape = "dot") %>%
        visGroups(groupname = "QG", color = list(background = "#FF9580", border = "#FF6A4D", highlight = "yellow"), shape = "dot") %>%
        visGroups(groupname = "QGDT", color = list(background = "#CC99FF", border = "#9933FF", highlight = "yellow"), shape = "triangle") %>%
        visGroups(groupname = "Others", color = list(background = "#FFFFE6", border = "#DBC658", highlight = "yellow"), shape = "dot") %>%
        visPhysics(enabled = T, solver = "forceAtlas2Based", forceAtlas2Based = list(gravitationalConstant = -20, avoidOverlap = 0)) %>%
        visLegend() %>%
        visLayout(randomSeed = 111)
    }

    visSave(final_net, file = OutputImage1, background = "white");
  },
  warning = function(w)
  {
    warn(file_logger, "Following warning encountered while generating network visualization.");
    print(w);
    base::sink(); cat("The program: gene_drug_association.R failed while generating network visualization. Check log: ", log_file, " for more details.", "\n", sep = "");
    q(save = "no", status = 1);
  },
  error = function(e)
  {
    error(file_logger, "Failed at network visualization step due to following reason: ");
    print(e);
    base::sink(); cat("The program: gene_drug_association.R failed while generating network visualization. Check log: ", log_file, " for more details.", "\n", sep = "");
    q(save = "no", status = 1);
  }
)

save(list = ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));


###--------------------------------------------------------------------------###
###                          SP matrix Visualization
###--------------------------------------------------------------------------###
info(file_logger, "Generating shortest path matrix visualization ...");

tryCatch(
  expr = 
  {
    ### checking for infinite values, remove them (for heatmap only)
    temp <- c()
    for (x in 1:nrow(dm1))
    {
      if(any(is.infinite(dm1[x,])))
      {
        temp = c(temp,x)
      }
    }

    if (length(temp) == 0)
    {
      heatmaply(dm1, key.title="Shortest Path", xlab = "Target.Drugname", ylab = "QueryGene", file = OutputImage2)
    } else
    {
      dm2 <- dm1[-temp, ]
      heatmaply(dm2, key.title="Shortest Path", xlab = "Target.Drugname", ylab = "QueryGene", file = OutputImage2)
    }
  },
  warning = function(w)
  {
    warn(file_logger, "Following warning encountered while generating shortest path matrix visualization.");
    print(w);
    base::sink(); cat("The program: gene_drug_association.R failed while generating shortest path matrix visualization. Check log: ", log_file, " for more details.", "\n", sep = "");
    q(save = "no", status = 1);
  },
  error = function(e)
  {
    error(file_logger, "Failed at shortest path matrix visualization step due to following reason: ");
    print(e);
    base::sink(); cat("The program: gene_drug_association.R failed while generating shortest path matrix visualization. Check log: ", log_file, " for more details.", "\n", sep = "");
    q(save = "no", status = 1);
  }
)

save(list = ls(), file = paste(output_dir, "/", "GeneDrugAssociation.RData", sep = ""));

info(file_logger, "Processing ends");
info(file_logger, "R session info:");

## Store R session information -----------------------------
print(sessionInfo());

## Store citation of the major packages involved -----------
current_Rsession <- sessionInfo();
print(current_Rsession);
if( length(current_Rsession$otherPkgs) > 0 )
{
	other_packages <- unlist(lapply(X = current_Rsession$otherPkgs, FUN = function(x){x$Package}));
	for(pkg in other_packages)
	{
		print(citation(package = pkg), bibtex = FALSE);
	}
}

base::sink();
cat("The program: gene_drug_association.R executed successfully!", "\n", sep = "");




