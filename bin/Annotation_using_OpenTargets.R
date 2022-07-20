#!/usr/bin/Rscript --vanilla

#########################################
## Program: Annotation_using_OpenTargets.R
## Purpose: Annotate a set of input genes against 'Open Targets' database using GraphQL queries.
##
## Usage: /path/to/Rscript --vanilla /path/to/Annotation_using_OpenTargets.R --genelist_file <path/to/gene list> --gene_id_column <column number for gene identifiers, e.g., 1> --output_dir <path/to/output directory> --output_prefix <prefix for output files such as figures> --log_file <filename for log>
##
## Developer: Srikant Verma (srikant_verma@persistent.com)
## Version: 1.0: Created on June 9, 2022
## Version: 1.1: Updated on June 12, 2022
## Version: 1.2: Updated on June 28, 2022
#########################################

## Version: 1.0: Created on June 9, 2022
##  It performs the following major tasks:
#   1. Queries 'Open Targets' database for known drug information using GraphQL queries
#
#
## Version: 1.1: Updated on June 12, 2022
#	Enhacements:
#		1. Improved processing time for annotation (fetching known drugs information) by making a single API call for entire genelist. Earlier, I was making API call for each gene separately.
#	Bug fixes:
#		1. variable 'reattempt' was not getting updated in 'error' section of tryCatch using the assignment operator '<-'. Now, we are using superassignment operator (<<-) instead to do the job.
#
#
## Version: 1.2: Updated on June 28, 2022
#	Enhancements/minor issue fix:
#		1. Send more error messages to STDOUT (which were previously directed to log files only).
#		2. Handle NULL values to parameters (such as --output_prefix) when it is provided as command line arguments to these parameters
#		3. Removd save(list=ls(), file = ...) from warning and error blocks of tryCatch() because it only saves warning/error message.
#	Bug fixes:
#		1. Fixed Ensembl gene id pattern in identifyGeneIdType() as described here: https://asia.ensembl.org/info/genome/stable_ids/index.html. Stable Ensembl IDs have the following format: ENS[species prefix][feature type prefix][a unique eleven digit number].
#		2. Fixed Uniport gene symbol pattern in identifyGeneIdType() as described here: https://www.genenames.org/about/old-guidelines/
#		3. Fixed Entrez gene id pattern in identifyGeneIdType()
#
# Known issues:

## Future releases:

## Notes:

## References:
#   1. Platform Browser API: https://api.platform.opentargets.org/api/v4/graphql/browser
#	2. Crash course: https://blog.opentargets.org/summer-school-crash-course-part-1/
#

##########################################

cat("For help: /path/to/Rscript --vanilla /path/to/Annotation_using_OpenTargets.R --help", "\n", sep = "");

## =========================================================
## Take user specified parameters
## =========================================================
library("optparse");	# command line option parser

args_list <- list(
    optparse::make_option(opt_str = c("--genelist_file"),	type = "character",	default = NULL,									help = "List of genes in CSV or TSV format. Default: NULL"),
	optparse::make_option(opt_str = c("--gene_id_column"),	type = "integer",	default = 1,									help = "Column number in --genelist_file which contains gene identifiers. Default: 1"),
	optparse::make_option(opt_str = c("--output_dir"),		type = "character",	default = "./",									help = "Output directory path. Default: ./"),
	optparse::make_option(opt_str = c("--output_prefix"),	type = "character",	default = NULL,									help = "Prefix for output files. Default: NULL"),
	optparse::make_option(opt_str = c("--log_file"),		type = "character",	default = "Annotation_with_OpenTargets.log",	help = "Log file. Default: Annotation_with_OpenTargets.log")
)
arguments <- optparse::parse_args(optparse::OptionParser(option_list = args_list));


## =========================================================
## Load required libraries
## =========================================================
library("typed");		    # To set types for a variable, and for arguments of a function
library("log4r");		    # Logging library
library("data.table");	    # Read input data fast
library("dplyr");			# For data manipulation
library("biomaRt");			# Provides annotation for genes
library("httr");			# For HTTP requests
library("jsonlite");        # To process JSON data

## =========================================================
## Set imp. variables
## =========================================================
declare("genelist_file", Character(length = 1, null_ok = FALSE), value = arguments$genelist_file, const = TRUE);
declare("gene_id_column", Integer(length = 1, null_ok = FALSE), value = arguments$gene_id_column, const = TRUE);
declare("output_dir", Character(length = 1, null_ok = FALSE), value = arguments$output_dir, const = TRUE);

if(!is.null(arguments$output_prefix) && (arguments$output_prefix == 'NULL'))
{
	arguments$output_prefix <- NULL;
}
declare("output_prefix", Character(length = 1, null_ok = TRUE), value = arguments$output_prefix, const = FALSE);
output_prefix <- ifelse(test = is.null(output_prefix), yes = "", no = paste(output_prefix, "_", sep = ""));

declare("log_file", Character(length = 1, null_ok = FALSE), value = arguments$log_file, const = FALSE);

declare("genelist_df", Data.frame());								# Gene list data - frame
declare("annotation_using_OpenTargets_df", Data.frame());


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

# =========================================================
## Custom functions
## ========================================================

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
									traceback();
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
				traceback();
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

	return(gene_info_df);
}



## Query 'Open Targets' database to list drugs with investigational or approved indications targeting the query gene according to their curated mechanism of action
queryTargetEndpointForKnownDrugs <- Any() ? function(target_gene_id = Character(), file_logger = Any())
{
	declare("base_url", Character());
	declare("query_string_to_determine_size", Character());
	declare("final_query_string", Character());
	declare("variables", List());
	declare("post_body", List());
	declare("r", Any());
	declare("request_status_code", Integer());
	declare("reattempt", Integer());

	# Set base URL of GraohQL API endpoint
	base_url <- "https://api.platform.opentargets.org/api/v4/graphql";

	debug(file_logger, paste("Processing gene id: ", target_gene_id, " ...", sep = ""));

	## This query string helps to know the number of entries to retrieve
	query_string_to_determine_size <- "
		query find_target_information($ensemblId: String!) {
			target(ensemblId: $ensemblId) {
				id
				approvedSymbol
				knownDrugs {
					count
				}
			}
		}
	"

	# Set 'variables' object of arguments to be passed to endpoint
	variables <- list("ensemblId" = target_gene_id);

	# Construct POST request body object with query string and variables
	post_body <- list(query = query_string_to_determine_size, variables = variables);

	# Perform REST request
	debug(file_logger, "POST request to determine the number of entries to retrieve ... ");
	r <- NULL;
	request_status_code <- as.integer(0);
	reattempt <- as.integer(0);
	while(request_status_code != 200)
	{
		tryCatch(
			expr =
			{
				if(reattempt > 0)
				{
					debug(file_logger, paste("Re-attempting POST request #: ", reattempt, sep = ""));
				}
				r <- POST(url = base_url, body = post_body, encode = 'json');
				request_status_code <- r$status_code;
			},
			error = function(e)
			{
				error(file_logger, "Following error encountered while retrieving data from 'Open Targets'.");
				print(e);

				if(reattempt < 10)
				{
					debug(file_logger, "Will retry ...");
					reattempt <<- as.integer(reattempt + 1);
				}else
				{
					fatal(file_logger, "Giving up now! Program will terminate now!");
					base::sink();
					print(e);
					cat("The program: Annotation_using_OpenTargets.R failed while retrieving data from the database.", "\n", sep = "");
					q(save = "no", status = 1);
				}
			}
		)
	}


	if(is.null(content(r)$data$target$knownDrugs$count))
	{
		debug(file_logger, paste("Total entries: ", 0, sep = ""));
		return(r);
	}else
	{
		debug(file_logger, paste("Total entries: ", content(r)$data$target$knownDrugs$count, sep = ""));
	}

	## This query string helps to download all data
	final_query_string <- "
		query find_target_information($ensemblId: String!, $size: Int!) {
			target(ensemblId: $ensemblId) {
				id
				approvedSymbol
				knownDrugs (size: $size) {
					uniqueDrugs
					uniqueTargets
					uniqueDiseases
					count
					rows {
						approvedName
						approvedSymbol
						drugId
						prefName
						drugType
						mechanismOfAction
						diseaseId
						label
						phase
						status
						urls {
							url
							name
						}
					}
				}
			}
		}
	"

	# Update 'variables' object of arguments to be passed to endpoint
	variables <- list("ensemblId" = target_gene_id, size = content(r)$data$target$knownDrugs$count);

	# Update POST request body object with query string and variables
	post_body <- list(query = final_query_string, variables = variables);

	# Perform updated REST request
	debug(file_logger, "POST request to fetch all relevant data ... ");
	r <- NULL;
	request_status_code <- as.integer(0);
	reattempt <- as.integer(0);
	while(request_status_code != 200)
	{
		tryCatch(
			expr =
			{
				if(reattempt > 0)
				{
					debug(file_logger, paste(paste("Re-attempting POST request #: ", reattempt, sep = "")));
				}
				r <- POST(url = base_url, body = post_body, encode = 'json');
				request_status_code <- r$status_code;
			},
			error = function(e)
			{
				error(file_logger, "Following error encountered while retrieving data from 'Open Targets'.");
				print(e);

				if(reattempt < 10)
				{
					debug(file_logger, "Will retry ...");
					reattempt <<- as.integer(reattempt + 1);
				}else
				{
					fatal(file_logger, "Giving up now! Program will terminate now!");
					base::sink();
					print(e);
					cat("The program: Annotation_using_OpenTargets.R failed while retrieving data from the database.", "\n", sep = "");
					q(save = "no", status = 1);
				}
			}
		)
	}
	return(r);
}


## Annotate genes with information on known drugs using 'Open Targets' database. Here, each query consists of a single gene only.
annotateUsingOpenTargetsViaMultipleAPIcalls <- Data.frame() ? function(ensembl_gene_ids = Character(), file_logger = Any())
{
    ## Store query results for all genes of interest

	declare("result_list", List(null_ok = TRUE));
	declare("annotation_df", Data.frame());

	for(gene_id in ensembl_gene_ids)
	{
		declare("r", Any());
		r <- queryTargetEndpointForKnownDrugs(target_gene_id = gene_id, file_logger);
		result_list[[gene_id]] <- content(r)$data;
	}

	annotation_df <- data.frame(id = NA, approved_name = NA, approved_symbol = NA, drug_id = NA, pref_name = NA, drug_type = NA, mechanism_of_action = NA, disease_id = NA, label = NA, phase = NA, status = NA, urls = NA);

	info(file_logger, "Here is the summary on annotation: ");

	for(gene_id in ensembl_gene_ids)
	{
		declare("current_id", Character(null_ok = TRUE));
		declare("current_approved_symbol", Character(null_ok = TRUE));
		declare("total_entries", Integer(null_ok = TRUE));
		declare("unique_drugs", Integer(null_ok = TRUE));
		declare("unique_diseases", Integer(null_ok = TRUE));

		current_id <- result_list[[gene_id]]$target$id;
		if(is.null(current_id))
		{
			cat("Could not find any entry for gene: ", gene_id, " in Open Targets database.", "\n", sep = "");
			next
		}
		current_approved_symbol <- result_list[[gene_id]]$target$approvedSymbol;
		if(is.null(current_approved_symbol)){ current_approved_symbol <- "" };
		total_entries <- result_list[[gene_id]]$target$knownDrugs$count;
		if(is.null(total_entries)){ total_entries <- as.integer(0)};
		unique_drugs <- result_list[[gene_id]]$target$knownDrugs$uniqueDrugs;
		if(is.null(unique_drugs)){ unique_drugs <- as.integer(0)};
		unique_diseases <- result_list[[gene_id]]$target$knownDrugs$uniqueDiseases;
		if(is.null(unique_diseases)){ unique_diseases <- as.integer(0)};

		cat("gene id: ", gene_id, "\t", "symbol: ", current_approved_symbol, "\n", sep = "");
		cat("#unique drugs: ", unique_drugs, "\n", sep = "");
		cat("#unique diseases: ", unique_diseases, "\n", sep = "");
		cat("Total entries: ", total_entries, "\n", sep = "");

		declare("temp_df", Data.frame());

		if(total_entries)
		{
			temp_df <- data.frame(
				id = result_list[[gene_id]]$target$id,
				approved_name = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ approved_name <- result_list[[gene_id]]$target$knownDrugs$rows[[x]]$approvedName; ifelse(test = {is.null(approved_name)}, yes = "", no = approved_name) }))),
				approved_symbol = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ approved_symbol <- result_list[[gene_id]]$target$knownDrugs$rows[[x]]$approvedSymbol; ifelse(test = {is.null(approved_symbol)}, yes = "", no = approved_symbol) }))),
				drug_id = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ drug_id <- result_list[[gene_id]]$target$knownDrugs$rows[[x]]$drugId; ifelse(test = {is.null(drug_id)}, yes = "", no = drug_id) }))),
				pref_name = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ pref_name <- result_list[[gene_id]]$target$knownDrugs$rows[[x]]$prefName; ifelse(test = {is.null(pref_name)}, yes = "", no = pref_name) }))),
				drug_type = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ drug_type <- result_list[[gene_id]]$target$knownDrugs$rows[[x]]$drugType; ifelse(test = {is.null(drug_type)}, yes = "", no = drug_type) }))),
				mechanism_of_action = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ mechanism_of_action <- result_list[[gene_id]]$target$knownDrugs$rows[[x]]$mechanismOfAction; ifelse(test = {is.null(mechanism_of_action)}, yes = "", no = mechanism_of_action) }))),
				disease_id = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ disease_id <- result_list[[gene_id]]$target$knownDrugs$rows[[x]]$diseaseId; ifelse(test = {is.null(disease_id)}, yes = "", no = disease_id) }))),
				label = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ label <- result_list[[gene_id]]$target$knownDrugs$rows[[x]]$label; ifelse(test = {is.null(label)}, yes = "", no = label) }))),
				phase = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ phase <- result_list[[gene_id]]$target$knownDrugs$rows[[x]]$phase; ifelse(test = {is.null(phase)}, yes = "", no = phase) }))),
				status = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ status <- result_list[[gene_id]]$target$knownDrugs$rows[[x]]$status; ifelse(test = {is.null(status)}, yes = "", no = status) }))),
				urls = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ urls <- paste(unlist(lapply(X = result_list[[gene_id]]$target$knownDrugs$rows[[x]]$urls, FUN = function(xx){ xx$url })), collapse = ",") }))),
				stringsAsFactors = FALSE
			);
		}else
		{
			temp_df <- data.frame(id = current_id, approved_name = current_approved_symbol, approved_symbol = "", drug_id = "", pref_name = "", drug_type = "", mechanism_of_action = "", disease_id = "", label = "", phase = "", status = "", urls = "", stringsAsFactors = FALSE);
		}
		annotation_df <- rbind(annotation_df, temp_df);
	}

	annotation_df <- annotation_df[-1,];

	return(annotation_df);
}


## Annotate genes with information on known drugs using 'Open Targets' database. Here, a single query is used for entire genelist.
annotateUsingOpenTargetsViaSingleAPIcall <- Data.frame() ? function(ensembl_gene_ids = Character(), file_logger = Any())
{
	declare("base_url", Character());
	declare("query_string_to_determine_size", Character());
	declare("final_query_string", Character());
	declare("variables", List());
	declare("post_body", List());
	declare("r", Any());
	declare("request_status_code", Integer());
	declare("reattempt", Integer());

	# Set base URL of GraohQL API endpoint
	base_url <- "https://api.platform.opentargets.org/api/v4/graphql";


	## Determine the max. number of entries to retrieve for each gene  -----------------------------

	debug(file_logger, "POST request to determine the max. number of entries to retrieve for each gene ... ");

	## This query string helps to know the max. number of entries that will suffice to get all relevant data for each gene
	query_string_to_determine_size <- "
		query find_target_information($ensemblIds: [String!]!) {
			targets(ensemblIds: $ensemblIds) {
				id
				approvedSymbol
				knownDrugs {
					count
				}
			}
		}
	"

	# Set 'variables' object of arguments to be passed to endpoint
	variables <- list("ensemblIds" = ensembl_gene_ids);

	# Construct POST request body object with query string and variables
	post_body <- list(query = query_string_to_determine_size, variables = variables);

	# Perform REST request
	r <- NULL;
	request_status_code <- as.integer(0);
	reattempt <- as.integer(0);
	while(request_status_code != 200)
	{
		tryCatch(
			expr =
			{
				if(reattempt > 0)
				{
					debug(file_logger, paste("Re-attempting POST request #: ", reattempt, sep = ""));
				}
				r <- POST(url = base_url, body = post_body, encode = 'json');
				request_status_code <- r$status_code;
			},
			error = function(e)
			{
				error(file_logger, "Following error encountered while retrieving data from 'Open Targets'.");
				print(e);
				traceback();

				if(reattempt < 10)
				{
					debug(file_logger, "Will retry ...");
					reattempt <<- as.integer(reattempt + 1);
				}else
				{
					fatal(file_logger, "Giving up now! Program will terminate now!");
					base::sink();
					print(e);
					cat("The program: Annotation_using_OpenTargets.R failed while retrieving data from the database.", "\n", sep = "");
					q(save = "no", status = 1);
				}
			}
		)
	}

	max_page_size <- max(unlist(invisible(lapply(X = 1:length(content(r)$data$targets), FUN = function(x){temp_size <- content(r)$data$targets[[x]]$knownDrugs$count; ifelse(test = {is.null(temp_size)}, yes = 0, no = temp_size)}))));

	debug(file_logger, paste("Maximum number of entries for a gene is: ", max_page_size, sep = ""));

	if(max_page_size == 0)
	{
		info(file_logger, "There is no information on known drugs for any of the input genes. Hence program will end here.");
		base::sink(); cat("There is no information on known drugs for any of the input genes. Hence program will end here.", "\n", sep = "");
		q(save = "no", status = 0);

	}

	## Fetch all relevant data for the input set of genes -----------------------------
	debug(file_logger, "POST request to fetch all relevant data ... ");

	## This query string helps to download all data
	final_query_string <- "
		query find_target_information($ensemblIds: [String!]!, $size: Int!) {
			targets(ensemblIds: $ensemblIds) {
				id
				approvedSymbol
				knownDrugs (size: $size) {
					uniqueDrugs
					uniqueTargets
					uniqueDiseases
					count
					rows {
						approvedName
						approvedSymbol
						drugId
						prefName
						drugType
						mechanismOfAction
						diseaseId
						label
						phase
						status
						urls {
							url
							name
						}
					}
				}
			}
		}
	"

	# Update 'variables' object of arguments to be passed to endpoint
	variables <- list("ensemblIds" = ensembl_gene_ids, size = max_page_size);

	# Update POST request body object with query string and variables
	post_body <- list(query = final_query_string, variables = variables);

	# Perform updated REST request
	r <- NULL;
	request_status_code <- as.integer(0);
	reattempt <- as.integer(0);
	while(request_status_code != 200)
	{
		tryCatch(
			expr =
			{
				if(reattempt > 0)
				{
					debug(file_logger, paste(paste("Re-attempting POST request #: ", reattempt, sep = "")));
				}
				r <- POST(url = base_url, body = post_body, encode = 'json');
				request_status_code <- r$status_code;
			},
			error = function(e)
			{
				error(file_logger, "Following error encountered while retrieving data from 'Open Targets'.");
				print(e);
				traceback();

				if(reattempt < 10)
				{
					debug(file_logger, "Will retry ...");
					reattempt <<- as.integer(reattempt + 1);
				}else
				{
					fatal(file_logger, "Giving up now! Program will terminate now!");
					base::sink();
					print(e);
					cat("The program: Annotation_using_OpenTargets.R failed while retrieving data from the database.", "\n", sep = "");
					q(save = "no", status = 1);
				}
			}
		)
	}

	info(file_logger, paste(length(content(r)$data$targets), "/", length(ensembl_gene_ids), " genes could be mapped to 'Open Targets database.'", sep = ""));

	## Store query results for all genes of interest in a data frame -------------------
	annotation_df <- data.frame(id = NA, approved_name = NA, approved_symbol = NA, drug_id = NA, pref_name = NA, drug_type = NA, mechanism_of_action = NA, disease_id = NA, label = NA, phase = NA, status = NA, urls = NA);

	info(file_logger, "Here is the summary on annotation: ");

	for(i in 1:length(content(r)$data$targets))
	{

		declare("current_id", Character(null_ok = TRUE));
		declare("current_approved_symbol", Character(null_ok = TRUE));
		declare("total_entries", Integer(null_ok = TRUE));
		declare("unique_drugs", Integer(null_ok = TRUE));
		declare("unique_diseases", Integer(null_ok = TRUE));

		declare("current_target_ann", List(null_ok = TRUE));
		current_target_ann <- content(r)$data$targets[[i]];
		#names(current_target_ann);	#[1] "id"             "approvedSymbol" "knownDrugs"

		current_id <- current_target_ann$id;
		current_approved_symbol <- current_target_ann$approvedSymbol;
		if(is.null(current_approved_symbol)){ current_approved_symbol <- "" };
		total_entries <- current_target_ann$knownDrugs$count;
		if(is.null(total_entries)){ total_entries <- as.integer(0)};
		unique_drugs <- current_target_ann$knownDrugs$uniqueDrugs;
		if(is.null(unique_drugs)){ unique_drugs <- as.integer(0)};
		unique_diseases <- current_target_ann$knownDrugs$uniqueDiseases;
		if(is.null(unique_diseases)){ unique_diseases <- as.integer(0)};

		cat("gene id: ", current_id, "\t", "symbol: ", current_approved_symbol, "\n", sep = "");
		cat("#unique drugs: ", unique_drugs, "\n", sep = "");
		cat("#unique diseases: ", unique_diseases, "\n", sep = "");
		cat("Total entries: ", total_entries, "\n", sep = "");

		declare("temp_df", Data.frame());
		if(total_entries)
		{
			temp_df <- data.frame(
				id = current_target_ann$id,
				approved_name = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ approved_name <- current_target_ann$knownDrugs$rows[[x]]$approvedName; ifelse(test = {is.null(approved_name)}, yes = "", no = approved_name) }))),
				approved_symbol = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ approved_symbol <- current_target_ann$knownDrugs$rows[[x]]$approvedSymbol; ifelse(test = {is.null(approved_symbol)}, yes = "", no = approved_symbol) }))),
				drug_id = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ drug_id <- current_target_ann$knownDrugs$rows[[x]]$drugId; ifelse(test = {is.null(drug_id)}, yes = "", no = drug_id) }))),
				pref_name = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ pref_name <- current_target_ann$knownDrugs$rows[[x]]$prefName; ifelse(test = {is.null(pref_name)}, yes = "", no = pref_name) }))),
				drug_type = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ drug_type <- current_target_ann$knownDrugs$rows[[x]]$drugType; ifelse(test = {is.null(drug_type)}, yes = "", no = drug_type) }))),
				mechanism_of_action = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ mechanism_of_action <- current_target_ann$knownDrugs$rows[[x]]$mechanismOfAction; ifelse(test = {is.null(mechanism_of_action)}, yes = "", no = mechanism_of_action) }))),
				disease_id = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ disease_id <- current_target_ann$knownDrugs$rows[[x]]$diseaseId; ifelse(test = {is.null(disease_id)}, yes = "", no = disease_id) }))),
				label = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ label <- current_target_ann$knownDrugs$rows[[x]]$label; ifelse(test = {is.null(label)}, yes = "", no = label) }))),
				phase = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ phase <- current_target_ann$knownDrugs$rows[[x]]$phase; ifelse(test = {is.null(phase)}, yes = "", no = phase) }))),
				status = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ status <- current_target_ann$knownDrugs$rows[[x]]$status; ifelse(test = {is.null(status)}, yes = "", no = status) }))),
				urls = invisible(unlist(lapply(X = 1:total_entries, FUN = function(x){ urls <- paste(unlist(lapply(X = current_target_ann$knownDrugs$rows[[x]]$urls, FUN = function(xx){ xx$url })), collapse = ",") }))),
				stringsAsFactors = FALSE
			);
		}else
		{
			temp_df <- data.frame(id = current_id, approved_name = current_approved_symbol, approved_symbol = "", drug_id = "", pref_name = "", drug_type = "", mechanism_of_action = "", disease_id = "", label = "", phase = "", status = "", urls = "", stringsAsFactors = FALSE);
		}
		annotation_df <- rbind(annotation_df, temp_df);
	}
	annotation_df <- annotation_df[-1,];

	return(annotation_df);
}

## =======================================================
## Read input list of genes
## =======================================================
tryCatch(
    expr =
    {
        ## Check read access to input file ---------------------
		if(file.access(names = genelist_file, mode = 4) != 0)
		{
			error(file_logger, paste("Unable to read --genelist_file file: ", genelist_file, sep = ""));
			base::sink(); cat("The program: Annotation_using_OpenTargets.R failed because it was unable to read --genelist_file. Check log for more details: ", log_file, "\n", sep = "");
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
			cat("The program: Annotation_using_OpenTargets.R failed while reading input list of genes. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}

		## Store genelist -------------------------------------
		genelist <- as.character(genelist_df[, gene_id_column]);
    },
    warning = function(w)
    {
        warn(file_logger, "Caught a warning while reading input list of genes. Something unexpected must have happened.");
		print(w);
		base::sink();
		print(w);
		cat("The program: Annotation_using_OpenTargets.R failed at 'Read input list of genes'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
    },
    error = function(e)
    {
        error(file_logger, "An error occured while reading input list of genes.");
		print(e);
		base::sink();
		print(e);
		cat("The program: Annotation_using_OpenTargets.R failed at 'Read input list of genes'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
    }
)

save(list=ls(), file = paste(output_dir, "/", "AnnotationOpenTargets.RData", sep = ""));


## =====================================================================================
##  Check gene identifier type, and perform conversion to Ensembl gene id, if required
## =====================================================================================
tryCatch(
	expr =
	{
		debug(file_logger, paste("Total #genes in the current study: ", length(genelist), sep = ""));

		gene_id_type <- identifyGeneIdType(ids = genelist);
		info(file_logger, paste("Inferred gene identifier type in the current input list of genes is: ", gene_id_type, sep = ""));

		if(gene_id_type == "Unknown")
		{
			save(list=ls(), file = paste(output_dir, "/", "AnnotationOpenTargets.RData", sep = ""));
			error(file_logger, "Unable to identify gene id type of input genes. Make sure that all ids are one of ensembl_gene_id (e.g., ENSG00000196611), entrezgene_id (e.g., 4312) or uniprot_gn_symbol (e.g., MMP1). Mixed types not allowed.");
			base::sink();
			cat("Unable to identify gene id type of input genes. Make sure that all ids are one of ensembl_gene_id (e.g., ENSG00000196611), entrezgene_id (e.g., 4312) or uniprot_gn_symbol (e.g., MMP1). Mixed types not allowed.", "\n", sep = "");
			cat("The program: Annotation_using_OpenTargets.R failed while identifying type of gene ids. Check log for more details: ", log_file, "\n", sep = "");
			q(save = "no", status = 1);
		}

		## Remove version information from all gene ids (in case of Ensembl, e.g., ENSG00000223972.5) --------
		if(gene_id_type == "ensembl_gene_id")
		{
			genelist <- gsub(genelist, pattern = "\\..*", replacement = "");
		}


		if(gene_id_type != "ensembl_gene_id")
		{
			info(file_logger, "'Open Targets' requires Ensembl gene identifiers. Therefore, we will map the current identifier to Ensembl gene id using biomaRt.");

			## Create a mapping for all genes to Ensembl gene ids, if required --------------------
			declare("ensembl", Any());
			declare("gene_info_df", Data.frame());
			declare("mapped_ensembl_gene_ids", Character());

			## Get annotation of all genes using bioMart ---------------------------
			info(file_logger, "Preparing annotation table for all genes using bioMart ...");
			gene_info_df <- annotateWithBioMart(biomart_database_name = "ENSEMBL_MART_ENSEMBL", species = "Homo sapiens", gene_id_type, gene_ids = genelist, file_logger);

			if(nrow(gene_info_df) == 0)
			{
				save(list=ls(), file = paste(output_dir, "/", "AnnotationOpenTargets.RData", sep = ""));
				error(file_logger, "Unable to annotate the input list of genes.");
				base::sink();
				cat("Unable to annotate the input list of genes.", "\n", sep = "");
				cat("The program: Annotation_using_OpenTargets.R failed while annotating gene ids for all genes under study using bioMart. Check log for more details: ", log_file, "\n", sep = "");
				q(save = "no", status = 1);
			}

			info(file_logger, "Summary on gene annotation as obtained by biomaRt:");
			cat("#Ensembl Gene Ids: ", length(unique(gene_info_df$ensembl_gene_id[ !(is.na(gene_info_df$ensembl_gene_id) | (gene_info_df$ensembl_gene_id == "")) ])), "\n");
			cat("#Entrez Gene Ids: ", length(unique(gene_info_df$entrezgene_id[ !(is.na(gene_info_df$entrezgene_id) | (gene_info_df$entrezgene_id == "")) ])), "\n");
			cat("#Uniprot Gene Symbols: ", length(unique(gene_info_df$uniprot_gn_symbol[ !(is.na(gene_info_df$uniprot_gn_symbol) | (gene_info_df$uniprot_gn_symbol == "")) ])), "\n");

			## There could be scenarios where there is no 1-1 mapping between Ensembl and Entrez.
			info(file_logger, paste("Gene annotation table has ", nrow(gene_info_df), " entries.", sep = ""));

			write.table(x = gene_info_df, file = paste(output_dir, "/", output_prefix, "gene_info_from_BioMart.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);

			if(gene_id_type == "entrezgene_id")
			{
				mapped_ensembl_gene_ids <- as.character(gene_info_df$ensembl_gene_id);
				names(mapped_ensembl_gene_ids) <- gene_info_df$entrezgene_id;

				#Removing duplicate or NA entries in 'ensembl_gene_id'
				mapped_ensembl_gene_ids <- mapped_ensembl_gene_ids[!(is.na(mapped_ensembl_gene_ids) | (mapped_ensembl_gene_ids == ""))];
				mapped_ensembl_gene_ids <- mapped_ensembl_gene_ids[!duplicated(mapped_ensembl_gene_ids)];

				#There might still be duplicate entrez gene ids in names(mapped_ensembl_gene_ids)
				mapped_ensembl_gene_ids <- mapped_ensembl_gene_ids[!duplicated(names(mapped_ensembl_gene_ids))];

				info(file_logger, paste("The 1-1 entrezgene_id-ensembl_gene_id mapping has ", length(mapped_ensembl_gene_ids), " entries", sep = ""));
			}else if(gene_id_type == "uniprot_gn_symbol")
			{
				mapped_ensembl_gene_ids <- as.character(gene_info_df$ensembl_gene_id);
				names(mapped_ensembl_gene_ids) <- gene_info_df$uniprot_gn_symbol;

				#Removing duplicate or NA entries in 'ensembl_gene_id'
				mapped_ensembl_gene_ids <- mapped_ensembl_gene_ids[!(is.na(mapped_ensembl_gene_ids) | (mapped_ensembl_gene_ids == ""))];
				mapped_ensembl_gene_ids <- mapped_ensembl_gene_ids[!duplicated(mapped_ensembl_gene_ids)];

				#There might still be duplicate uniprot gene symbols in names(mapped_ensembl_gene_ids)
				mapped_ensembl_gene_ids <- mapped_ensembl_gene_ids[!duplicated(names(mapped_ensembl_gene_ids))];

				info(file_logger, paste("The 1-1 uniprot_gn_symbol-ensembl_gene_id mapping has ", length(mapped_ensembl_gene_ids), " entries", sep = ""));
			}else
			{
				save(list=ls(), file = paste(output_dir, "/", "AnnotationOpenTargets.RData", sep = ""));
				error(file_logger, "Unsupported gene identifier encountered! Currently only supported gene identifiers are: ensembl_gene_id, entrezgene_id and uniprot_gn_symbol. Therefore, program will terminate here.");
				base::sink();
				cat("Unsupported gene identifier encountered! Currently only supported gene identifiers are: ensembl_gene_id, entrezgene_id and uniprot_gn_symbol. Therefore, program will terminate here.", "\n", sep = "");
				cat("The program: Annotation_using_OpenTargets.R failed at 'Prepare annotation data using biomaRt'. Check log for more details: ", log_file, "\n", sep = "");
				q(save = "no", status = 1);
			}

			## Update input genelist --------------------
      		genelist <- mapped_ensembl_gene_ids[genelist];

			## Remove NA
			genelist <- genelist[!is.na(names(genelist))];

			## Remove genes entries ""
			genelist <- genelist[genelist != ""];

			info(file_logger, paste("Updated query genelist has number of genes: ", length(genelist), sep = ""));
		}
	},
	warning = function(w)
	{
		warn(file_logger, "Caught a warning while preparing annotation data for the genes. Something unexpected must have happened.");
		print(w);
		base::sink();
		print(w);
		cat("The program: Annotation_using_OpenTargets.R failed at 'Check gene identifier type, and perform conversion to Ensembl gene id, if required'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	},
	error = function(e)
	{
		error(file_logger, "An error occured while preparing annotation data for the genes.");
		print(e);
		base::sink();
		print(e);
		cat("The program: Annotation_using_OpenTargets.R failed at 'Check gene identifier type, and perform conversion to Ensembl gene id, if required'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	}
)

save(list=ls(), file = paste(output_dir, "/", "AnnotationOpenTargets.RData", sep = ""));


## =====================================================================================
## Annotation using 'Open Targets' database
## =====================================================================================
tryCatch(
	expr =
	{
		#annotation_using_OpenTargets_df <- annotateUsingOpenTargetsViaMultipleAPIcalls(ensembl_gene_ids = genelist, file_logger);
		annotation_using_OpenTargets_df <- annotateUsingOpenTargetsViaSingleAPIcall(ensembl_gene_ids = genelist, file_logger)

		write.table(annotation_using_OpenTargets_df, file = paste(output_dir, "/", output_prefix, "annotation_with_OpenTargets.tsv", sep = ""), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE);
	},
	warning = function(w)
	{
		warn(file_logger, "Caught a warning while annotating genes with 'Open Targets'. Something unexpected must have happened. See below the warning message: ");
		print(w);
		base::sink();
		print(w);
		cat("The program: Annotation_using_OpenTargets.R failed at 'Annotation using 'Open Targets' database'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	},
	error = function(e)
	{
		error(file_logger, "An error occured while annotating genes with 'Open Targets'. See below the error message: ");
		print(e);
		base::sink();
		print(e);
		cat("The program: Annotation_using_OpenTargets.R failed at 'Annotation using 'Open Targets' database'. Check log for more details: ", log_file, "\n", sep = "");
		q(save = "no", status = 1);
	}
)


save(list=ls(), file = paste(output_dir, "/", "AnnotationOpenTargets.RData", sep = ""));


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
cat("The program: Annotation_using_OpenTargets.R executed successfully!", "\n", sep = "");
