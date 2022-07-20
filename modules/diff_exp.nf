process DIFF_EXP {
	container 'leenab/diff-expression:04'
	memory '32 GB'
	echo true
	publishDir "${params.outdir}", mode:'copy', overwrite: true

	input:
		path 'count.csv'
		path 'info.csv'
	output:
		path ("${params.output_dir}/*"), emit: diff_exp_all_data
		path ("${params.output_dir}/**/*.csv"), emit: diff_exp_data
		path ("${params.output_dir}/${params.output_prefix}_counts_after_filter.csv"), emit: ref_gene_file

	script:
	def usecpmcutoff = params.use_CPM_cutoff ? ' --use_CPM_cutoff TRUE': ' --use_CPM_cutoff FALSE'

    def otherknownfactorcols
    if(params.other_known_factors_columns == null) {
        otherknownfactorcols = ''
    } else {
        otherknownfactorcols = " --other_known_factors_columns $params.other_known_factors_columns"
    }

    def fulldesign
    if(params.full_design == null) {
        fulldesign = ''
    }
    else {
        fulldesign =  " --full_design $params.full_design"
    }

    def reduceddesign
    if (params.reduced_design == null) {
        reduceddesign = ''
    } else {
        reduceddesign = " --reduced_design $params.reduced_design"
    }

    def topnforviz
    if (params.top_n_for_viz == null) {
        topnforviz = ''
    } else {
        topnforviz = " --top_n_for_viz $params.top_n_for_viz"
    }

	"""
		Differential_gene_expression_analysis_on_counts.R --DGE_tool $params.DGE_tool \\
		--count_file count.csv --gene_id_column $params.gene_id_column --count_column $params.count_column \\
		--counts_threshold $params.counts_threshold \\
		--sample_coverage_threshold $params.sample_coverage_threshold  --normalization_method $params.normalization_method \\
		--pseudocounts $params.pseudocounts --sample_info_file info.csv --sample_id_column $params.sample_id_column \\
		--primary_factor_column $params.primary_factor_column --contrasts $params.contrasts \\
		--DESeq2_dispersion_method $params.DESeq2_dispersion_method \\
		--DESeq2_test_type $params.DESeq2_test_type --stat_sig_measure $params.stat_sig_measure \\
		--stat_cutoff $params.stat_cutoff --log2_fold_change_cutoff $params.log2_fold_change_cutoff \\
		--output_dir $params.output_dir --output_prefix $params.output_prefix --log_file $params.log_file \\
		$usecpmcutoff \\
		$otherknownfactorcols \\
		$fulldesign \\
		$reduceddesign \\
		$topnforviz
	"""
}