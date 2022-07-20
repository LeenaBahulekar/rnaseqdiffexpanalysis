process PATHWAY_ENRICHMENT {
	container 'leenab/diff-expression:04'
	memory '32 GB'
	echo true
	publishDir "${params.outdir}/${diff_exp_dir}/pathway_enrichment", mode:'copy', overwrite: true
	
	input:
		tuple val(diff_exp_dir), path ('diff_exp_cnt.csv')
		path ref_gene_file
	
	output:
		path ("results/*"), emit: pathway_op 
		
	script:
	def runSPIA = params.run_SPIA ? ' --run_SPIA TRUE' : ' --run_SPIA FALSE'
	"""
		Functional_enrichment.R --genelist_file diff_exp_cnt.csv \\
		--gene_id_column $params.pathway_gene_id_column --log2fc_column $params.pathway_log2fc_column \\
		--reference_genelist_file  $ref_gene_file --ref_gene_id_column $params.ref_gene_id_column \\
		--species "$params.species" --data_sources $params.data_sources \\
		--threshold $params.threshold --output_dir $params.pathway_result_dir \\
		--output_prefix $params.pathway_output_prefix --log_file $params.pathway_log_file \\
		$runSPIA
	"""

}