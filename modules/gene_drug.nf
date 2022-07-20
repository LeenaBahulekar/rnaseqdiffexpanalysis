process GENE_DRUG {
	container 'leenab/diff-expression:04'
	memory '32 GB'
	echo true
	publishDir "${params.outdir}/${diff_exp_dir}/gene_drug", mode:'copy', overwrite: true
	
	input:
		tuple val(diff_exp_dir), path ('diff_exp_cnt.csv')
		path ppi
		path adjacency
	
	output:
		path ("${params.gene_result_dir}/*"), emit: gene_drug_op 
		
	script:
	"""
		gene_drug_association.R --genelist_file diff_exp_cnt.csv --ppi_file $ppi \\
		--drug_target_file $adjacency --gene_id_column $params.genedrug_gene_id_column \\
		--log2fc_column $params.log2fc_column --output_dir $params.gene_result_dir \\
		--log_file $params.gene_drug_log_file
	"""
}
