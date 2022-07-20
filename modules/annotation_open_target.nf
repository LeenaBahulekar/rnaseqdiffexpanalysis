process ANNOTATION_OPEN_TARGET {
	container 'leenab/diff-expression:04'
	memory '32 GB'
	echo true
	publishDir "${params.result_dir}/${diff_exp_dir}/annotation_open_target", mode:'copy', overwrite: true

	input:
		tuple val(diff_exp_dir), path ('diff_exp_cnt.csv')

	output:
		path ("results/*"), emit: annotation_op

	script:
	"""
		Annotation_using_OpenTargets.R --genelist_file diff_exp_cnt.csv \\
		--gene_id_column $params.annotation_gene_id_column --output_dir $params.annotation_result_dir \\
		--log_file $params.annotation_log_file --output_prefix $params.annotation_output_prefix
	"""
}