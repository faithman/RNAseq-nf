#!/usr/bin/env nextflow

params.fqs = "$baseDir/test.tsv"
params.transcriptome = "$baseDir/test_data/c.elegans.cdna.ncrna.fa"
params.multiqc = "$baseDir/multiqc"
params.fragment_len = '250'
params.fragment_sd = '50'
params.bootstrap = '100'
params.experiment = "$baseDir/experiment_info.txt"
params.email = ""

File fq_file = new File(params.fqs)

log.info """\
         R N A S E Q - N F   P I P E L I N E  
         ===================================
         transcriptome: ${ params.transcriptome }
         fqs          : ${ params.fqs }
         output       : ${ params.output }
         fragment_len : ${ params.fragment_len } 
         fragment_sd  : ${ params.fragment_sd }
         bootstrap    : ${ params.bootstrap }
         experiment   : ${ params.experiment }
         email        : ${ params.email }

         """
         .stripIndent()


transcriptome_file = file(params.transcriptome)
multiqc_file = file(params.multiqc)
exp_file = file(params.experiment)
/*
 * Make sure files exist
 */

if( !transcriptome_file.exists() ) exit 1, "Missing transcriptome file: ${transcriptome_file}"

if( !exp_file.exists() ) exit 1, "Missing Experiment parameters file: ${exp_file}"

Channel
    .fromPath(params.fqs, checkIfExists: true)
    .ifEmpty {exit 1, "sample sheet not found"}
    .splitCsv(header:true, sep: "\t")
    .map { row -> row.fq1 = row.fq1;row}
    .map { row -> row.fq2 = row.fq2;row}
    .map { row -> [row.sample, file(row.fq1), file(row.fq2)]}
    .view()
    .into {read_1_ch; read_3_ch}


process kal_index {

    input:
        file transcriptome_file

    output:
        file "transcriptome.index" into transcriptome_index

    script:
        //
        // Kallisto mapper index
        //
        """
        kallisto index -i transcriptome.index ${transcriptome_file}
        """
}

process kal_mapping {

    publishDir "${params.output}", mode: 'copy'

    tag "reads: $SM"

    input:
        file index from transcriptome_index
        tuple SM, path(fq1), path(fq2) from read_1_ch

    output:
        file "kallisto_${SM}" into kallisto_out_dirs

    script:
    //
    // Kallisto tools mapper
    //
        """
        mkdir kallisto_${SM}
        kallisto quant --bootstrap ${params.bootstrap} -i ${index} -o kallisto_${SM} ${fq1} ${fq2}
        """
}


process fastqc {

    tag "${ SM }"

    input:
	tuple SM, path(fq1), path(fq2) from read_3_ch

    output:
        file("${SM}_log") into fastqc_ch

    script:
        """
        mkdir -p ${SM}_log
        fastqc -o ${SM}_log -f fastq -q ${fq1} ${fq2}
        """
}


process sleuth {

	publishDir "${params.output}", mode: 'copy'

    input:
        file 'kallisto/*' from kallisto_out_dirs.collect()   
        file exp_file

    output: 
        file 'geneMode.so'
	file 'sleuth_significant.tsv'
	file 'gene_table.tsv'
	file 'full_sleuth_results.tsv'

    script:
        //
        // Setup sleuth R dependancies and environment
        //
     
        """
        sleuth.R kallisto ${exp_file}
        """
}

workflow.onComplete {
    summary = """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    """
    println summary
    def outlog = new File("${params.output}/log.txt")
    outlog.newWriter().withWriter {
        outlog << param_summary
        outlog << summary
    }
    // mail summary
    if (params.email) {
        ['mail', '-s', 'SEmRNA-seq-nf', params.email].execute() << summary
    }
}

