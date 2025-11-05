nextflow.enable.dsl = 2

// ───────────── PARAMETERS ─────────────
params.input_dir   = "${projectDir}/pp_data/"
params.output_dir  = "${projectDir}/results"
params.tax_train   = "${projectDir}/silva_nr_v138_train_set.fa.gz"
params.tax_species = "${projectDir}/silva_species_assignment_v138.fa.gz"

params.preprocess_r = "${projectDir}/poopie_pipeline.R"
params.summary_r    = "${projectDir}/generate_summary_single.R"
params.biomarker_r  = "${projectDir}/biomarker_single.R"
params.report_py    = "${projectDir}/poopie_report.py"
params.know_b       = "${projectDir}/pp_report.json"

// ───────────── PROCESS: PREPROCESS ─────────────
process PREPROCESS {
    tag "$sample_id"
    publishDir "${params.output_dir}/preprocess", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "results/rds/ps_rel.rds", emit: ps_rds
    path "results/json/sample_summary.json", emit: json_out

    script:
    """
    echo "[INFO] Running preprocessing for ${sample_id}"
    echo "Forward: ${reads[0]}"
    echo "Reverse: ${reads[1]}"

    mkdir -p results/rds results/json

    Rscript ${params.preprocess_r} \
      --read1 ${reads[0]} \
      --read2 ${reads[1]} \
      --output results \
      --taxonomy_train ${params.tax_train} \
      --taxonomy_species ${params.tax_species} \
      --threads 4
    """
}

// ───────────── PROCESS: SUMMARY ─────────────
process SUMMARY {
    tag "$sample_id"
    publishDir "${params.output_dir}/summary", mode: 'copy'

    input:
    path ps_rds

    output:
    path "results/reports/genus_abundance.csv", emit: genus_csv

    script:
    """
    echo "[INFO] Generating genus summary..."
    mkdir -p results/reports

    Rscript ${params.summary_r} \
      --phyloseq_rds ${ps_rds} \
      --output_dir results/reports
    """
}

// ───────────── PROCESS: BIOMARKERS ─────────────
process BIOMARKERS {
    tag "$sample_id"
    publishDir "${params.output_dir}/biomarkers", mode: 'copy'

    input:
    path ps_rds

    output:
    path "results/biomarkers_csv/"
    path "results/biomarkers_plots/"

    script:
    """
    echo "[INFO] Running biomarker discovery..."
    mkdir -p results/biomarkers_csv results/biomarkers_plots

    Rscript ${params.biomarker_r} \
      --rds ${ps_rds} \
      --sample ${sample_id} \
      --output results
    """
}

// ───────────── PROCESS: REPORT ─────────────
process REPORT {
    tag "$sample_id"
    publishDir "${params.output_dir}/pdf", mode: 'copy'

    input:
    path json_file
    path genus_csv

    output:
    path "*.pdf"
    path "*.txt"

    script:
    """
    echo "[INFO] Generating final PDF report for ${sample_id}"
    mkdir -p results/pdf

    python3 ${params.report_py} \
      --kb ${params.know_b} \
      --input_json ${json_file} \
      --genus_csv ${genus_csv} \
      --sample ${sample_id} \
      --output results/pdf/${sample_id}_Report.pdf
    """
}

// ───────────── WORKFLOW ─────────────
workflow {

    // Detect paired FASTQs
    reads = Channel
        .fromFilePairs("${params.input_dir}/*_R{1,2}_001.fastq.gz", flat: true)
        .map { id, files -> tuple(id, files) }
        .tap { id, files ->
            println "[FASTQ PAIR] ${id}"
            println "  R1: ${files[0]}"
            println "  R2: ${files[1]}"
        }

    // Run pipeline
    preprocess_ch = PREPROCESS(reads)
    summary_ch    = SUMMARY(preprocess_ch.ps_rds)
    biomarkers_ch = BIOMARKERS(preprocess_ch.ps_rds)
    REPORT(preprocess_ch.json_out, summary_ch.genus_csv)
}

