nextflow.enable.dsl=2

// ===============================
// PARAMETERS
// ===============================
params.input_dir   = "s3://poopie-data/"
params.output_dir  = "${projectDir}/results"
params.tax_train   = "s3://poopie-data/silva_nr_v138_train_set.fa.gz"
params.tax_species = "s3://poopie-data/silva_species_assignment_v138.fa.gz"
params.preprocess_r = "${projectDir}/poopie_pipeline.R"
params.summary_r    = "${projectDir}/generate_summary_single.R"
params.biomarker_r  = "${projectDir}/biomarker_single.R"
params.report_py    = "${projectDir}/poopie_report.py"
params.kb_json      = "${projectDir}/pp_report.json"

params.sample_id    = "MS205-N715-A-S505-A_S92_L001"

params.supabase_url     = System.getenv('SUPABASE_URL') ?: 'https://tbyenonhykkizfdbcpnz.supabase.co'
params.supabase_key     = System.getenv('SUPABASE_KEY') ?: 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InRieWVub25oeWtraXpmZGJjcG56Iiwicm9sZSI6ImFub24iLCJpYXQiOjE3NTc0MjY5ODEsImV4cCI6MjA3MzAwMjk4MX0.XbS2XgZTYDjoa6SrY4QrwMBVXxW315lYG2AKe4sheOU'
params.supabase_bucket  = System.getenv('SUPABASE_BUCKET') ?: 'reports'


// ===============================
// PROCESS DEFINITIONS
// ===============================
process PREPROCESS {
    tag "$sample_id"
    publishDir "${params.output_dir}/preprocess", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path tax_train
    path tax_species

    output:
    tuple val(sample_id), path("rds/ps_rel.rds"), emit: ps_rds
    tuple val(sample_id), path("json/full_microbiome_summary.json"), emit: json_out

    script:
    """
    echo "[INFO] Running preprocessing for ${sample_id}"
    echo "[DEBUG] FASTQ files staged:"
    ls -lh *.fastq.gz

    mkdir -p results/rds results/json

    Rscript ${params.preprocess_r} \\
        --input . \\
        --output . \\
        --sample_id ${sample_id} \\
        --taxonomy_train ${tax_train} \\
        --taxonomy_species ${tax_species} \\
        --threads 4
    """
}





process SUMMARY {
    tag "$sample_id"
    publishDir "${params.output_dir}/summary", mode: 'copy'

    input:
    tuple val(sample_id), path(ps_rds)

    output:
    tuple val(sample_id), path("results/reports/genus_abundance.csv"), emit: genus_csv

    script:
    """
    mkdir -p results/reports
    Rscript ${params.summary_r} \\
        --phyloseq_rds ${ps_rds} \\
        --sample_id ${sample_id} \\
        --output_dir results/reports
    """
}

process BIOMARKERS {
    tag "$sample_id"
    publishDir "${params.output_dir}/biomarkers", mode: 'copy'

    input:
    tuple val(sample_id), path(ps_rds)

    output:
    tuple val(sample_id), path("results/biomarkers_csv/*.csv"), emit: biomarker_csv
    tuple val(sample_id), path("results/biomarkers_plots/*.png"), emit: biomarker_plots

    script:
    """
    mkdir -p results/biomarkers_csv results/biomarkers_plots
    echo "[INFO] Running biomarker discovery for ${sample_id}"

    CLEAN_ID=\$(echo ${sample_id} | sed 's/_R[12]_001//')
    Rscript ${params.biomarker_r} \\
        --rds ${ps_rds} \\
        --sample \${CLEAN_ID} \\
        --output results
    """
}

process REPORT {
    tag "$sample_id"
    publishDir "${params.output_dir}/pdf", mode: 'copy'

    input:
    tuple val(sample_id), path(json_file)

    output:
    tuple val(sample_id), path("results/pdf/*.pdf"), emit: report_pdfs
    tuple val(sample_id), path("results/pdf/*.json"), emit: report_jsons

    script:
    """
    mkdir -p results/pdf
    echo "[INFO] Generating timestamped report for ${sample_id}"

    python3 ${params.report_py} \\
        --kb ${params.kb_json} \\
        --input_json ${json_file} \\
        --sample ${sample_id} \\
        --output results/pdf
    """
}

process UPLOAD_SUPABASE {
    tag "$sample_id"
    publishDir "${params.output_dir}/supabase_upload", mode: 'copy'

    input:
    tuple val(sample_id), path(pdf), path(json)

    output:
    path "upload_done.txt", emit: upload_flag

    environment = [
        'SUPABASE_URL': params.supabase_url,
        'SUPABASE_KEY': params.supabase_key,
        'SUPABASE_BUCKET': params.supabase_bucket
    ]

    script:
    """
    echo "[INFO] Uploading files for ${sample_id} to Supabase..."
    python3 - <<PY
import os
from supabase import create_client, Client
url  = os.environ["SUPABASE_URL"]
key  = os.environ["SUPABASE_KEY"]
bucket = os.environ["SUPABASE_BUCKET"]
supabase: Client = create_client(url, key)
for f in ["${pdf}", "${json}"]:
    dest_name = os.path.basename(f)
    with open(f, "rb") as file_data:
        supabase.storage.from_(bucket).upload(dest_name, file_data)
print("✅ Upload complete for ${sample_id}")
PY
    echo "done" > upload_done.txt
    """
}


// ===============================
// WORKFLOW
// ===============================
workflow {

    // Define taxonomy databases
    tax_train_ch   = Channel.fromPath(params.tax_train, checkIfExists: true)
    tax_species_ch = Channel.fromPath(params.tax_species, checkIfExists: true)

    // Detect paired FASTQ files
    paired_fastqs_ch = Channel
        .fromPath("${params.input_dir}/*_{R1,R2}_001.fastq.gz", checkIfExists: true)
        .map { file ->
            def sid = file.name.replaceAll(/_R[12]_001\\.fastq\\.gz$/, '')
            tuple(sid, file)
        }
        .groupTuple()

    paired_fastqs_ch.view { "DEBUG: Paired FASTQs -> ${it}" }

    // Run processes
    preprocess_ch = PREPROCESS(paired_fastqs_ch, tax_train_ch, tax_species_ch)
    summary_ch    = SUMMARY(preprocess_ch.ps_rds)
    biomarker_ch  = BIOMARKERS(preprocess_ch.ps_rds)
    report_ch     = REPORT(preprocess_ch.json_out)

    upload_input = report_ch.report_pdfs
        .combine(report_ch.report_jsons)
        .map { pdf, json -> tuple(params.sample_id, pdf, json) }

    UPLOAD_SUPABASE(upload_input)
}
