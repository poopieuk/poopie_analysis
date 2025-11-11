nextflow.enable.dsl=2

// ===============================
// PARAMETERS
// ===============================
params.input_dir   = "s3://poopie-data/inputs"
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
params.supabase_key     = System.getenv('SUPABASE_KEY') ?: 'YOUR_SUPABASE_KEY'
params.supabase_bucket  = System.getenv('SUPABASE_BUCKET') ?: 'reports'

// ===============================
// PROCESS DEFINITIONS
// ===============================
process PREPROCESS {
    tag "$sample_id"
    publishDir "${params.output_dir}/preprocess", mode: 'copy'

    input:
    val(sample_id)
    path input_file

    output:
    tuple val(sample_id), path("rds/ps_rel.rds"), emit: ps_rds
    tuple val(sample_id), path("json/full_microbiome_summary.json"), emit: json_out

    script:
    """
    echo "[INFO] Running preprocessing for ${sample_id}"
    echo "File staged: ${input_file}"
    mkdir -p results/rds results/json

    Rscript ${params.preprocess_r} \\
        --input ${input_file} \\
        --output . \\
        --sample_id ${sample_id} \\
        --taxonomy_train ${params.tax_train} \\
        --taxonomy_species ${params.tax_species} \\
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

    // List input FASTQ files from S3
    input_files_ch = Channel.fromPath("${params.input_dir}/*.fastq.gz", checkIfExists: true)

    // Create (sample_id, file) tuples
    preprocess_in = input_files_ch.map { file -> tuple(params.sample_id, file) }

    // ✅ Correct call
    preprocess_ch = PREPROCESS(preprocess_in)

    // Downstream
    summary_ch    = SUMMARY(preprocess_ch.ps_rds)
    biomarker_ch  = BIOMARKERS(preprocess_ch.ps_rds)
    report_ch     = REPORT(preprocess_ch.json_out)

    upload_input = report_ch.report_pdfs
        .combine(report_ch.report_jsons)
        .map { pdf, json -> tuple(params.sample_id, pdf, json) }

    UPLOAD_SUPABASE(upload_input)
}
