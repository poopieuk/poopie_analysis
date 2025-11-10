nextflow.enable.dsl=2

// ===============================
// PARAMETERS
// ===============================
params.input_dir   = "s3://poopie-data/"          // Input FASTQs from S3
params.output_dir  = "${projectDir}/results"
params.tax_train   = "$s3://poopie-data/silva_nr_v138_train_set.fa.gz"
params.tax_species = "$s3://poopie-data/silva_species_assignment_v138.fa.gz"
params.preprocess_r = "${projectDir}/poopie-pipeline/poopie_pipeline.R"
params.summary_r    = "${projectDir}/poopie-pipeline/generate_summary_single.R"
params.biomarker_r  = "${projectDir}/poopie-pipeline/biomarker_single.R"
params.report_py    = "${projectDir}/poopie-pipeline/poopie_report.py"
params.kb_json      = "${projectDir}/poopie-pipeline/pp_report.json"

params.sample_id    = "MS205-N715-A-S505-A_S92_L001"

// Supabase upload configuration
params.supabase_url     = "${SUPABASE_URL ?: 'https://tbyenonhykkizfdbcpnz.supabase.co'}"
params.supabase_key     = "${SUPABASE_KEY ?: 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InRieWVub25oeWtraXpmZGJjcG56Iiwicm9sZSI6ImFub24iLCJpYXQiOjE3NTc0MjY5ODEsImV4cCI6MjA3MzAwMjk4MX0.XbS2XgZTYDjoa6SrY4QrwMBVXxW315lYG2AKe4sheOU'}"
params.supabase_bucket  = "${SUPABASE_BUCKET ?: 'poopie-reports'}"

// ===============================
// PROCESS: PREPROCESS
// ===============================
process PREPROCESS {
    tag "$sample_id"
    publishDir "${params.output_dir}/preprocess", mode: 'copy'

    input:
    tuple val(sample_id), path(dummy)

    output:
    tuple val(sample_id), path("rds/ps_rel.rds"), emit: ps_rds
    tuple val(sample_id), path("json/full_microbiome_summary.json"), emit: json_out

    script:
    """
    echo "[INFO] Running preprocessing for ${sample_id}"
    mkdir -p results/rds results/json

    Rscript ${params.preprocess_r} \\
        --input ${params.input_dir} \\
        --output . \\
        --sample_id ${sample_id} \\
        --taxonomy_train ${params.tax_train} \\
        --taxonomy_species ${params.tax_species} \\
        --threads 4
    """
}

// ===============================
// PROCESS: SUMMARY
// ===============================
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

// ===============================
// PROCESS: BIOMARKERS
// ===============================
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

// ===============================
// PROCESS: REPORT
// ===============================
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

// ===============================
// PROCESS: UPLOAD_SUPABASE
// ===============================
process UPLOAD_SUPABASE {
    tag "$sample_id"
    publishDir "${params.output_dir}/upload_logs", mode: 'copy'

    input:
    tuple val(sample_id), path(pdf_file), path(json_file)

    env SUPABASE_URL = params.supabase_url
    env SUPABASE_KEY = params.supabase_key
    env SUPABASE_BUCKET = params.supabase_bucket

    output:
    path "upload_done_${sample_id}.txt"

    script:
    """
    echo "[INFO] Uploading final outputs for ${sample_id} to Supabase Storage..."

    for FILE in ${pdf_file} ${json_file}; do
      BASENAME=\$(basename "\${FILE}")
      MIME_TYPE=\$(file --mime-type -b "\${FILE}")
      echo "[UPLOAD] -> \${SUPABASE_BUCKET}/\${BASENAME} (\${MIME_TYPE})"

      curl -s -X POST "\${SUPABASE_URL}/storage/v1/object/\${SUPABASE_BUCKET}/\${BASENAME}" \\
        -H "Authorization: Bearer \${SUPABASE_KEY}" \\
        -H "Content-Type: \${MIME_TYPE}" \\
        --data-binary "@\${FILE}" || echo "[WARN] Upload failed for \${BASENAME}"
    done

    echo "✅ Uploaded PDF and JSON for ${sample_id} to Supabase bucket: \${SUPABASE_BUCKET}" > upload_done_${sample_id}.txt
    """
}

// ===============================
// WORKFLOW
// ===============================
workflow {
    single_sample_ch = Channel.of(tuple(params.sample_id, file('dummy.txt')))

    preprocess_ch = PREPROCESS(single_sample_ch)
    summary_ch    = SUMMARY(preprocess_ch.ps_rds)
    biomarker_ch  = BIOMARKERS(preprocess_ch.ps_rds)
    report_ch     = REPORT(preprocess_ch.json_out)

    // Combine PDF + JSON → Supabase upload
    upload_input = report_ch.report_pdfs
        .combine(report_ch.report_jsons)
        .map { pdf, json -> tuple(params
   
    UPLOAD_SUPABASE(upload_input)
}
