nextflow.enable.dsl = 2

// ===============================
// PARAMETERS
// ===============================
params.input_dir   = "s3://source-bioscience-data/raw-data/22250580300021_1"
params.output_dir  = "${projectDir}/results"
params.tax_train   = "s3://poopie-data/silva_nr_v138_train_set.fa.gz"
params.tax_species = "s3://poopie-data/silva_species_assignment_v138.fa.gz"

params.preprocess_r = "${projectDir}/poopie_pipeline.R"
params.summary_r    = "${projectDir}/generate_summary_single.R"
params.biomarker_r  = "${projectDir}/biomarker_single.R"
params.report_py    = "${projectDir}/poopie_report.py"
params.kb_json      = "${projectDir}/pp_report.json"



params.supabase_url     = System.getenv('SUPABASE_URL') ?: 'https://tbyenonhykkizfdbcpnz.supabase.co'
params.supabase_key     = System.getenv('SUPABASE_KEY') ?: 'eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InRieWVub25oeWtraXpmZGJjcG56Iiwicm9sZSI6ImFub24iLCJpYXQiOjE3NTc0MjY5ODEsImV4cCI6MjA3MzAwMjk4MX0.XbS2XgZTYDjoa6SrY4QrwMBVXxW315lYG2AKe4sheOU'
params.supabase_bucket  = System.getenv('SUPABASE_BUCKET') ?: 'reports'

// ===============================
// PROCESS DEFINITIONS
// ===============================

// ---------- PREPROCESS ----------
process PREPROCESS {
    tag "$sample_id"
    publishDir "${params.output_dir}/preprocess/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path tax_train
    path tax_species
    path preprocess_r

    output:
    tuple val(sample_id), path("rds/ps_rel.rds"), emit: ps_rds
    tuple val(sample_id), path("json/full_microbiome_summary.json"), emit: json_out

    script:
    """
    echo "[INFO] Running preprocessing for ${sample_id}"
    echo "[DEBUG] FASTQ files staged:"
    ls -lh *.fastq.gz

    mkdir -p results/rds results/json

    CLEAN_ID=\$(basename \$(ls *.fastq.gz | head -n1) | sed 's/_R[12]\\.fastq\\.gz//')
    echo "[DEBUG] Clean sample ID -> \${CLEAN_ID}"
Rscript -e "install.packages('optparse', repos='https://cloud.r-project.org')"
    Rscript ${preprocess_r} \\
        --input . \\
        --output . \\
        --taxonomy_train ${tax_train} \\
        --taxonomy_species ${tax_species} \\
        --threads 4
    """

}

// ---------- SUMMARY ----------
process SUMMARY {
    tag "$sample_id"
    publishDir "${params.output_dir}/summary", mode: 'copy'

    input:
    tuple val(sample_id), path(ps_rds)
    path summary_r

    output:
    tuple val(sample_id), path("results/reports/*.csv"), emit: genus_csv

    script:
    """
    mkdir -p results/reports
    Rscript ${summary_r} \\
        --phyloseq_rds ${ps_rds} \\
        --sample_id ${sample_id} \\
        --output_dir results/reports
    """
}

// ---------- BIOMARKERS ----------
process BIOMARKERS {
    tag "$sample_id"
    publishDir "${params.output_dir}/biomarkers", mode: 'copy'

    input:
    tuple val(sample_id), path(ps_rds)
    path biomarker_r

    output:
    tuple val(sample_id), path("results/biomarkers_csv/*.csv"), emit: biomarker_csv
    tuple val(sample_id), path("results/biomarkers_plots/*.png"), emit: biomarker_plots

    script:
    """
    mkdir -p results/biomarkers_csv results/biomarkers_plots
    echo "[INFO] Running biomarker discovery for ${sample_id}"

    CLEAN_ID=\$(echo ${sample_id} | sed 's/_R[12]//')
    Rscript ${biomarker_r} \\
        --rds ${ps_rds} \\
        --sample \${CLEAN_ID} \\
        --output results
    """
}

// ---------- REPORT ----------
process REPORT {
    tag "$sample_id"
    publishDir "${params.output_dir}/pdf", mode: 'copy'

    input:
    tuple val(sample_id), path(json_file)
    path report_py
    path kb_json

    output:
    tuple val(sample_id), path("results/*.pdf"), emit: report_pdfs
    tuple val(sample_id), path("results/*_${sample_id}.txt"), emit: report_txts


    script:
    """
    mkdir -p results

    echo "[INFO] Generating timestamped report for ${sample_id}"

    python3 ${report_py} \
        --kb ${kb_json} \
        --input_json ${json_file} \
        --sample ${sample_id} \
        --output results/${sample_id}.pdf

    # Required for join()
    echo "Summary for ${sample_id}" > results/${sample_id}.txt
    """
}


// ---------- UPLOAD_SUPABASE ----------
process UPLOAD_SUPABASE {
    tag "$sample_id"
    publishDir "${params.output_dir}/supabase_upload", mode: 'copy'

    input:
    tuple val(sample_id), path(pdf), path(txt)

    output:
    path "upload_done.txt", emit: upload_flag

    script:
    """
    echo "[INFO] Uploading files for ${sample_id} to Supabase..."
    pip install --quiet supabase
    export SUPABASE_URL="${params.supabase_url}"
    export SUPABASE_KEY="${params.supabase_key}"
    export SUPABASE_BUCKET="${params.supabase_bucket}"

    python3 - << 'PY'
import os, mimetypes
from supabase import create_client, Client
from storage3.exceptions import StorageApiError
url=os.getenv("SUPABASE_URL"); key=os.getenv("SUPABASE_KEY"); bucket=os.getenv("SUPABASE_BUCKET")
supabase: Client = create_client(url, key)
files = [r"${pdf}", r"${txt}"]
for f in files:
    if not os.path.exists(f):
        print(f"⚠️ Skipping missing file: {f}"); continue
    dest=os.path.basename(f)
    mime=mimetypes.guess_type(f)[0] or "application/octet-stream"
    print(f"⬆️ Uploading {dest} ({mime})")
    try:
        with open(f,"rb") as fd:
            supabase.storage.from_(bucket).upload(dest, fd, file_options={"content-type": mime})
        print(f"✅ Uploaded {dest}")
    except StorageApiError as e:
        print(f"⚠️ Upload failed for {dest}: {e}")
print("🎉 All uploads completed!")
PY
    echo "done" > upload_done.txt
    """
}

// ===============================
// WORKFLOW
// ===============================
workflow {

    // --- static reference files as channels ---
    preprocess_r_ch = Channel.fromPath(params.preprocess_r, checkIfExists: true)
    summary_r_ch    = Channel.fromPath(params.summary_r, checkIfExists: true)
    biomarker_r_ch  = Channel.fromPath(params.biomarker_r, checkIfExists: true)
    report_py_ch    = Channel.fromPath(params.report_py, checkIfExists: true)
    kb_json_ch      = Channel.fromPath(params.kb_json, checkIfExists: true)

    tax_train_ch    = Channel.fromPath(params.tax_train, checkIfExists: true)
    tax_species_ch  = Channel.fromPath(params.tax_species, checkIfExists: true)

    // --- group paired FASTQs ---
    paired_fastqs_ch = Channel.fromFilePairs(
    "${params.input_dir}/*_R{1,2}_001.fastq.gz"
    )

    // 🔍 Automatically find R1/R2 FASTQs inside subfolders
paired_fastqs_ch = Channel
    .fromFilePairs("${params.input_dir}/*_{R1,R2}.fastq.gz", flat: true, checkIfExists: true)
    .map {sample_id, reads ->
        tuple(sample_id, reads)
    }

paired_fastqs_ch.view { "DEBUG: Paired FASTQs -> ${it}" }


    paired_fastqs_ch.view { "DEBUG: Paired FASTQs -> ${it}" }

    // --- run the workflow chain ---
    preprocess_ch = PREPROCESS(paired_fastqs_ch, tax_train_ch, tax_species_ch, preprocess_r_ch)
    summary_ch    = SUMMARY(preprocess_ch.ps_rds, summary_r_ch)
    biomarker_ch  = BIOMARKERS(preprocess_ch.ps_rds, biomarker_r_ch)
    report_ch     = REPORT(preprocess_ch.json_out, report_py_ch, kb_json_ch)

    // --- combine PDF and TXT outputs for upload --- 
     upload_input_ch = report_ch.report_pdfs .join(report_ch.report_txts) .map { sample_id, pdf, txt -> [sample_id, pdf, txt] }


    upload_input_ch.view { "DEBUG Upload tuple -> ${it}" }

    upload_input_ch | UPLOAD_SUPABASE
}
