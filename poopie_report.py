"""
Your Poopie Report — Dynamic Microbiome Report Generator (v8 full)
--------------------------------------------------------------
✅ Full narrative version restored (all educational paragraphs)
✅ Improved pastel gauges for realistic, professional look
✅ PDF only (no Word)
"""

# ========= Imports =========

import shutil
import json
import pandas as pd
import argparse
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
import re

def normalize_id(s):
    """Normalize sample IDs for consistent matching between R output and Nextflow."""
    if s is None:
        return "unknown_sample"
    s = str(s).strip().lower()
    s = s.replace("-", "_").replace(" ", "_")
    s = re.sub(r"\.f(ast)?q(\.gz)?$", "", s)   # remove .fastq/.fq/.gz
    s = re.sub(r"_r[12]_001$", "", s)          # remove _R1_001/_R2_001
    s = re.sub(r"_s\d+_l\d+", "", s)           # remove _S81_L001, etc.
    s = re.sub(r"_l\d+$", "", s)               # remove trailing lane like _L001
    s = re.sub(r"_$", "", s)                   # cleanup trailing _
    return s





from fuzzywuzzy import fuzz
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, PageBreak, HRFlowable, Flowable
)
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import cm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus.tableofcontents import TableOfContents

import argparse

parser = argparse.ArgumentParser(description="Generate Poopie microbiome report")
parser.add_argument("--genus_csv", help="Path to genus_abundance.csv", required=False)

parser.add_argument("--input_json", help="Path to full_microbiome_summary_clean.json")
parser.add_argument("--kb", help="Knowledge base JSON")
parser.add_argument("--sections", help="Sections JSON")
parser.add_argument("--output", help="Path to output PDF")

parser.add_argument("--sample", help="Sample ID", default=None)
parser.add_argument("--all", action="store_true", help="Generate report for all samples")

args = parser.parse_args()

print(f"[INFO] Input JSON: {args.input_json}")
print(f"[INFO] Knowledge base: {args.kb}")
print(f"[INFO] Sections JSON: {args.sections}")

print(f"[INFO] Output file: {args.output}")
print(f"[INFO] Sample: {args.sample}")
print(f"[INFO] All samples: {args.all}")

# --- Load knowledge base and section JSONs safely ---
try:
    with open(args.kb, "r") as f:
        kb_data = json.load(f)
        kb = kb_data  # alias for backward compatibility

    kb_normalized = {k.lower().replace(" ", "_"): v for k, v in kb_data.items()}
    print(f"[INFO] Loaded {len(kb_normalized)} entries from KB.")
except Exception as e:
    print(f"⚠️ Warning: Could not load KB file ({args.kb}): {e}")
    kb_normalized = {}

try:
    with open(args.sections, "r") as f:
        section_texts = json.load(f)
    print(f"[INFO] Loaded sections JSON successfully.")
except Exception as e:
    print(f"⚠️ Warning: Could not load sections file ({args.sections}): {e}")
    section_texts = {"sections": {}}

# --- Default thresholds for abundance and health scoring ---
thresholds = {
    "low": 0.1,
    "medium": 1.0,
    "high": 5.0,
    "very_high": 10.0
}

# Optional: override thresholds if provided in KB JSON
if "thresholds" in kb_normalized:
    thresholds.update(kb_normalized["thresholds"])


class CleanGauge(Flowable):
    """Universal flat gauge for scores or abundance (Chuckling Goat–style)."""

    def __init__(self, value, width=450, height=12, max_value=10.0, label=None, mode="score"):
        super().__init__()
        self.value = max(0, min(value, max_value))
        self.width = width
        self.height = height
        self.max_value = max_value
        self.mode = mode  # "score" or "abundance"
        self.label = label or (f"{value:.1f}/10" if mode == "score" else f"{value:.2f}%")

    def _color_for_value(self, v):
        # Consistent red–orange–green mapping
        ratio = v / self.max_value
        if ratio < 0.4:
            return colors.HexColor("#E57373")  # red
        elif ratio < 0.7:
            return colors.HexColor("#FDD835")  # amber
        else:
            return colors.HexColor("#81C784")  # green

    def draw(self):
        c = self.canv
        fill_width = (self.value / self.max_value) * self.width

        # background
        c.setFillColor(colors.HexColor("#F2F2F2"))
        c.roundRect(0, 0, self.width, self.height, self.height / 2, stroke=0, fill=1)

        # fill
        c.setFillColor(self._color_for_value(self.value))
        c.roundRect(0, 0, fill_width, self.height, self.height / 2, stroke=0, fill=1)

        # subtle border
        c.setStrokeColor(colors.HexColor("#CCCCCC"))
        c.setLineWidth(0.3)
        c.roundRect(0, 0, self.width, self.height, self.height / 2, stroke=1, fill=0)

        # label
        c.setFont("Helvetica", 8)
        c.setFillColor(colors.HexColor("#333333"))
        c.drawString(self.width + 8, 1, self.label)


# ========= Dynamic Scoring Helper =========
def compute_dynamic_score(abundance_df, target_genera, ideal_total, invert=False):
    """
    Computes a biologically scaled 0–10 score.
    - ideal_total = % abundance considered "optimal" (e.g. 10 means 10%)
    - invert=True for pathogens (where less is better)
    """
    if abundance_df.empty:
        return 0.0

    total = abundance_df[abundance_df["Genus"].isin(target_genera)]["Abundance"].sum()
    ratio = total / ideal_total

    if invert:
        # For pathogens: higher abundance means lower score
        score = max(0.0, 10.0 - (ratio * 10))
    else:
        # For beneficials: higher abundance means higher score
        ratio = min(ratio, 1.0)
        score = math.pow(ratio, 0.5) * 10  # smoother scaling

    return round(score, 1)


# ========= Parser & Data Load =========
#parser = argparse.ArgumentParser(description="Generate personalized microbiome report (PDF only).")
#parser.add_argument("--sample", default="MS205-N705-A-S507-A_S109_L001_R1_001.fastq.gz")
#parser.add_argument("--all", action="store_true")
#args = parser.parse_args()

# --- Load data ---
#try:
 #   with open("/Users/ahamed/Downloads/results/json/full_microbiome_summary_clean.json", "r") as f:
 #       data_all = json.load(f)
#except FileNotFoundError:
 #   sys.exit("❌ Could not find /Users/ahamed/Downloads/results/json/full_microbiome_summary_clean.json")
# --- Load data from CLI arg ---
try:
    with open(args.input_json, "r") as f:
        data_all = json.load(f)
    print(f"✅ Loaded input JSON from {args.input_json}")
except Exception as e:
    sys.exit(f"❌ Failed to load input JSON ({args.input_json}): {e}")

# --- Normalize JSON structure early ---
if isinstance(data_all, list):
    data_all = {
        "abundance": data_all,
        "diversity": data_all,
        "metadata": data_all,
        "scores": {}
    }

# --- Helper function to ensure consistent Sample columns ---
def normalize_df(section):
    df = pd.DataFrame(section)
    if df.empty:
        return pd.DataFrame(columns=["Sample", "Sample_norm"])
    # Fallback if 'Sample' missing
    if "Sample" not in df.columns and "_row" in df.columns:
        df["Sample"] = df["_row"]
    if "Sample" not in df.columns:
        df["Sample"] = "Unknown"
    # Normalize
    df["Sample"] = df["Sample"].astype(str).str.strip()
    df["Sample_norm"] = df["Sample"].str.lower().str.replace(r'[^a-z0-9]+', '_', regex=True)
    return df

# --- Create base DataFrames safely ---
abundance_df = normalize_df(data_all.get("abundance", []))
diversity_df = normalize_df(data_all.get("diversity", []))
meta_df      = normalize_df(data_all.get("metadata", []))

# --- Normalize sample IDs inside all dataframes for matching ---
import re

def normalize_id(x):
    """Normalize sample IDs exactly like the JSON normalization."""
    if isinstance(x, str):
        x = x.lower().replace("-", "_").replace(".", "_")
        x = re.sub(r'[^a-z0-9]+', '_', x)
        return x.strip('_')
    return x

for df in [abundance_df, diversity_df, meta_df]:
    if "Sample" in df.columns:
        df["Sample_norm"] = df["Sample"].apply(normalize_id)


# --- Normalize Sample IDs in all dataframes ---
for df_name, df in [("abundance", abundance_df), ("diversity", diversity_df)]:
    if "Sample" in df.columns:
        df["Sample_norm"] = df["Sample"].apply(normalize_id)
        print(f"[DEBUG] Normalized {len(df)} rows in {df_name}_df")
    else:
        print(f"[WARN] No 'Sample' column found in {df_name}_df")

# --- Optional: Load genus_abundance.csv ---
genus_df = None
if args.genus_csv:
    try:
        genus_df = pd.read_csv(args.genus_csv, index_col=0)  # rows = Sample, cols = Genera
        genus_df.index.name = "Sample"
        print(f"[INFO] Loaded genus_abundance.csv with shape {genus_df.shape} from {args.genus_csv}")
        # Convert wide (samples x genera) -> long (Sample, Genus, Abundance)
        genus_long = (
            genus_df
            .reset_index()
            .melt(id_vars=["Sample"], var_name="Genus", value_name="Abundance")
            .fillna(0)
        )
        # Ensure types
        genus_long["Sample"] = genus_long["Sample"].astype(str)
        genus_long["Genus"] = genus_long["Genus"].astype(str)
        genus_long["Abundance"] = pd.to_numeric(genus_long["Abundance"], errors="coerce").fillna(0.0)

        # ✅ Use genus CSV as primary abundance source
        abundance_df = genus_long
        print("[INFO] Using genus_abundance.csv as primary abundance source.")
    except Exception as e:
        print(f"⚠️ Warning: Could not read genus_abundance.csv ({args.genus_csv}): {e}")
else:
    print("[INFO] No genus_abundance.csv provided — using abundance from JSON only.")

scores = data_all.get("scores", {})

print(f"[INFO] DataFrames -> abundance:{len(abundance_df)}, diversity:{len(diversity_df)}, metadata:{len(meta_df)}")


# ========= Normalize Genus Function =========
def normalize_genus(name: str) -> str:


    """Normalize genus names for consistent lookup."""

    if not isinstance(name, str):
        return ""
    return name.strip().lower().replace("/", "_").replace("-", "_").replace(" ", "_")



# ========= Load JSON Knowledgebase =========
# --- Load Knowledge Base JSON safely ---
try:
    with open(args.kb, "r") as f:
        kb_data = json.load(f)
    print(f"✅ Loaded knowledge base with {len(kb_data)} entries.")
except Exception as e:
    print(f"❌ Failed to load knowledge base from {args.kb}: {e}")
    sys.exit(1)
# --- Load Input JSON ---
try:
    with open(args.input_json, "r") as f:
        data_all = json.load(f)
    print(f"✅ Loaded input JSON with {len(data_all)} records.")
except Exception as e:
    print(f"❌ Failed to load input JSON: {e}")
    sys.exit(1)

# --- Load Sections JSON ---
sections = None
if getattr(args, "sections_json", None):
    try:
        with open(args.sections_json) as f:
            sections = json.load(f)
        print(f"✅ Loaded sections JSON: {args.sections_json}")
    except Exception as e:
        print(f"⚠️ Warning: Could not load sections file ({args.sections_json}): {e}")
else:
    print("⚠️ Warning: No sections JSON provided — continuing without it.")

#
#try:
#    with open("pp_report.json", "r") as f:
#        kb_raw = json.load(f)



#    kb = kb_raw.get("bacteria", {})
#    thresholds = kb_raw.get("thresholds", {"high": 10.0, "moderate": 1.0})

#    try:
#        with open("section_explanations.json", "r") as f:
#            section_texts = json.load(f)
#    except FileNotFoundError:
#        sys.exit("❌ Could not find section_explanations.json")

    # ✅ Now we can safely normalize
    kb_normalized = {normalize_genus(k): v for k, v in kb.items()}

# ========= Load Knowledge Base (pp_report.json) =========
try:
    kb_path = args.kb if hasattr(args, "kb") and args.kb else "pp_report.json"
    with open(kb_path, "r") as f:
        kb_raw = json.load(f)
    print(f"✅ Loaded knowledge base from {kb_path}")
except FileNotFoundError:
    sys.exit(f"❌ Could not find knowledge base JSON at {kb_path}")
except Exception as e:
    sys.exit(f"❌ Failed to load knowledge base: {e}")

# Extract bacteria info and thresholds from pp_report.json
if "bacteria" in kb_raw:
    kb = kb_raw["bacteria"]
else:
    # Fallback if file is already a simple genus:info dict
    kb = kb_raw
    print("[WARN] 'bacteria' key not found; using top-level JSON structure.")

thresholds = kb_raw.get("thresholds", {"high": 10.0, "moderate": 1.0})

# Normalize all genus keys for consistent lookup
kb_normalized = {normalize_genus(k): v for k, v in kb.items()}
print(f"[INFO] Loaded {len(kb_normalized)} normalized genus entries from KB.")



#try:
#    with open("section_explanations.json", "r") as f:
#        section_texts = json.load(f)
#except FileNotFoundError:
#        sys.exit("❌ Could not find section_explanations.json")

    # ✅ Now we can safely normalize
#kb_normalized = {normalize_genus(k): v for k, v in kb.items()}


#except FileNotFoundError:
  #  sys.exit("❌ Could not find pp_report.json")


#abundance_df = pd.DataFrame(data_all["abundance"])


#diversity_df = pd.DataFrame(data_all["diversity"])
# Handle JSON list or dict structures gracefully
if isinstance(data_all, list):
    df = pd.DataFrame(data_all)
    abundance_df = df
    diversity_df = df  # same data, you can later subset columns if needed
elif isinstance(data_all, dict):
    abundance_df = pd.DataFrame(data_all.get("abundance", []))
    diversity_df = pd.DataFrame(data_all.get("diversity", []))

    # Handle missing metadata gracefully
    if "metadata" in data_all and data_all["metadata"]:
        meta_df = pd.DataFrame(data_all["metadata"])
        print(f"✅ Loaded metadata with {len(meta_df)} entries.")
    else:
        print("⚠️ No metadata found in JSON — using empty dataframe.")
        meta_df = pd.DataFrame()

    #meta_df = pd.DataFrame(data_all["metadata"])
    meta_df = pd.DataFrame(data_all.get("metadata", []))
    if meta_df.empty:
        meta_df = pd.DataFrame(columns=["Sample", "Sample_norm"])
    else:
        meta_df["Sample_norm"] = meta_df["Sample"].apply(normalize_id)


else:
    raise ValueError("Unexpected JSON format — expected list or dict with 'abundance' and 'diversity' keys.")

#meta_df = pd.DataFrame(data_all["metadata"])
#scores = data_all.get("scores", {})
# --- Handle 'scores' safely depending on JSON structure ---
if isinstance(data_all, dict):
    scores = data_all.get("scores", {})
elif isinstance(data_all, list):
    # fallback: no scores section in flat list JSON
    scores = {}
else:
    scores = {}



# ========= Helpers =========
def classify(v):
    if v >= 7.5:
        return "Balanced"
    elif v >= 6.0:
        return "Moderate"
    else:
        return "Needs Support"


def divider():
    return HRFlowable(width="100%", color=colors.lightgrey, thickness=0.7, spaceAfter=8, spaceBefore=4)


def on_page(canvas_obj, doc):
    """Draw footer page numbers, skipping title + TOC pages."""
    page_num = canvas_obj.getPageNumber()

    # Skip numbering for title and TOC pages
    if page_num <= 2:
        return

    canvas_obj.setFont("Helvetica", 9)
    canvas_obj.setFillColor(colors.grey)
    canvas_obj.drawCentredString(A4[0] / 2, 1.2 * cm, f"Page {page_num - 2}")


# ========= JSON Knowledgebase Lookup =========
from fuzzywuzzy import process


def get_bacteria_info(genus, abundance, thresholds=None, kb_normalized=None):

    if thresholds is None:
        thresholds = {"high": 10, "low": 1}
    if kb_normalized is None:
        kb_normalized = {}
    thresholds_high = thresholds.get("high", 10)
    thresholds_mod = thresholds.get("moderate", 1)
    # Normalize genus name for lookup
    key = genus.lower().replace(" ", "_").replace("-", "_")

    # 1️⃣ Direct match
    if key in kb_normalized:
        entry = kb_normalized[key]
    else:
        # 2️⃣ Fuzzy match with safety checks
        from fuzzywuzzy import process
        result = process.extractOne(key, kb_normalized.keys()) if kb_normalized else None
        if result:
            match, score = result
            if match and score > 80:
                entry = kb_normalized[match]
            else:
                return None
        else:
            return None

    if abundance >= thresholds_high:
        level_text = entry.get("high_abundance", "")
    elif abundance >= thresholds_mod:
        level_text = entry.get("moderate_abundance", "")
    else:
        level_text = entry.get("low_abundance", "")

    return {
        "overview": entry.get("overview", ""),
        "role": entry.get("role", ""),
        "level_text": level_text,
        "recommendations": entry.get("recommendations", []),
    }

    """Fetch bacterial information from the knowledge base (with fuzzy matching)."""
    genus_norm = normalize_genus(genus)
    info = kb_normalized.get(genus_norm) if kb_normalized else None

    # 🔍 Fuzzy match if exact name not found
    if not info and kb_normalized:
        for kb_key in kb_normalized.keys():
            if genus_norm in kb_key or kb_key in genus_norm:
                info = kb_normalized[kb_key]
                print(f"[INFO] Matched genus '{genus}' ↔ '{kb_key}' in KB.")
                break

    if not info:
        return None

    # Determine level based on abundance
    high_thr = thresholds.get("high", 10.0) if thresholds else 10.0
    mod_thr = thresholds.get("moderate", 1.0) if thresholds else 1.0
    level = (
        "high" if abundance >= high_thr
        else "moderate" if abundance >= mod_thr
        else "low"
    )

    level_text = f"{genus.capitalize()} is present at {abundance:.2f}%, which is considered {level} abundance."

    return {
        "overview": info.get("overview", "No overview available."),
        "role": info.get("role", "No role information."),
        "recommendations": info.get("recommendations", []),
        "level": level,
        "level_text": level_text
        
    }







def normalize_genus(name: str) -> str:
    if not isinstance(name, str):
        return ""
    return name.strip().lower().replace("/", "_").replace("-", "_").replace(" ", "_")


for needed_col, df in [("Sample", diversity_df), ("Sample", abundance_df), ("Genus", abundance_df)]:
    if needed_col not in df.columns:
        df[needed_col] = ""


diversity_df["Sample_norm"] = diversity_df["Sample"].apply(normalize_id)
abundance_df["Sample_norm"] = abundance_df["Sample"].apply(normalize_id)
abundance_df["Genus_norm"] = abundance_df["Genus"].apply(normalize_genus)


# --- Safety: ensure meta_df exists even if JSON doesn't define it ---
if 'meta_df' not in locals():
    print("[WARN] meta_df not defined — creating placeholder.")
    meta_df = pd.DataFrame({"Sample": []})


meta_df["Sample_norm"] = meta_df["Sample"].apply(normalize_id)
prebiotic_targets = {
    "Inulin": ["Bifidobacterium", "Faecalibacterium", "Roseburia"],
    "GOS": ["Bifidobacterium", "Lactobacillus"],
    "PHGG": ["Faecalibacterium", "Roseburia", "Subdoligranulum"],
    "Resistant Starch": ["Ruminococcus", "Faecalibacterium", "Eubacterium"]
}


# ========= Main Report Function =========
def build_microbiome_resilience_radar(summary_dict, output_path):
    """
    Builds a 3-axis radar (spider) chart showing balance of good, variable, and bad bacteria.
    Saves a small PNG for embedding into the report.
    """
    categories = ['Good', 'Variable', 'Bad']
    values = [
        summary_dict['Good_Score'],
        summary_dict['Variable_Score'],
        summary_dict['Bad_Score']
    ]

    # Close the polygon by repeating the first value
    values += values[:1]
    angles = np.linspace(0, 2 * np.pi, len(categories) + 1)

    fig, ax = plt.subplots(figsize=(3.5, 3.5), subplot_kw=dict(polar=True))

    # Plot data
    ax.plot(angles, values, color='#81C784', linewidth=2)
    ax.fill(angles, values, color='#A5D6A7', alpha=0.4)

    # Labels
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=10, fontweight='bold')
    ax.set_yticks([2, 4, 6, 8, 10])
    ax.set_yticklabels(['2', '4', '6', '8', '10'], color='gray', size=8)
    ax.set_ylim(0, 10)
    ax.grid(color='lightgrey', linestyle='dotted', linewidth=0.7)

    plt.title("Microbiome Resilience Radar", fontsize=12, pad=15, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=120, bbox_inches='tight')
    plt.close(fig)

from reportlab.platypus import Table, TableStyle

from reportlab.platypus import Table, TableStyle, Paragraph
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle

from reportlab.platypus import Table, TableStyle, Paragraph, Spacer
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm

from reportlab.platypus import Table, TableStyle, Paragraph, Spacer
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm

def make_biomarker_table(df_grouped, table_title, styles, genus_subset=None, kb_normalized=None, thresholds=None):
    """
    Builds a clean, wrapped biomarker summary table showing all detected bacteria.
    """

    thresholds = thresholds or {"high": 5.0, "moderate": 1.0}
    high_thr = thresholds.get("high", 5.0)
    mod_thr = thresholds.get("moderate", 1.0)

    # --- Styles ---
    header_style = ParagraphStyle(
        "hdr",
        parent=styles["Body"],
        alignment=1,
        fontSize=9,
        leading=11,
        textColor=colors.white
    )
    cell_style = ParagraphStyle(
        "cell",
        parent=styles["Body"],
        fontSize=8,
        leading=10.5,
        wordWrap="CJK",
        spaceAfter=2,
    )

    # --- Normalize input list ---
    def norm(s): return str(s).strip().lower().replace("/", "_").replace("-", "_")
    genus_subset_norm = [norm(g) for g in genus_subset] if genus_subset else []

    # --- Filter dataframe but KEEP order ---
    df_grouped["Genus_norm"] = df_grouped["Genus"].apply(norm)

    if genus_subset_norm:
        rows_df = df_grouped[df_grouped["Genus_norm"].apply(lambda g: any(g.startswith(b) for b in genus_subset_norm))]
    else:
        rows_df = df_grouped.copy()

    # --- Convert to list of dicts (preserves order fully) ---
    biomarker_rows = rows_df.to_dict(orient="records")

    # --- Debug output ---
    dbg_list = [r["Genus"] for r in biomarker_rows]
    print(f"🧫 {table_title}: {len(dbg_list)} genera → {dbg_list}")

    # --- Build table header ---
    data = [[
        Paragraph("<b>Bacterium</b>", header_style),
        Paragraph("<b>Abundance (%)</b>", header_style),
        Paragraph("<b>Why it matters</b>", header_style),
        Paragraph("<b>Status</b>", header_style)
    ]]

    # --- Populate data ---
    for r in biomarker_rows:
        genus = r["Genus"]
        abundance = float(r["Abundance"])
        key = norm(genus)


        info = kb_normalized.get(key, {}) if kb_normalized else {}

        # Normalize genus name before lookup
        genus_name_norm = normalize_genus(key)
        info = kb_normalized.get(genus_name_norm, {}) if kb_normalized else {}

        if not info:
            print(f"[WARN] No summary for genus: {key} (lookup: {genus_name_norm})")


        # Determine description based on abundance thresholds
        if abundance >= high_thr:
            why_text = info.get("high_abundance", "High presence — supports gut balance and diversity.")
        elif abundance >= mod_thr:
            why_text = info.get("moderate_abundance", "Moderate presence — contributes to a balanced microbiome.")
        else:
            why_text = info.get("low_abundance", "Low presence — may benefit from dietary or prebiotic support.")

        # Determine status + color
        if abundance >= high_thr:
            status = "Balanced"; bg = colors.Color(0.85, 0.95, 0.85)
        elif abundance >= mod_thr:
            status = "Moderate"; bg = colors.Color(1.00, 0.96, 0.80)
        else:
            status = "Needs Support"; bg = colors.Color(1.00, 0.90, 0.90)

        # Append row
        data.append([
            Paragraph(f"<b><i>{genus}</i></b>", cell_style),
            Paragraph(f"{abundance:.2f}", cell_style),
            Paragraph(why_text, cell_style),
            Paragraph(f"<para backColor='{bg.hexval()}'><b>{status}</b></para>", cell_style)
        ])

    # --- Build table ---
    table = Table(data, colWidths=[3*cm, 3*cm, 8*cm, 3*cm])
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#003366")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
        ("BOX", (0, 0), (-1, -1), 0.3, colors.grey),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))

    # Zebra striping for clarity
    for i in range(1, len(data)):
        if i % 2 == 0:
            table.setStyle(TableStyle([("BACKGROUND", (0, i), (-1, i), colors.whitesmoke)]))

    return [
        Paragraph(f"<b>{table_title}</b>", styles["Section"]),
        Spacer(1, 6),
        table,
        Spacer(1, 12),
    ]

    return [table, Spacer(1, 10)]


from reportlab.graphics.shapes import Drawing, String
from reportlab.graphics.charts.piecharts import Pie
from reportlab.lib import colors


from reportlab.graphics.shapes import Drawing, String
from reportlab.graphics.charts.piecharts import Pie
from reportlab.lib import colors

def make_biomarker_pie(good_bugs, bad_bugs, variable_bugs):
    """
    Creates a pastel pie chart summarizing good / bad / variable bacteria proportions.
    Centered under title.
    """
    total = sum([good_bugs, bad_bugs, variable_bugs])
    if total == 0:
        total = 1

    data = [
        round((good_bugs / total) * 100, 1),
        round((bad_bugs / total) * 100, 1),
        round((variable_bugs / total) * 100, 1),
    ]
    labels = [
        f"Good Bugs ({data[0]}%)",
        f"Bad Bugs ({data[1]}%)",
        f"Variable Bugs ({data[2]}%)"
    ]

    pastel_colors = [
        colors.HexColor("#B6E2A1"),  # soft green
        colors.HexColor("#F5B7B1"),  # coral pink
        colors.HexColor("#F9E79F")   # soft yellow
    ]

    # --- Drawing container ---
    d = Drawing(400, 250)  # increased width slightly for center balance

    pie = Pie()
    pie.x = 100             # move slightly right
    pie.y = 40              # move a bit lower under title
    pie.width = 200
    pie.height = 160
    pie.data = data
    pie.labels = labels
    pie.sideLabels = True
    pie.slices.strokeWidth = 0.5
    pie.slices.strokeColor = colors.white

    for i, color in enumerate(pastel_colors):
        pie.slices[i].fillColor = color
        pie.slices[i].popout = 0 if data[i] < max(data) else 8  # highlight largest group

    # --- Add title centered above chart ---
    d.add(String(200, 220, "Microbial Biomarker Composition",
                 fontSize=11, fontName="Helvetica-Bold",
                 textAnchor="middle", fillColor=colors.HexColor("#333333")))

    d.add(pie)
    return d



def generate_report(sample_id):

    global prebiotic_targets
    from datetime import datetime
    import os

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pdf_path = os.path.join(args.output, f"{args.sample}_report_{timestamp}.pdf")
    txt_path = os.path.join(args.output, f"{args.sample}_summary_{timestamp}.txt")

    # Utility: normalize sample names safely
    import re
    # --- Normalize sample columns in all dataframes ---
    import re

    def normalize_id(x):
        """Normalize sample IDs exactly like JSON normalization."""
        if isinstance(x, str):
            x = x.lower().replace("-", "_").replace(".", "_")
            x = re.sub(r'[^a-z0-9]+', '_', x)
            return x.strip('_')
        return x

    if "Sample" in abundance_df.columns:
        abundance_df["Sample_norm"] = abundance_df["Sample"].apply(normalize_id)
    if "Sample" in diversity_df.columns:
        diversity_df["Sample_norm"] = diversity_df["Sample"].apply(normalize_id)
    if "Sample" in meta_df.columns:
        meta_df["Sample_norm"] = meta_df["Sample"].apply(normalize_id)

    def normalize_id(x):
        """Normalize sample IDs exactly like JSON normalization."""
        if isinstance(x, str):
            x = x.lower().replace("-", "_").replace(".", "_")
            x = re.sub(r'[^a-z0-9]+', '_', x)
            return x.strip('_')
        return x

    sample_id_norm = normalize_id(sample_id)
    print(f"\n📘 Generating report for: {sample_id_norm}")

    sample_abundance = abundance_df[abundance_df["Sample_norm"] == sample_id_norm].copy()
    sample_diversity = diversity_df[diversity_df["Sample_norm"] == sample_id_norm].copy()

    if sample_diversity.empty:
        print(f"[WARN] Diversity data not found for '{sample_id}' (normalized as '{sample_id_norm}').")
    else:
        print(f"✅ Found diversity data for '{sample_id_norm}'")

    print("\n🧠 DEBUG: Checking sample matching inside generate_report()")
    print("Available normalized sample IDs in abundance_df:")
    print(abundance_df["Sample_norm"].unique())
    print(f"Looking for: {sample_id_norm}")

    if sample_abundance.empty:
        print(f"❌ No rows found for '{sample_id_norm}' — mismatch likely in normalization or dataframe overwrite.")
        import sys;
        sys.exit(1)
    else:
        print(f"✅ Found {len(sample_abundance)} abundance rows for '{sample_id_norm}'")

    # define filename early

    pdf_name = f"Your_Poopie_Report_{sample_id.replace('.', '_')}.pdf"

   # pdf_name = f"{timestamp}_{sample_id}.pdf"


    # ===== Initialize styles, story, and TOC early =====
    styles = getSampleStyleSheet()
    styles["Title"].fontName = "Helvetica-Bold"
    styles["Title"].fontSize = 20
    styles.add(ParagraphStyle(name="Body", fontName="Helvetica", fontSize=10.5, leading=14))
    styles.add(ParagraphStyle(name="Section", fontName="Helvetica-Bold", fontSize=13, spaceAfter=6))
    if "TOCHeading" not in styles:
        styles.add(ParagraphStyle(name="TOCHeading", parent=styles["Section"]))
    from reportlab.platypus import Image

    # --- Title Page (Cover) ---
    def add_title_page(story, styles):
        """Creates a friendly, minimal title page with the poop-smiley logo."""
        story.append(Spacer(1, 4 * cm))
        story.append(Spacer(1, 4 * cm))
        # 💩 Add poop logo (centered). Make sure 'poop_logo.png' exists in the same directory.
        story.append(Paragraph(
            "<font size=18><b>Poopie</b></font>", styles["Title"]
        ))

        story.append(Spacer(1, 1 * cm))
        story.append(Paragraph(
            "<font size=14><b>Happy Gut. Happy Health</b></font>", styles["Title"]
        ))
        story.append(Spacer(1, 0.5 * cm))
        # 🌿 Title text
        story.append(Paragraph(
            "<font size=12><b>Your Personalized Microbiome Report</b></font>", styles["Title"]
        ))
        story.append(Spacer(1, 0.5 * cm))


        story.append(Spacer(1, 6 * cm))
        story.append(Spacer(1, 4 * cm))

        # Optional footer on title page
        story.append(Paragraph(
            "<font size=9 color='#888888'>© 2025 Poopie.uk — Confidential & Educational Use Only</font>",
            styles["Body"]
        ))

        story.append(PageBreak())

    story = []
    # --- Cover Page ---
    add_title_page(story, styles)

    # ===== Table of Contents setup =====
    from reportlab.platypus.tableofcontents import TableOfContents



    toc = TableOfContents()
    toc.levelStyles = [
        ParagraphStyle(name="TOCLevel1", fontName="Helvetica", fontSize=10.5, leftIndent=0, firstLineIndent=0,
                       spaceBefore=4, spaceAfter=2, leading=12),
        ParagraphStyle(name="TOCLevel2", fontName="Helvetica", fontSize=9.5, leftIndent=15, firstLineIndent=0,
                       spaceBefore=2, spaceAfter=1, leading=11),
    ]
    toc.dotsMinLevel = 0
    toc.dotsMaxLevel = 2



    # --- Initialize TOC tracking ---
    toc_started = False
    doc = None
    intro_page_number = None
    #canvas_obj = None

    # --- Hook to auto-add sections to Table of Contents ---
    def after_flowable(flowable):
        """Register section headings for TOC after 'Introduction' and mark page anchor."""
        nonlocal toc_started, intro_page_number

        from reportlab.platypus import Paragraph as RLParagraph
        if isinstance(flowable, RLParagraph):
            style_name = getattr(flowable.style, "name", "")
            text = flowable.getPlainText().strip()

            # Start TOC tracking when Introduction section appears
            if style_name == "Section" and text.lower() == "introduction":
                toc_started = True
                if intro_page_number is None:
                    intro_page_number = getattr(doc.canv, "_pageNumber", doc.page)

            # Add TOC entries for sections after Introduction
            if toc_started and style_name == "Section" and text.lower() != "table of contents":
                page_num = getattr(doc.canv, "_pageNumber", doc.page)
                toc.addEntry(level=0, text=text, pageNum=page_num)



    # --- Helper for consistent section headings ---

        # Title Page

        # TOC Page
    story.append(Paragraph("<b>Table of Contents</b>", styles["Section"]))
    story.append(divider())
    story.append(toc)
    story.append(PageBreak())

    # --- Helper for creating section headings ---
    def add_section(title_text, story, styles):
        story.append(Paragraph(title_text, styles["Section"]))
        story.append(divider())

        #toc.addEntry(level=level, text=title_text, pageNum=len(story))

    # You can now continue building your PDF below:
    # (Add intro, diversity, probiotic sections, etc.)



    sample_diversity = diversity_df[diversity_df["Sample_norm"] == sample_id_norm]

    # --- Safe sample matching ---
    def normalize_id(x):
        """Normalize sample names for robust comparison."""
        if isinstance(x, str):
            return x.lower().replace("-", "_").replace(".", "_")
        return x

    def normalize_id(x):
        if isinstance(x, str):
            return x.lower().replace("-", "_").replace(".", "_")
        return x

    # Normalize both JSON and sample_id names before matching
    diversity_df["Sample_norm"] = diversity_df["Sample"].apply(normalize_id)
    sample_norm = normalize_id(sample_id)

    sample_diversity = diversity_df[diversity_df["Sample_norm"] == sample_norm]

    # --- Overview diversity indices (safe) ---
    try:
        if sample_diversity.empty:
            raise IndexError
        row = sample_diversity.iloc[0]
    except IndexError:
        print(f"[WARN] Diversity data missing for '{sample_id}'. Using defaults.")
        row = pd.Series({"Observed": 0, "Shannon": 0.0, "Simpson": 0.0})

    #sample_diversity = diversity_df[diversity_df["Sample_norm"] == sample_id_norm]

    # --- Collapse sample abundance to genus level ---
    sample_abundance = (
        abundance_df[abundance_df["Sample_norm"] == sample_id_norm]
        .copy()
    )

    # Standardize genus names (remove species-level info)
    sample_abundance["Genus"] = sample_abundance["Genus"].str.split().str[0]

    # Collapse duplicates to get total genus-level abundance
    sample_abundance = (
        sample_abundance.groupby("Genus", as_index=False)["Abundance"].sum()
        .sort_values(by="Abundance", ascending=False)
        .reset_index(drop=True)
    )

    sample_meta = meta_df[meta_df["Sample_norm"] == sample_id_norm]

    if sample_abundance.empty:
        print(f"⚠️ No abundance rows found for sample '{sample_id_norm}'")
        sys.exit("❌ No matching sample in abundance data.")

    # --- Overview diversity indices ---


    #row = sample_diversity.iloc[0]

    observed = int(row.get("Observed", 0))
    shannon = float(row.get("Shannon", 0.0))
    simpson = float(row.get("Simpson", 0.0))

    if shannon >= 5.0:
        diversity_state = "high diversity — a rich and balanced gut ecosystem"
    elif shannon >= 3.5:
        diversity_state = "moderate diversity — healthy but with room for more variety"
    else:
        diversity_state = "low diversity — may benefit from more plant fibre and diversity"

    # ========= Intro Text =========
    intro_text = """
        <b>Welcome to your microbiome report</b><br/><br/>
        Thank you for taking the time to explore your gut health. Your gut microbiome is a living ecosystem made up of trillions of bacteria that support digestion, immunity, energy and mood. Testing your microbiome helps you understand which species are living inside you and how they may be working together.<br/><br/>
        Your gut microbiome is the community of trillions of bacteria, fungi, and other microorganisms that live inside your digestive tract.
        Testing it helps you understand how your gut is functioning and what may be influencing your overall well-being. <br/><br/>

        <b>1. Discover your unique balance.</b><br/>
        Everyone's gut ecosystem is different. Microbiome testing shows which helpful and opportunistic bacteria are present and how they compare to balanced patterns seen in healthy populations.
        <br/><b>2. Understand digestion and comfort.</b><br/>
        Some bacteria produce short-chain fatty acids that nourish the gut lining, while others can contribute to gas or irritation if they grow too much. A microbiome test helps identify where your gut might need more support.
        <br/><b>3. Support nutrition and lifestyle choices.</b><br/>
        The results can guide food and lifestyle recommendations - such as adding more fibre, prebiotics, or fermented foods - tailored to the actual composition of your gut bacteria.
        <br/><b>4. Track changes over time.</b><br/>
        Testing periodically allows you to see how your gut community responds to diet, stress, travel, or medication. It is a way to monitor progress and maintain balance.
        <br/><b>5. Contribute to preventive health.</b><br/>
        While microbiome testing isn't diagnostic, understanding your gut ecosystem supports informed discussions with nutrition or healthcare professionals about digestion, immunity, and overall wellness.<br/><br/>

        <b>In short:</b>
        This test has been developed by leading scientists, using the latest advances in gut microbiome science. It provides a detailed look at your **gut wellness** by examining the balance, diversity, and activity of the microorganisms that live in your digestive system.
        <br/>Your microbiome is assessed across several key **health categories**, each showing how well a different part of your gut ecosystem is functioning:<br/><br/>
         • Diversity <br/>Measures how many different types of microbes live in your gut and how evenly they are balanced. Greater diversity is usually linked to resilience and stable digestion.<br/><br/>
         • Probiotics <br/>Looks at beneficial bacteria such as *Bifidobacterium* and *Lactobacillus* that help protect your gut lining, support digestion, and keep harmful species under control.<br/><br/>
         • Prebiotics <br/>Assesses how well your gut bacteria are being nourished by the fibres and plant foods in your diet. These “good bacteria foods” help maintain balance and fuel healthy activity.<br/><br/>
         • Postbiotics <br/>Examines the beneficial compounds your microbes produce, such as short-chain fatty acids, which calm the gut, feed the intestinal lining, and support immunity.<br/><br/>
         • Pathogens & Opportunists <br/>Identifies bacteria that can become problematic if they overgrow. Detecting these early helps you and your practitioner support balance before they cause symptoms.<br/><br/>
         • Mental Wellness <br/>Focuses on bacteria associated with the production of neurotransmitters like serotonin and GABA, which influence mood, focus, and emotional balance.<br/><br/>
         • Gut–Skin Connection <br/>Looks at bacteria linked to skin health, hydration, and inflammation — highlighting how gut balance can influence issues such as dryness, breakouts, or sensitivity.<br/><br/>

    <b>What Each Section Means</b><br/><br/>

    Each category in your microbiome report focuses on a unique aspect of gut wellness. Together, they provide a complete picture of how your digestive and microbial systems support your overall health.<br/><br/>

    Diversity <br/> 
    Diversity describes how many different types of microbes live in your gut and how evenly they share space. A rich and balanced gut community supports flexibility, comfort, and resilience against imbalance. Lower diversity can sometimes reflect recent antibiotic use, stress, or a low-fibre diet, while higher diversity often points to a well-nourished microbiome.<br/><br/>

    Probiotics  <br/>
    Probiotics are the beneficial bacteria that help maintain a healthy gut lining and protect you from unwanted microbes. They play key roles in digestion, nutrient absorption, and immune balance. Natural probiotic species such as *Bifidobacterium* and *Lactobacillus* thrive when you include fermented foods like yogurt, kefir, kimchi, or sauerkraut in your diet.<br/><br/>

    Prebiotics <br/>
    Prebiotics are special fibres and plant compounds that feed your good bacteria. They pass through the digestive tract and reach the colon, where microbes ferment them into nutrients that strengthen your gut barrier and reduce irritation. Examples include inulin, GOS, PHGG, and resistant starch — all found in fruits, vegetables, legumes, and whole grains.<br/><br/>

    Postbiotics  <br/>
    Postbiotics are the beneficial substances produced when your microbes break down prebiotic fibres. These include short-chain fatty acids such as butyrate, acetate, and propionate, which nourish gut cells, calm inflammation, and support immunity. A well-fed microbiome naturally produces more of these healing compounds.<br/><br/>

    Pathogens & Opportunists<br/>
    This section highlights microbes that can cause discomfort or imbalance if they grow in excess. They are normal members of the gut community but are best kept in check by strong beneficial bacteria and healthy diet patterns. Your report flags these so that your practitioner can support rebalancing if needed. <br/><br/>

    Mental Wellness  <br/>
    Certain gut microbes help produce neurotransmitters such as serotonin and GABA, which influence mood, stress response, and focus. An imbalanced microbiome can sometimes contribute to fatigue or low mood, while a balanced one helps maintain emotional wellbeing and clarity. The brain has a direct effect on the stomach and intestines. For example, the very thought of eating can release the stomach's juices before food gets there. This connection goes both ways. A troubled intestine can send signals to the brain, just as a troubled brain can send signals to the gut. Therefore, a person's stomach or intestinal distress can be the cause or the product of anxiety, stress, or depression. That's because the brain and the gastrointestinal (GI) system are intimately connected.<br/><br/>

    This is especially true in cases where a person experiences gastrointestinal upset with no obvious physical cause. For such functional GI disorders, it is difficult to try to heal a distressed gut without considering the role of stress and emotion.<br/><br/>

    Gut–Skin Connection  <br/>
    Your gut and skin are closely linked through the immune and detoxification systems. When the gut barrier is balanced, the skin often appears clearer and calmer. This section examines bacteria that influence inflammation, hydration, and skin health from within.<br/><br/>

    Together, these sections show how your microbiome interacts with every part of your wellbeing — digestion, energy, mood, and even skin.  <br/><br/>
    Understanding each one helps you and your Nutritional Therapist create the most effective, personalised plan for lasting gut balance.
        Remember that your microbiome is dynamic and can shift with diet, lifestyle and stress. Use these insights as a guide to support balance and resilience, not as a diagnosis.<br/><br/>
        <font color="#777777" size="9">This report is provided for educational and wellness purposes only. It is not intended to diagnose or treat any medical condition.</font>

        """

    diversity_text = f"""
        <b>Your Gut Microbiome Diversity</b><br/><br/>
        Your gut microbiome is like a bustling ecosystem — filled with trillions of tiny organisms that help digest food, produce nutrients, and keep your gut lining strong. Just like a healthy garden or rainforest needs many different plants to stay balanced, your gut also thrives when it has a rich and varied community of microbes. When we analyse your sample, we look at how many different kinds of bacteria are present and how evenly they live together. Together, these factors describe your gut microbiome diversity — an important reflection of gut balance and resilience.<br/><br/>

        1. <b>Richness</b> — How many kinds of microbes are there?<br/><br/>

        Richness tells us how many different bacterial types we detect in your gut. A gut with higher richness has more variety — this often means it’s more adaptable and better at handling changes in diet or stress. Lower richness can occur after antibiotics, illness, or limited-fibre diets, when some helpful bacteria temporarily drop in number. Think of richness as variety on your plate — the more ingredients, the more balanced and nourishing your meal. <br/><br/>

        2. <b>Evenness</b> — How balanced is your gut community?<br/><br/>

        Evenness looks at how evenly your bacteria share space. If one or two types dominate, the gut may be less stable. When the population is more evenly shared, your microbiome works in harmony — different microbes play different roles to support comfort and digestion. A balanced microbiome is like a well-tended garden — no single plant overgrows the others, and everything flourishes together.<br/><br/>

        3. <b>Overall Diversity</b> — Your microbiome’s balance score<br/><br/>

        Diversity combines both richness and evenness into one measure. It reflects the overall health and adaptability of your gut ecosystem. We often express this using the Shannon Diversity Index — a number that increases when your gut has both many species and a balanced distribution between them.<br/><br/>

        High diversity → usually linked with strong gut resilience and better digestive comfort.<br/>

        Moderate diversity → common and healthy in many people.<br/>

        Low diversity → may suggest reduced flexibility, often seen after antibiotics, stress, or restricted diets.<br/><br/>

        In simple terms, your diversity score is like an ecosystem rating — it shows how varied and balanced your gut community is.<br/><br/>

        <b>The Science Behind the Numbers</b><br/><br/>

        We use DNA sequencing of your stool sample to identify the bacteria living in your gut.
Each read is matched to known microbial genomes, allowing us to estimate which species (or genera) are present and how abundant they are.<br/>

        <br/>Then, we calculate several key metrics:<br/><br/>

        <b>Observed Species (Richness):</b><br/>

        Counts how many different bacteria are detected in your sample.
More species = greater microbial variety.<br/><br/>
        <b>Shannon Diversity Index:</b><br/>

        Balances both richness and evenness — it increases when your gut has many species and they are fairly balanced in abundance.
Higher values (typically 4–6) indicate a healthy, resilient microbiome.<br/><br/>
<b>Simpson Index (Evenness):</b><br/>

        Focuses on how evenly the species share space.
If one or two species dominate, this number drops.<br/><br/>
        These indices are then scaled to a 0–10 range to make them easy to interpret in your personalized report.<br/><br/>

        <b>Supporting a Diverse Gut</b><br/><br/>

        Diversity is dynamic — it can improve over time with consistent habits. To nurture a rich and balanced microbiome:<br/><br/>

        Eat a wide variety of plant-based foods (aim for 20–30 different types per week).<br/>

        Include fermented foods like yogurt, kefir, kimchi, or kombucha.<br/>

        Stay hydrated and manage stress — both strongly influence microbial balance.<br/><br/>

        A diverse gut is a resilient gut — when your microbes thrive, they support your digestion, mood, and overall wellbeing.<br/><br/>


        
        """

    # ========= Probiotic Detection & Diversity Score =========
    # ========= Probiotics =========


    # Dynamic probiotic score
    probiotic_targets = ["Lactobacillus", "Bifidobacterium", "Faecalibacterium", "Roseburia", "Adlercreutzia","Akkermansia","Lactococcus","Christensenella","Coprococcus","Barnesiella","Oxalobacter"]
    probiotic_score = compute_dynamic_score(sample_abundance, probiotic_targets, ideal_total=10)
    #story.append(Paragraph(f"<b>Probiotic Balance Score:</b> {probiotic_score}", styles["Body"]))
    #story.append(CleanGauge(probiotic_score, width=420, height=14, mode="score", label=f"{probiotic_score:.1f}"))
    #story.append(Spacer(1, 12))

    probiotic_section_text = f"""
            Probiotics are beneficial live microorganisms that support your digestive and immune systems when present in healthy numbers. 
            They are often called <i>"friendly bacteria"</i> because they help maintain balance inside your gut ecosystem, where many species — 
            both helpful and harmful — naturally coexist.<br/><br/>

            Your gut is like a living community, and probiotics act as peacekeepers. 
            They help crowd out unwanted bacteria, strengthen the gut barrier, and produce compounds that soothe inflammation. 
            You can naturally increase probiotic levels through fermented foods such as <i>kefir, yoghurt, kimchi, sauerkraut</i>, and <i>kombucha</i>.<br/><br/>

            These beneficial microbes are key players in many body systems. They:
            <ul>
            <li>Promote smooth digestion and nutrient absorption</li>
            <li>Produce essential postbiotics like butyrate and B-vitamins</li>
            <li>Help balance the immune system and reduce gut inflammation</li>
            <li>Communicate with the brain to support mood, energy, and focus</li>
            <li>Contribute to clear skin and overall metabolic health</li>
            </ul>

            We identify probiotic bacteria in your sample using <b>16S rRNA sequencing technology</b>, 
            which pinpoints key beneficial genera such as <i>Lactobacillus</i>, <i>Bifidobacterium</i>, <i>Faecalibacterium</i>, 
            <i>Roseburia</i>, <i>Streptococcus</i>, and <i>Akkermansia</i>.<br/><br/>

            """
    # --- Detect probiotic genera present in this sample ---
    probiotic_genera = ["Lactobacillus", "Bifidobacterium", "Faecalibacterium", "Roseburia", "Ruminococcus",
                        "Coprococcus"]
    detected_probiotics = list(
        sample_abundance[sample_abundance["Genus"].isin(probiotic_genera)]["Genus"].unique()
    )

    # --- Compute probiotic diversity score dynamically ---
    probiotic_diversity_score = compute_dynamic_score(sample_abundance, probiotic_genera, ideal_total=10)

    print(f"🦠 Detected probiotics: {detected_probiotics}")
    print(f"🌿 Probiotic Diversity Score: {probiotic_diversity_score}/10")

    prebiotic_section_text = f"""
                <b>Understanding Prebiotics</b><br/><br/>
                Your gut microbes need the right kind of “food” to stay active and balanced. Prebiotics are special types of dietary fibres and plant compounds that feed the beneficial bacteria in your gut — helping them grow, multiply, and produce nutrients that keep your digestive system healthy.<br/><br/>

                While probiotics are live bacteria, <b>prebiotics are what those good bacteria eat.</b> They pass through your stomach undigested and reach the large intestine, where microbes break them down to produce short-chain fatty acids (such as butyrate) that nourish the gut lining, support immunity, and calm inflammation.<br/><br/>

                <b>1. Why prebiotics matter</b><br/><br/>

        <i>Fuel for good bacteria</i>: Prebiotics selectively feed beneficial species like Faecalibacterium, Bifidobacterium, and Roseburia. <br/><br/>

        <i>Gut barrier support</i>: The compounds they help create strengthen the intestinal wall and reduce irritation.<br/><br/>

        <i>Digestive comfort</i>: A well-fed microbiome promotes regularity and reduces gas or bloating over time.<br/><br/>

        <i>Whole-body impact</i>: Balanced microbial activity influences mood, metabolism, and immune function.<br/><br/>

        <b>3. How we assess prebiotic response</b>

        In your report, we look at the presence of bacteria that thrive on specific prebiotics — for example, whether your sample shows strong growth of Bifidobacterium (responds to GOS) or Faecalibacterium (responds to resistant starch). This helps guide which fibres your gut may benefit from most.<br/><br/>

        <b>4. Supporting your gut with prebiotics</b><br/><br/>

        A balanced intake of prebiotics helps your beneficial bacteria flourish —
        supporting digestion, energy, mood, and long-term gut health.<br/><br/>

                """
    postbiotic_section_text = f"""
            You may have heard of probiotics (the good bacteria) and prebiotics (the food that nourishes them). But there’s also a third — and equally important — part of this picture: postbiotics. Postbiotics are the beneficial substances produced when your gut bacteria digest fibre and other nutrients. They include compounds like short-chain fatty acids (butyrate, acetate, propionate), enzymes, vitamins (such as B12 and K), peptides, and other small molecules that directly support gut and whole-body health. Think of postbiotics as the healthy “by-products” of a well-fed microbiome — the proof that your beneficial bacteria are working for you.<br/><br/>
            <b>1. Why postbiotics matter</b><br/><br/>

        <i>Nourish the gut lining</i>: Compounds like butyrate act as fuel for intestinal cells, keeping the gut barrier strong and reducing inflammation.<br/><br/>

        <i>Balance the immune system</i>: Postbiotics help the immune system stay alert but not overreactive, supporting long-term digestive comfort.<br/><br/>

        <i>Support metabolism and mood</i>: Some postbiotic molecules communicate with other organs — influencing energy balance and even brain health.<br/><br/>

        <i>Promote resilience</i>: A gut that makes plenty of postbiotics tends to recover faster from stress, illness, or antibiotic use.<br/><br/>

        <b>2. How postbiotics are formed</b><br/><br/>

        When beneficial bacteria (such as Faecalibacterium, Roseburia, or Bifidobacterium) ferment prebiotic fibres, they produce short-chain fatty acids and other helpful compounds.<br/><br/>

        These postbiotics then:<br/>

        Strengthen the intestinal wall<br/>

        Lower the gut’s pH (making it less friendly for harmful microbes)<br/>

        Help regulate bowel movements and nutrient absorption<br/><br/>

        In essence, a healthy postbiotic profile means your microbes are active and well-fed.<br/><br/>

        <b>3. What influences postbiotic production</b><br/><br/>

        Your postbiotic levels depend on:<br/>

        Diet quality and fibre variety (the more plant fibres, the more fuel for bacteria)<br/>

        Microbial diversity (a rich and balanced microbiome produces more beneficial compounds)<br/>

        Hydration and gut transit time (enough water keeps fermentation processes healthy)<br/>

        Even moderate changes — such as increasing vegetables, whole grains, and legumes — can boost postbiotic activity in just a few weeks.<br/><br/>

        <b>4. Supporting healthy postbiotic production</b><br/><br/>


        When your bacteria are happy and well-fed, they create postbiotics that in turn nourish you — supporting digestion, energy, and immune balance.<br/><br/>
            """

    # ========= Section Setup =========
    styles = getSampleStyleSheet()
    styles["Title"].fontName = "Helvetica-Bold"
    styles["Title"].fontSize = 20
    styles.add(ParagraphStyle(name="Body", fontName="Helvetica", fontSize=10.5, leading=14))
    styles.add(ParagraphStyle(name="Section", fontName="Helvetica-Bold", fontSize=13, spaceAfter=6))



    prebiotic_genera = ["Faecalibacterium", "Roseburia", "Subdoligranulum", "Bifidobacterium"]
    prebiotic_score = compute_dynamic_score(sample_abundance, prebiotic_targets["Inulin"] +
                                            prebiotic_targets["GOS"] +
                                            prebiotic_targets["PHGG"] +
                                            prebiotic_targets["Resistant Starch"], ideal_total=8)
    #story.append(Paragraph(f"<b>Prebiotic Support Score:</b> {prebiotic_score}", styles["Body"]))
    #story.append(CleanGauge(prebiotic_score, width=420, height=14, mode="score", label=f"{prebiotic_score:.1f}"))
    #story.append(Spacer(1, 12))


    add_section("Introduction", story, styles)

    story.append(Paragraph(intro_text, styles["Body"]))
    story.append(PageBreak())
    story.append(PageBreak())
    # --- If genus_df available, add a Top-10 chart for the chosen sample ---
    if 'genus_df' in globals() and genus_df is not None:
        lookup_id = args.sample  # this should match Sample index in genus_df (pre-normalization)
        if lookup_id in genus_df.index:
            story.append(PageBreak())
            story.append(Paragraph("Genus Abundance Overview", styles["Section"]))
            story.append(divider())

            top10 = genus_df.loc[lookup_id].sort_values(ascending=False).head(10)
            story.append(Paragraph(f"Top 10 most abundant genera in sample {lookup_id}:", styles["Body"]))

            fig, ax = plt.subplots(figsize=(6, 3))
            top10[::-1].plot(kind="barh", ax=ax)  # default style/colors
            ax.set_xlabel("Relative Abundance (%)")
            ax.set_ylabel("Genus")
            plt.tight_layout()
            img_path = f"/tmp/top10_{lookup_id.replace('/', '_')}.png"
            plt.savefig(img_path, dpi=120)
            plt.close(fig)

            from reportlab.platypus import Image
            story.append(Image(img_path, width=400, height=250))
            story.append(Spacer(1, 12))

    story.append(Paragraph("Gut Microbiome Diversity", styles["Section"]))
    story.append(divider())
    story.append(Paragraph(diversity_text, styles["Body"]))
    # === NEW SECTION: Microbial Diversity Overview (Chuckling Goat–style) ===
    story.append(PageBreak())
    story.append(Paragraph("Microbial Diversity Overview", styles["Section"]))
    story.append(divider())

    div_score_text = f"""
        A healthy gut microbiome is defined by both richness (how many different microbes you have) "
        "and evenness (how evenly they are distributed). Together, they reflect how resilient and adaptable "
        "your gut ecosystem is. Below are your current diversity measures and their interpretations, "
        "presented on a 0–10 scale for easy reference. Based on your scores microbiome diversity, you have {diversity_state}
    """
    story.append(Paragraph(div_score_text, styles["Body"]))
    story.append(Spacer(1, 12))

    # --- Richness (Simplified: Shannon only, no observed count) ---
    story.append(Paragraph(f"<b>Richness (Shannon Index):</b> {shannon}/10", styles["Body"]))
    story.append(
        CleanGauge(shannon, width=450, height=14, max_value=10.0, mode="shannon",
                   label=f"{shannon}/10")
    )
    story.append(Spacer(1, 12))

    # --- Evenness Gauge ---
    evenness_value = simpson
    evenness_score = round(evenness_value * 10, 1)
    story.append(Paragraph(f"<b>Evenness (Simpson Index):</b> {evenness_score}/10", styles["Body"]))
    story.append(
        CleanGauge(evenness_score, width=450, height=14, max_value=10.0, mode="score", label=f"{evenness_score}/10")
    )
    story.append(Spacer(1, 20))

    # --- Diversity Score (Overall) ---
    diversity_score = round((shannon / 6) * 10, 1)
    diversity_score = min(10, max(0, diversity_score))
    story.append(Paragraph(f"<b>Microbial Diversity Score:</b> {diversity_score}/10", styles["Body"]))
    story.append(
        CleanGauge(diversity_score, width=450, height=16, max_value=10.0, mode="score",
                   label=f"{diversity_score}/10")
    )
    story.append(Spacer(1, 12))
    # --- Interpretation Text ---
    if diversity_score < 4:
        interp_text = (
            "Your gut diversity is currently below average, meaning a few bacterial groups dominate while others are underrepresented. "
            "This can reduce your gut’s flexibility and resilience. To support improvement, increase the variety of plant-based foods "
            "in your diet — aim for 20–30 different plants per week."
        )
    elif diversity_score < 7:
        interp_text = (
            "Your microbiome shows moderate diversity, indicating a reasonably balanced ecosystem. "
            "Adding more types of fibre and fermented foods can further strengthen this balance and resilience."
        )
    else:
        interp_text = (
            "Excellent diversity — your microbiome appears balanced and resilient, with strong representation of beneficial species. "
            "Continue maintaining a varied, fibre-rich diet to sustain this healthy ecosystem."
        )

    story.append(Paragraph(interp_text, styles["Body"]))
    story.append(Spacer(1, 20))

    story.append(Paragraph(
        "High microbial diversity is generally linked to improved digestive comfort, stronger immune function, and better overall wellbeing. "
        "Consistency and variety in fibre intake are the strongest ways to keep your gut ecosystem thriving.",
        styles["Body"]
    ))
    #story.append(PageBreak())

    #story.append(Spacer(1, 10))
    #story.append(CleanGauge(shannon, width=450, height=14, label=f"{shannon}/10"))
    story.append(PageBreak())
    story.append(Paragraph("Probiotics", styles["Section"]))
    story.append(divider())
    story.append(Paragraph(probiotic_section_text, styles["Body"]))

    # ========= Probiotic Section =========

    story.append(Paragraph(
        "Below are the probiotic or probiotic-like "
        "genera detected in your sample, along with their relative abundance levels.",
        styles["Body"]
    ))
    story.append(Spacer(1, 12))

    # --- Filter probiotics detected in this sample ---
    sample_probiotics = (
        sample_abundance[sample_abundance["Genus"].isin(probiotic_genera)]
        .copy()
        .groupby("Genus", as_index=False)["Abundance"]
        .sum()
        .sort_values(by="Abundance", ascending=False)
    )

    if sample_probiotics.empty:
        story.append(Paragraph(
            "No key probiotic genera (such as Lactobacillus or Bifidobacterium) were detected in this sample. "
            "This is common in adult microbiomes not recently exposed to fermented foods or probiotic supplements.",
            styles["Body"]
        ))
    else:
        max_abundance = max(sample_probiotics["Abundance"].max(), 5.0)
        for _, row in sample_probiotics.iterrows():
            genus = row["Genus"]
            genus_norm = normalize_genus(genus)
            raw_abundance = float(row["Abundance"])
            label = f"{raw_abundance:.3f}%" if raw_abundance < 1 else f"{raw_abundance:.1f}%"

            # --- Header for each genus
            story.append(Paragraph(f"<b><i>{genus}</i></b> — Estimated abundance: {label}", styles["Body"]))
            story.append(CleanGauge(
                raw_abundance, width=450, height=14, max_value=max_abundance,
                mode="abundance", label=label
            ))
            story.append(Spacer(1, 6))

            # --- Extract JSON-based description dynamically
            info = get_bacteria_info(genus_norm, raw_abundance)

            # --- Extract JSON-based description dynamically ---
            genus_norm = normalize_genus(genus)

            info = get_bacteria_info(
                genus=genus_norm,
                abundance=raw_abundance,
                thresholds=thresholds,
                kb_normalized=kb_normalized
            )


            if info:
                story.append(Paragraph(f"<b>Overview:</b> {info['overview']}", styles["Body"]))
                story.append(Paragraph(f"<b>Role:</b> {info['role']}", styles["Body"]))
                story.append(Paragraph(f"<b>Interpretation:</b> {info['level_text']}", styles["Body"]))
                if info["recommendations"]:
                    story.append(Paragraph("<b>How to support:</b>", styles["Body"]))
                    story.append(
                        Paragraph("<br/>".join([f"• {r}" for r in info["recommendations"]]), styles["Body"])
                    )
                story.append(Spacer(1, 10))
            else:
                story.append(Paragraph(
                    f"<font color='grey'><i>No detailed summary available yet for {genus}.</i></font>",
                    styles["Body"]
                ))
                story.append(Spacer(1, 10))


    story.append(Spacer(1, 12))


    story.append(Spacer(1, 12))


    # --- Probiotic Diversity Score Section ---
    if probiotic_diversity_score < 4:
        interp = "below average probiotic diversity"
    elif probiotic_diversity_score < 7:
        interp = "moderate probiotic diversity"
    else:
        interp = "excellent probiotic diversity"

    story.append(Paragraph("<b>Probiotic Diversity Score</b>", styles["Body"]))
    story.append(
        CleanGauge(probiotic_diversity_score, width=450, height=14, mode="score",
                   label=f"{probiotic_diversity_score}/10"))
    story.append(Spacer(1, 6))

    story.append(Paragraph(
        f"Your <b>Probiotic Diversity Score</b> is <b>{probiotic_diversity_score}/10</b>, indicating {interp}.<br/><br/>"
        "This score reflects the number and balance of beneficial probiotic species detected in your gut. "
        "A higher score suggests stronger microbial harmony and resilience, while a lower score means there’s room "
        "to nurture more beneficial bacteria.<br/><br/>"
        "To enhance probiotic diversity, include foods such as <i>kefir, yoghurt, sauerkraut, kimchi,</i> "
        "and fibre-rich prebiotics that feed these helpful microbes. Consistency over time helps probiotics establish and flourish.",
        styles["Body"]
    ))
    story.append(PageBreak())
    # --- Calculate Prebiotic Response Scores ---
    prebiotic_scores = {}
    max_ref_abundance = 10.0  # Normalization baseline (adjust as needed)

    for prebiotic, genera in prebiotic_targets.items():
        total_abundance = sample_abundance[sample_abundance["Genus"].isin(genera)]["Abundance"].sum()
        # Normalize abundance (log scale gives better balance)
        scaled_score = min(10, round((math.log1p(total_abundance) / math.log1p(max_ref_abundance)) * 10, 1))
        prebiotic_scores[prebiotic] = scaled_score

    # ========= Prebiotics Section =========
    story.append(Paragraph("Prebiotics", styles["Section"]))
    story.append(divider())

    story.append(Paragraph(prebiotic_section_text, styles["Body"]))

    story.append(Paragraph(
        "Prebiotics are types of dietary fibre that feed your beneficial gut bacteria. "
        "They help increase short-chain fatty acid production and encourage the growth of probiotic species like "
        "<i>Bifidobacterium</i> and <i>Faecalibacterium</i>. The more responsive your microbiome is to these fibres, "
        "the easier it becomes to maintain a healthy balance and regular digestion.",
        styles["Body"]
    ))
    story.append(Spacer(1, 12))

    prebiotic_data = (
        section_texts.get("sections", {}).get("Prebiotics", {
            "title": "Prebiotics",
            "text": "No Prebiotics section found in knowledge base; skipping this section."
        })
    )
    #story.append(Paragraph("<b>Prebiotics</b>", styles["Section"]))
    #story.append(divider())
    #story.append(Paragraph(prebiotic_data["intro"], styles["Body"]))
    story.append(Spacer(1, 10))

    # --- compute prebiotic score dynamically (you already had this part) ---
    prebiotic_score = compute_dynamic_score(sample_abundance,
                                            ["Faecalibacterium", "Roseburia", "Veillonella", "Bifidobacterium","Ruminococcus","Bacteroides","Lactobacillus","Eubacterium","Roseburia","Anaerostipes","Coprococcus","Prevotella","Phascolarctobacterium","Akkermansia",""],
                                            ideal_total=8)

    story.append(Paragraph(f"<b>Prebiotic Response Score:</b> {prebiotic_score:.1f}/10", styles["Body"]))
    story.append(CleanGauge(prebiotic_score, width=450, height=14, mode="score", label=f"{prebiotic_score:.1f}/10"))
    story.append(Spacer(1, 12))

    prebiotic_targets = ["Faecalibacterium", "Roseburia", "Subdoligranulum", "Bifidobacterium", "Akkermansia"]
    detected_prebiotics = [g for g in prebiotic_targets if g in sample_abundance["Genus"].unique()]

    if detected_prebiotics:
        detected_list = ", ".join([f"<i>{g}</i>" for g in detected_prebiotics])
        story.append(Paragraph(
            f"Your prebiotic response score reflects the presence of {detected_list}. "
            f"These beneficial bacteria are known to ferment dietary fibres into short-chain fatty acids, "
            f"helping support digestion, gut lining integrity, and immune balance.",
            styles["Body"]
        ))
    else:
        story.append(Paragraph(
            "No key prebiotic-associated bacteria were detected in this sample. "
            "This may indicate reduced fibre fermentation activity; increasing plant variety may help.",
            styles["Body"]
        ))


    #story.append(divider())

   # story.append(Spacer(1, 12))
    story.append(Paragraph(
        "High scores indicate your microbes efficiently utilize these prebiotic fibres, "
        "producing beneficial compounds like butyrate. Lower scores suggest your gut community may respond better "
        "to gradual fibre increases or more variety in plant-based foods.",
        styles["Body"]
    ))


    # ========= Postbiotics Section =========

    # ========= Postbiotics Section =========
    story.append(PageBreak())
    story.append(Paragraph("Postbiotics", styles["Section"]))
    story.append(divider())
    story.append(Paragraph(postbiotic_section_text, styles["Body"]))
    story.append(Paragraph(
        "Postbiotics are beneficial compounds such as short-chain fatty acids produced when your microbes digest fibre. "
        "They nourish your gut lining, reduce inflammation, and support energy metabolism. ",
        styles["Body"]
    ))
    story.append(Spacer(1, 10))

    postbiotic_targets = ["Faecalibacterium", "Roseburia", "Flavonifractor", "Coprococcus", "Akkermansia","Agathobacter","Anaerostipes","Butyricicoccus","Butyricimonas","Butyrivibrio","Pseudoflavonifractor","Lachnospira","Bifidobacterium","Blautia","Dorea","Lactobacillus","Prevotella","Ruminococcus","Alistipes","Veillonella","Sutterella","Megasphaera"]
    postbiotic_score = compute_dynamic_score(sample_abundance, postbiotic_targets, ideal_total=8)

    story.append(Paragraph(f"<b>Postbiotic Production Score:</b> {postbiotic_score:.1f}/10", styles["Body"]))
    story.append(CleanGauge(postbiotic_score, width=450, height=14, mode="score", label=f"{postbiotic_score:.1f}/10"))
    story.append(Spacer(1, 12))

    detected_postbiotics = [g for g in postbiotic_targets if g in sample_abundance["Genus"].unique()]
    if detected_postbiotics:
        detected_list = ", ".join([f"<i>{g}</i>" for g in detected_postbiotics])
        story.append(Paragraph(
            f"Your postbiotic production score is supported by {detected_list}. "
            f"These microbes help produce beneficial short-chain fatty acids such as butyrate and acetate, "
            f"which nourish gut cells and maintain microbial balance.",
            styles["Body"]
        ))
    else:
        story.append(Paragraph(
            "No major postbiotic-producing bacteria were detected. "
            "This could reflect low dietary fibre intake or recent antibiotic use.",
            styles["Body"]
        ))

    story.append(Spacer(1, 12))
    story.append(divider())

    story.append(Paragraph(
        "Your results show healthy postbiotic activity — especially in butyrate and propionate production. "
        "These compounds nourish intestinal cells and reduce inflammation. Maintaining fibre diversity and probiotic support "
        "will keep postbiotic output stable over time.",
        styles["Body"]
    ))
    story.append(PageBreak())

    # ========= Mental Wellness Section =========
    story.append(Paragraph("Mental Wellness", styles["Section"]))
    story.append(divider())

    story.append(Paragraph(
        "The gut and brain are deeply interconnected through the vagus nerve, immune signals, and microbial metabolites such as serotonin and GABA. "
        "Around 90% of serotonin is produced in the gut — meaning your microbiome plays a central role in emotional balance, focus, and stress resilience. "
        "Gut bacteria associated with mental health include a decrease in beneficial, anti-inflammatory bacteria like Faecalibacterium and Coprococcus in "
        "conditions like depression and anxiety, and an increase in pro-inflammatory bacteria like Eggerthella. The gut microbiome influences mental health "
        "through pathways like the vagus nerve, immune system, and the production of mood-regulating chemicals like serotonin. Probiotics containing species "
        "like Lactobacillus and Bifidobacterium have shown promise for improving mood in some studies.The brain has "
        "a direct effect on the stomach and intestines. For example, the very thought of eating can release the stomach's juices before food gets there. "
        "This connection goes both ways. A troubled intestine can send signals to the brain, just as a troubled brain can send signals to the gut. Therefore, "
        "a person's stomach or intestinal distress can be the cause or the product of anxiety, stress, or depression. That's because the brain and the "
        "gastrointestinal (GI) system are intimately connected. This is especially true in cases where a person experiences gastrointestinal upset with no "
        "obvious physical cause. For such functional GI disorders, it is difficult to try to heal a distressed gut without considering the role of stress and "
        "emotion.",
        styles["Body"]
    ))


    story.append(Spacer(1, 10))

    mental_bacteria = ["Bifidobacterium", "Lactobacillus", "Faecalibacterium", "Coprococcus"]
    mental_score = compute_dynamic_score(sample_abundance, mental_bacteria, ideal_total=7)
    detected_mental_bacteria = [g for g in mental_bacteria if g in sample_abundance["Genus"].unique()]

    if mental_score < 4:
        mental_text = f"Your mental health score is {mental_score}. Low representation of mood-supporting bacteria. Add more fermented foods and fibre to improve your gut–brain axis."
    elif mental_score < 7:
        mental_text = f"Your mental health score is {mental_score}. Moderate balance — continue supporting it with probiotics and stress management."
    else:
        mental_text = f"Your mental health score is {mental_score}. Excellent representation of mood-supporting bacteria — your gut–brain axis is thriving."

    story.append(Paragraph(mental_text, styles["Body"]))
    story.append(Spacer(1, 12))
    story.append(Paragraph(f"<b>Gut–Brain Axis Score:</b> {mental_score}/10", styles["Body"]))
    story.append(CleanGauge(mental_score, width=450, height=16, mode="score", label=f"{mental_score}/10"))
    story.append(Spacer(1, 12))
    if detected_mental_bacteria:
        detected_list = ", ".join([f"<i>{g}</i>" for g in detected_mental_bacteria])
        story.append(Paragraph(
            f"Your mental health score reflects the presence of {detected_list}. ",
            styles["Body"]
        ))
    else:
        story.append(Paragraph(
            "No key prebiotic-associated bacteria were detected in this sample. "
            "This may indicate reduced fibre fermentation activity; increasing plant variety may help.",
            styles["Body"]
        ))
    story.append(Spacer(1, 12))
    story.append(Paragraph(
        "Faecalibacterium is a beneficial gut bacterium that plays a role in mental health by producing the short-chain "
        "fatty acid (SCFA) butyrate, which supports brain health. A higher abundance of Faecalibacterium is linked to reduced "
        "depressive symptoms and increased mental resilience, potentially by regulating inflammation and supporting neurogenesis"
        " in the hippocampus. Conversely, lower levels of this bacterium have been observed in individuals with depression and "
        "other mood disorders.  ",
        styles["Body"]
    ))
    story.append(Spacer(1, 12))
    story.append(Paragraph(
        "Coprococcus is a gut bacterium associated with mental health due to its role in producing anti-inflammatory "
        "compounds like butyrate, which can strengthen the gut lining and have beneficial effects on the brain via the "
        "gut-brain axis. Studies show that people with mental health conditions like depression often have lower levels "
        "of Coprococcus. The decrease in this beneficial bacteria, and the associated reduction in butyrate and other "
        "beneficial metabolites, may contribute to inflammation and other factors that negatively impact mental well-being. ",
        styles["Body"]
    ))
    story.append(PageBreak())

    # ========= Skin Health Section =========

    story.append(Paragraph("Skin Health", styles["Section"]))
    story.append(divider())

    story.append(Paragraph(
        "Your gut and skin are connected through immune and metabolic pathways. "
        "Balanced microbial diversity can help reduce dryness, redness, and inflammation. "
        "The microbiome plays an important role in a wide variety of skin disorders. Not only "
        "is the skin microbiome altered, but also surprisingly many skin diseases are accompanied "
        "by an altered gut microbiome. The microbiome is a key regulator for the immune system, as "
        "it aims to maintain homeostasis by communicating with tissues and organs in a bidirectional "
        "manner. Hence, dysbiosis in the skin and/or gut microbiome is associated with an altered "
        "immune response, promoting the development of skin diseases, such as atopic dermatitis, "
        "psoriasis, acne vulgaris, dandruff, and even skin cancer. Here, we focus on the associations "
        "between the microbiome, diet, metabolites, and immune responses in skin pathologies. "
        "This review describes an exhaustive list of common skin conditions with associated dysbiosis "
        "in the skin microbiome as well as the current body of evidence on gut microbiome dysbiosis, "
        "dietary links, and their interplay with skin conditions. An enhanced understanding of the local "
        "skin and gut microbiome including the underlying mechanisms is necessary to shed light on the "
        "microbial involvement in human skin diseases and to develop new therapeutic approaches.",
        styles["Body"]
    ))
    story.append(Spacer(1, 10))

    skin_bacteria = ["Akkermansia", "Bifidobacterium", "Faecalibacterium", "Roseburia"]
    skin_score = compute_dynamic_score(sample_abundance, skin_bacteria, ideal_total=7)
    detected_skin_bacteria = [g for g in skin_bacteria if g in sample_abundance["Genus"].unique()]

    if skin_score < 4:
        skin_text = f"Your skin health score is {skin_score}. Low levels of skin-supportive bacteria — increase prebiotics and polyphenol-rich foods."
    elif skin_score < 7:
        skin_text = f"Your skin health score is {skin_score}. Moderate balance — your gut–skin axis is stable. Keep supporting it with hydration and antioxidants."
    else:
        skin_text = f"Your skin health score is {skin_score}. Excellent balance — your gut–skin axis appears healthy and resilient."

    story.append(Paragraph(skin_text, styles["Body"]))
    story.append(Spacer(1, 12))
    story.append(Paragraph(f"<b>Gut–Skin Axis Score:</b> {skin_score}/10", styles["Body"]))
    story.append(CleanGauge(skin_score, width=450, height=16, mode="score", label=f"{skin_score}/10"))
    story.append(Spacer(1, 12))
    if detected_skin_bacteria:
        detected_list = ", ".join([f"<i>{g}</i>" for g in detected_skin_bacteria])
        story.append(Paragraph(
            f"Your skin health score reflects the presence of {detected_list}. "
            f"These beneficial bacteria are known to ferment dietary fibres into short-chain fatty acids, "
            f"helping support digestion, gut lining integrity, and immune balance. "
            f"Faecalibacterium prausnitzii, a key gut bacterium, influences skin health primarily through its anti-inflammatory effects and by "
            f"reinforcing the intestinal barrier, a connection known as the gut-skin axis. When the population of "
            f"F. prausnitzii is low, systemic inflammation can increase, potentially contributing to inflammatory "
            f"skin conditions like atopic dermatitis (eczema) and psoriasis. ",
            styles["Body"]
        ))
    else:
        story.append(Paragraph(
            "No key prebiotic-associated bacteria were detected in this sample. "
            "This may indicate reduced fibre fermentation activity; increasing plant variety may help.",
            styles["Body"]
        ))

    story.append(PageBreak())

    # ========= Pathogens Section =========

    story.append(Paragraph("Potential Pathogens", styles["Section"]))
    story.append(divider())

    story.append(Paragraph(
        "While most gut bacteria are helpful, some can become problematic if they grow too much. "
        "This section shows potentially pathogenic genera found in your sample.",
        styles["Body"]
    ))
    story.append(Spacer(1, 10))

    pathogen_section_text = """
    <b>Understanding Potential Pathogens & Opportunists</b><br/><br/>

    Your gut microbiome is a diverse ecosystem made up of both helpful and opportunistic microbes. 
    While most bacteria play beneficial roles in digestion and immune regulation, 
    a few species are considered <b>opportunistic</b> — meaning they can become troublesome 
    when they grow in excess or when the balance of beneficial bacteria is disrupted.<br/><br/>

    <b>1. What are opportunistic bacteria?</b><br/><br/>
    These are bacteria that normally exist in small, harmless amounts within a healthy gut. 
    However, under certain conditions — such as stress, illness, antibiotic use, or dietary imbalance — 
    they can multiply and start to produce by-products that irritate the gut lining, disturb digestion, 
    or compete with beneficial microbes for nutrients. Examples include certain species of 
    <i>Escherichia</i>, <i>Enterococcus</i>, <i>Streptococcus</i>, <i>Clostridium</i>, and <i>Campylobacter</i>.<br/><br/>

    <b>2. Why balance matters</b><br/><br/>
    In a balanced gut ecosystem, beneficial bacteria such as <i>Bifidobacterium</i> and <i>Faecalibacterium</i> 
    keep these opportunists under control. They do this by producing short-chain fatty acids (like butyrate) 
    that lower the gut’s pH, making it less hospitable for pathogens, and by directly competing for resources and space.<br/><br/>
    If this balance shifts — for example, after antibiotics or a highly processed diet — 
    these opportunistic microbes can temporarily gain an advantage. 
    This doesn’t necessarily mean infection or disease, but rather that your gut’s natural equilibrium 
    needs gentle support to return to harmony.<br/><br/>

    <b>3. Common reasons opportunistic bacteria rise</b><br/>
    • Recent antibiotic or medication use<br/>
    • High intake of refined carbohydrates or processed foods<br/>
    • Chronic stress, which alters gut motility and immune signaling<br/>
    • Low fibre intake, reducing “good bacteria” food sources<br/>
    • Poor sleep or travel, which can disrupt the microbiome rhythm<br/><br/>

    <b>4. How we assess them in your sample</b><br/><br/>
    We identify potential pathogens and opportunists through genetic sequencing (16S rRNA analysis). 
    The report shows their relative abundance as a percentage of your total microbial population. 
    Small traces are normal and expected — everyone has some. 
    Only significant overgrowths are flagged as “elevated” or “needs support.”<br/><br/>

    <b>5. Supporting microbial balance naturally</b><br/><br/>
    If your results show higher-than-expected levels of opportunistic bacteria, 
    the focus is on <b>rebalance, not eradication</b>. 
    You can encourage beneficial species and restore harmony through various methods as advised by your consultant<br/><br/>

    <b>6. Interpreting your results</b><br/><br/>
    Low levels of opportunistic bacteria are completely normal and part of a healthy gut ecosystem. 
    If one or more are moderately elevated, it may simply reflect recent diet, travel, or medication changes. 
    Only when several opportunistic species are consistently high does it indicate a more sustained imbalance 
    (often called <i>dysbiosis</i>).<br/><br/>

    <b>In summary</b><br/><br/>
    Your gut is designed to be self-regulating. When supported with the right nutrition, hydration, and stress management, 
    beneficial microbes naturally regain control and keep opportunistic species in check. 
    Your “Pathogen Balance Score” shows how well this equilibrium is currently maintained — 
    a higher score means a more resilient and well-regulated microbiome.
    """

    pathogen_targets = ["Escherichia", "Clostridium", "Campylobacter", "Enterococcus", "Streptococcus","Citrobacter","Enterobacter","Pseudomonas","Salmonella","Sphingomonas","Vibrio","Yersinia","Staphylococcus","Eggerthella","Fusobacterium","Clostridioides"]
    pathogen_data = (
        sample_abundance[sample_abundance["Genus"].isin(pathogen_targets)]
        .groupby("Genus", as_index=False)["Abundance"]
        .sum()
        .sort_values(by="Abundance", ascending=False)
    )
    pathogens_score = compute_dynamic_score(sample_abundance, pathogen_targets, ideal_total=2, invert=True)

    if pathogen_data.empty:
        story.append(Paragraph(
            "No significant pathogens were detected — your gut shows excellent microbial balance.",
            styles["Body"]
        ))
    else:
        max_abundance = max(pathogen_data["Abundance"].max(), 1.0)
        for _, row in pathogen_data.iterrows():
            genus, abundance = row["Genus"], float(row["Abundance"])
            label = f"{abundance:.3f}%" if abundance < 1 else f"{abundance:.1f}%"
            story.append(Paragraph(f"<b><i>{genus}</i></b> — {label}", styles["Body"]))


        #if not pathogen_data.empty:
          #  detected_pathogens = ", ".join([f"<i>{g}</i>" for g in pathogen_data["Genus"]])
          #  story.append(Paragraph(
          #      f"In your sample, the following opportunistic bacteria were detected in measurable amounts: {detected_pathogens}. "
         #       f"All are common members of the gut ecosystem and typically pose no issue unless they grow in excess.",
         #       styles["Body"]
         #   ))
        #else:
         #   story.append(Paragraph(
         #       "No notable opportunistic or pathogenic bacteria were detected — your gut shows excellent microbial resilience.",
         #       styles["Body"]
         #   ))
        #story.append(
        #CleanGauge(abundance, width=450, height=14, max_value=max_abundance, mode="abundance", label=label))
    story.append(Spacer(1, 10))
    story.append(Spacer(1, 6))
    story.append(Spacer(1, 10))
    story.append(Paragraph(f"<b>Pathogen Balance Score:</b> {pathogens_score}/10", styles["Body"]))
    story.append(Spacer(1, 10))
    story.append(CleanGauge(pathogens_score, width=450, height=14, mode="score", label=f"{pathogens_score}/10"))
    story.append(Spacer(1, 10))
    story.append(Paragraph(pathogen_section_text, styles["Body"]))
    story.append(PageBreak())

    good_bugs = [
        "Bifidobacterium breve", "Bifidobacterium bifidum", "Eubacterium rectale",
        "Butyricicoccus pullicaecorum", "Anaerostipes hadrus", "Anaerobutyricum hallii",
        "Coprococcus comes", "Subdoligranulum", "Oscillibacter", "Phascolarctobacterium faecium",
        "Barnesiella intestinihominis", "Christensenella minuta", "Lactobacillus", "Lactococcus lactis",
        "Streptococcus thermophilus", "Butyricimonas", "Oxalobacter formigenes", "Methanobrevibacter smithii",
        "Faecalibacterium","Roseburia","Bifidobacterium"
    ]

    bad_bugs = [
        "Klebsiella pneumoniae", "Escherichia coli", "Proteus mirabilis",
        "Citrobacter freundii", "Salmonella enterica", "Campylobacter jejuni",
        "Helicobacter pylori", "Fusobacterium nucleatum", "Pseudomonas aeruginosa",
        "Staphylococcus aureus", "Candida albicans", "Candida glabrata",
        "Blastocystis hominis", "Giardia lamblia", "Cryptosporidium",
        "Entamoeba histolytica", "Bilophila wadsworthia", "Desulfovibrio",
        "Shigella", "ARGs", "Biofilm","Clostridioides","Enterococcus", "Morganella", "Klebsiella oxytoca","Klebsiella"
    ]

    variable_bugs = [
        "Bacteroides fragilis", "Bacteroides ovatus", "Bacteroides thetaiotaomicron",
        "Bacteroides stercoris", "Bacteroides uniformis", "Alistipes putredinis",
        "Prevotella copri", "Parabacteroides distasonis", "Parabacteroides merdae",
        "Blautia obeum", "Blautia wexlerae", "Eggerthella lenta",
        "Enterocloster bolteae", "Enterocloster clostridioformis", "Veillonella",
        "Megasphaera", "Sutterella", "Dorea longicatena", "Romboutsia timonensis",
        "Intestinibacter bartlettii", "Intestinimonas massiliensis", "Flavonifractor plautii",
        "Gordonibacter pamelaeae", "Hungatella hathewayi", "Porphyromonas", "Ruminococcus", "Collinsella"
    ]

    def detect_biomarkers(sample_abundance, sample_id_norm=None):
        """
        Detects presence and balance of good, variable, and bad gut bacteria.
        Accepts either a full dataset or a pre-filtered subset.
        """
        if "Sample_norm" in sample_abundance.columns and sample_id_norm:
            sample_data = sample_abundance[sample_abundance["Sample_norm"] == sample_id_norm]
        else:
            sample_data = sample_abundance  # already filtered or grouped

        def get_abundance(genera_list):
            return sample_data[sample_data["Genus"].isin(genera_list)]["Abundance"].sum()

        total_good = get_abundance([g.split()[0] for g in good_bugs])
        print(total_good)
        print(good_bugs)
        total_bad = get_abundance([g.split()[0] for g in bad_bugs])
        total_variable = get_abundance([g.split()[0] for g in variable_bugs])

        # Compute relative ratios for easy visualization
        total_all = total_good + total_bad + total_variable
        if total_all == 0:
            total_all = 1  # avoid zero division

        good_score = round((total_good / total_all) * 10, 1)
        bad_score = round((10 - ((total_bad / total_all) * 10)), 1)
        variable_score = round((total_variable / total_all) * 10, 1)
        gut_balance_index = round((good_score + bad_score) / 2, 1)

        # Return both detailed and summary data
        summary = {
            "Good_Bugs (%)": total_good,
            "Variable_Bugs (%)": total_variable,
            "Bad_Bugs (%)": total_bad,
            "Good_Score": good_score,
            "Bad_Score": bad_score,
            "Variable_Score": variable_score,
            "Gut_Balance_Index": gut_balance_index
        }

        return summary

    # 🧹 Step 1: De-duplicate genus entries by summing their abundances
    sample_abundance_grouped = (
        sample_abundance.groupby("Genus", as_index=False)["Abundance"]
        .sum()
        .sort_values(by="Abundance", ascending=False)
    )

    story.append(Paragraph("Microbial Biomarker Summary", styles["Section"]))
    story.append(divider())

    # 🧮 Compute biomarker balance summary using grouped data
    bio_summary = detect_biomarkers(sample_abundance_grouped)

    # 🧾 Show core balance data
  #  story.append(Paragraph(
  #      f"<b>Good Bugs:</b> {bio_summary['Good_Bugs (%)']:.2f}%<br/>"
  #      f"<b>Variable Bugs:</b> {bio_summary['Variable_Bugs (%)']:.2f}%<br/>"
  #      f"<b>Bad Bugs:</b> {bio_summary['Bad_Bugs (%)']:.2f}%",
  #      styles["Body"]
  #  ))

    story.append(Spacer(1, 10))

    # 🧭 Build Microbiome Resilience Radar image
    import os
    radar_path = os.path.join("/tmp", f"radar_{sample_id_norm}.png")
    build_microbiome_resilience_radar(bio_summary, radar_path)
    # === Microbial Biomarker Pie Chart ===
    story.append(Paragraph(
        "This chart summarises your gut bacterial balance, showing the proportion of beneficial (Good), opportunistic (Bad), and variable (context-dependent) microbes detected in your sample.",
        styles["Body"]
    ))
    story.append(Spacer(1, 12))

    # Count total abundance in each group
    good_bugs_total = float(sample_abundance_grouped[sample_abundance_grouped["Genus"].isin(
        ["Faecalibacterium", "Coprococcus", "Roseburia", "Bifidobacterium", "Lactobacillus",
         "Akkermansia", "Subdoligranulum", "Ruminococcus"]
    )]["Abundance"].sum())

    bad_bugs_total = float(sample_abundance_grouped[sample_abundance_grouped["Genus"].isin(
        ["Campylobacter", "Escherichia/Shigella", "Enterococcus", "Clostridium", "Enterobacter", "Desulfovibrio"]
    )]["Abundance"].sum())

    variable_bugs_total = float(sample_abundance_grouped[sample_abundance_grouped["Genus"].isin(
        ["Prevotella", "Sutterella", "Dorea", "Collinsella", "Parabacteroides", "Alistipes"]
    )]["Abundance"].sum())

    # Add chart
    story.extend(filter(None, [
        make_biomarker_pie(good_bugs_total, bad_bugs_total, variable_bugs_total),
    ]))
    story.append(Spacer(1, 20))

    from reportlab.platypus import Image
    story.append(Spacer(1, 6))
    story.append(Image(radar_path, width=300, height=300))
    story.append(Spacer(1, 12))


    # 🧠 Interpret the balance scores
    gbi = bio_summary["Gut_Balance_Index"]
    good = bio_summary["Good_Bugs (%)"]
    bad = bio_summary["Bad_Bugs (%)"]
    variable = bio_summary["Variable_Bugs (%)"]

    if gbi >= 8:
        balance_text = (
            "Your gut microbiome shows <b>excellent resilience and stability</b>. "
            "Beneficial microbes dominate the ecosystem, helping to maintain optimal digestion, "
            "metabolic health, and immune modulation. Opportunistic species are well-controlled, "
            "suggesting strong microbial diversity and balance."
        )
    elif gbi >= 6:
        balance_text = (
            "Your gut balance is <b>moderately healthy</b>. Beneficial bacteria are present in good numbers, "
            "though there are some variable or opportunistic strains that could grow under stress, poor diet, or antibiotics. "
            "Maintaining fibre diversity and probiotic foods will help keep your microbiome stable."
        )
    else:
        balance_text = (
            "Your microbiome shows <b>reduced balance</b>, with higher levels of variable or opportunistic microbes. "
            "This pattern may indicate inflammation, low fibre intake, or recent antibiotic exposure. "
            "Introducing prebiotic and probiotic support may help rebalance your gut ecosystem."
        )

    story.append(Paragraph(balance_text, styles["Body"]))
    story.append(Spacer(1, 12))
    # 🧬 Detailed tables by category
    story.append(Spacer(1, 12))
    story.append(Paragraph("<b>Detailed Microbial Overview</b>", styles["Section"]))
    story.append(divider())

    story.extend(make_biomarker_table(
        sample_abundance_grouped,
        "Detected Good Bacteria",
        styles,
        genus_subset=["Faecalibacterium", "Coprococcus", "Roseburia", "Bifidobacterium",
                      "Lactobacillus", "Akkermansia", "Subdoligranulum", "Ruminococcus"],
        kb_normalized=kb_normalized,
        thresholds=thresholds
    ))

    story.extend(make_biomarker_table(
        sample_abundance_grouped,
        "Detected Opportunistic Bacteria",
        styles,
        genus_subset=["Campylobacter", "Escherichia/Shigella", "Enterococcus",
                      "Clostridium", "Enterobacter", "Desulfovibrio"],
        kb_normalized=kb_normalized,
        thresholds=thresholds
    ))

    story.extend(make_biomarker_table(
        sample_abundance_grouped,
        "Detected Variable Bacteria",
        styles,
        genus_subset=["Prevotella", "Sutterella", "Dorea", "Collinsella", "Parabacteroides", "Alistipes"],
        kb_normalized=kb_normalized,
        thresholds=thresholds
    ))
    story.append(Spacer(1, 6))
    story.append(Paragraph("<b>Gut Balance Index</b>", styles["Body"]))
    story.append(
        CleanGauge(
            bio_summary["Gut_Balance_Index"],
            width=450,
            height=16,
            label=f"{bio_summary['Gut_Balance_Index']}/10"
        )
    )
    story.append(Spacer(1, 12))
    story.append(Paragraph(
        "This Gut Balance Index reflects the ratio of beneficial to disruptive microbial species. "
        "Higher values indicate a microbiome that resists imbalance, supports SCFA production, and contributes to lower inflammation risk.",
        styles["Body"]
    ))
    story.append(Spacer(1, 12))

    # 🧬 Show detected bacteria categories (grouped, unique genera)
    def list_detected_biomarkers(abundance_df, bug_list, label):
        """Lists detected microbes from good/bad/variable categories (unique genera only)."""
        detected = abundance_df[abundance_df["Genus"].isin([g.split()[0] for g in bug_list])]
        if detected.empty:
            return f"No {label.lower()} bacteria detected."
        lines = [f"<i>{row['Genus']}</i> — {row['Abundance']:.2f}%" for _, row in detected.iterrows()]
        return "<br/>".join(lines)

    #story.append(Spacer(1, 8))
    #story.append(Paragraph("<b>Detected Good Bacteria:</b>", styles["Body"]))
    #story.append(Paragraph(list_detected_biomarkers(sample_abundance_grouped, good_bugs, "Good"), styles["Body"]))
    #story.append(Spacer(1, 8))

   # story.append(Paragraph("<b>Detected Opportunistic Bacteria:</b>", styles["Body"]))
    #story.append(Paragraph(list_detected_biomarkers(sample_abundance_grouped, bad_bugs, "Bad"), styles["Body"]))

    #story.append(Spacer(1, 8))
    #story.append(PageBreak())
    #story.append(Paragraph("<b>Detected Variable Bacteria:</b>", styles["Body"]))
    #story.append(Paragraph(list_detected_biomarkers(sample_abundance_grouped, variable_bugs, "Bad"), styles["Body"]))
    story.append(PageBreak())
    # ========= FINAL PDF BUILD =========
    # ========= FINAL PDF BUILD (Stable & Quiet) =========
    import os
    from reportlab.platypus import SimpleDocTemplate

    # --- Output path ---

    pdf_path = os.path.join(
        os.path.expanduser("~/Downloads"),
        f"Your_Poopie_Report_{sample_id.replace('.', '_')}.pdf"
    )

    from datetime import datetime

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")


    # ensure output directory exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    pdf_name = f"{timestamp}_{sample_id}.pdf"
    pdf_path = os.path.join(os.path.dirname(args.output), pdf_name)
    # --- Also create a .txt file with JSON output ---
    txt_name = f"{timestamp}_{sample_id}.txt"
    txt_path = os.path.join(os.path.dirname(args.output), txt_name)

    try:
        # Pretty-print JSON data to text file
        with open(txt_path, "w") as f:
            json.dump(data_all, f, indent=2)
        print(f"📝 TXT report saved to: {txt_path}")
    except Exception as e:
        print(f"❌ Failed to write TXT report: {e}")

    print(f"📝 Writing PDF to: {pdf_path}")

    # --- Page numbering only after Introduction ---
    def on_page_numbered(canvas, doc_obj):
        """Draw footer page numbers only from Introduction onward."""
        nonlocal intro_page_number  # refers to outer variable that marks where Introduction starts

        # Skip page numbering before the Introduction section
        if intro_page_number is None:
            return

        # Get the current page number from ReportLab's canvas object
        page_num = canvas.getPageNumber()

        # Only start showing numbers once we're past the intro start
        if page_num >= intro_page_number:
            display_num = page_num - intro_page_number + 1
            canvas.setFont("Helvetica", 9)
            canvas.setFillColor(colors.grey)
            canvas.drawCentredString(A4[0] / 2, 1.2 * cm, f"Page {display_num}")

    # --- Create PDF document ---
    doc = SimpleDocTemplate(
        pdf_path,
        pagesize=A4,
        rightMargin=2 * cm,
        leftMargin=2 * cm,
        topMargin=2 * cm,
        bottomMargin=1.8 * cm,
    )
    doc.afterFlowable = after_flowable
    # === just before PDF build ===
    print(f"✅ Unique genera detected: {len(sample_abundance['Genus'].unique())}")
    print(sample_abundance.head(10))

    print("⚙️ Starting PDF build ...")

    # ========= Reorder Pages =========
    # Move the "Microbial Biomarker Summary" section right after the Introduction

    def find_section_index(title):
        """Helper to find a section by title text"""
        for i, elem in enumerate(story):
            # Check if the flowable has text and matches title
            if hasattr(elem, "text") and title.lower() in elem.text.lower():
                return i
        return None

    intro_idx = find_section_index("Introduction")
    biomarker_idx = find_section_index("Microbial Biomarker Summary")

    if intro_idx is not None and biomarker_idx is not None:
        # Extract everything from "Microbial Biomarker Summary" to the next PageBreak (or end)
        end_idx = biomarker_idx + 1
        for j in range(biomarker_idx + 1, len(story)):
            if isinstance(story[j], PageBreak):
                end_idx = j + 1
                break

        # Slice out the section
        biomarker_section = story[biomarker_idx:end_idx]

        # Delete original section
        del story[biomarker_idx:end_idx]

        # Insert right after introduction
        insert_pos = intro_idx + 4  # after intro + divider
        story[insert_pos:insert_pos] = biomarker_section

        print("✅ Moved 'Microbial Biomarker Summary' section after Introduction.")
    else:
        print("⚠️ Could not find 'Microbial Biomarker Summary' or 'Introduction' section.")

    try:
        doc.multiBuild(story, onFirstPage=on_page_numbered, onLaterPages=on_page_numbered)

        print("✅ PDF generated successfully!")

        print(f"📄 Saved to: {pdf_path}")
    except Exception as e:
        print(f"❌ PDF build error: {e}")




# ========= Script Entry Point =========
if __name__ == "__main__":

    generate_report(args.sample)
