import os
import re
import gzip
import urllib.request
import time
import sys
import subprocess

# =========================
# AUTO-INSTALL DEPENDENCIES
# =========================
def ensure_package(pip_name, import_name=None):
    import_name = import_name or pip_name
    try:
        return __import__(import_name)
    except ImportError:
        print(f"[AUTO-INSTALL] Installing missing package: {pip_name}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pip_name])
        return __import__(import_name)

matplotlib = ensure_package("matplotlib")
import matplotlib.pyplot as plt
venn2 = ensure_package("matplotlib-venn", "matplotlib_venn").venn2

# =========================
# CONFIG (ROBUST PATHS)
# =========================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_DIR   = os.path.join(BASE_DIR, "data")
RAW_DIR    = os.path.join(DATA_DIR, "raw")
PROC_DIR   = os.path.join(DATA_DIR, "processed")

GENCODE_GTF    = os.path.join(RAW_DIR, "gencode.v49.chr.gtf.gz")
MIRBASE_FA     = os.path.join(RAW_DIR, "mature.fa")
MIRGENEDB_FILE = os.path.join(RAW_DIR, "mirgenedb.bed")

GENCODE_URL  = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.chr.gtf.gz"
MIRBASE_URL  = "https://www.mirbase.org/download/mature.fa"
MIRGENEDB_URL = "https://mirgenedb.org/static/data/hsa/hsa-all.bed"

# =========================
# UTILS
# =========================
def ensure_dirs():
    os.makedirs(RAW_DIR,  exist_ok=True)
    os.makedirs(PROC_DIR, exist_ok=True)
    print("[DEBUG] Working directory:", os.getcwd())

def debug_paths():
    print("\n[DEBUG PATHS]")
    print("  BASE_DIR     :", BASE_DIR)
    print("  RAW_DIR      :", RAW_DIR)
    print("  PROC_DIR     :", PROC_DIR)
    print("  GENCODE_GTF  :", GENCODE_GTF)
    print("  MIRBASE_FA   :", MIRBASE_FA)
    print("  MIRGENEDB_BED:", MIRGENEDB_FILE)
    print("-" * 50)

def download_file(url, output):
    if os.path.exists(output):
        print(f"[SKIP] Already exists: {output}")
        return
    print(f"[DOWNLOAD] {url}")
    try:
        urllib.request.urlretrieve(url, output)
        print(f"[DONE] Saved to: {output}")
    except Exception as e:
        print(f"[ERROR] Download failed: {e}")

def check_file(path, label=""):
    tag = f" ({label})" if label else ""
    if not os.path.exists(path):
        print(f"[ERROR] Missing file{tag}: {path}")
        return False
    size_mb = os.path.getsize(path) / (1024 * 1024)
    print(f"[OK] {path}{tag}  ({size_mb:.2f} MB)")
    return True

# =========================
# NORMALIZATION  (shared by all parsers)
# =========================
def normalize(name):
    """
    Converts any miRNA identifier to uppercase gene-level format.

    Handles all three database naming styles:
        GENCODE  : mir-21          →  MIR21
        miRBase  : hsa-miR-21-5p  →  MIR21
        miRGeneDB: Hsa-Mir-21_pre →  MIR21

    Key transforms applied (in order):
        1. Remove species prefix (hsa- / Hsa-)
        2. Normalise family keyword (Mir- / miR- / mir- → mir | Let- / let- → let)
        3. Remove arm annotation (-5p / -3p / _5p / _3p with optional *)
        4. Remove the hyphen between keyword and number (mir-21 → mir21)
        5. Uppercase everything
    """
    name = name.strip()

    # 1. Strip species prefix
    name = re.sub(r'^[Hh]sa-', '', name)

    # 2. Normalise family keyword to lowercase stub without trailing hyphen
    #    Handles: Mir-  miR-  mir-  →  mir
    #             Let-  let-        →  let
    name = re.sub(r'^[Mm]ir-', 'mir', name)
    name = re.sub(r'^[Ll]et-', 'let', name)

    # 3. Remove arm annotations wherever they appear
    #    Covers: -5p  -3p  _5p  _3p  with optional trailing *
    name = re.sub(r'[-_][35]p\*?$', '', name)

    # 4. Collapse the separator between keyword and number
    #    mir-21 → mir21   let-7g → let7g
    name = re.sub(r'^(mir|let)-', r'\1', name)

    # 5. Uppercase
    name = name.upper()
    return name

def collapse_family(name):
    """
    Collapse miRGeneDB paralog/family suffixes after normalize() has been applied.

    miRGeneDB encodes paralogs as  Hsa-Mir-8-P1a_pre, Hsa-Mir-8-P2a_pre …
    After normalize() those become  MIR8-P1A, MIR8-P2A …
    This function strips the -P<digits><letters> suffix so all paralogs
    map to the same gene-level token (MIR8).

    Also handles the LET-7 family style from miRGeneDB:
        LET-7-P2A1 → LET7   (the LET-7 hyphen is already gone after normalize,
                              but the -P2A1 suffix remains)
    """
    # Remove -P<number><optional letters> suffix (case-insensitive due to prior upper())
    name = re.sub(r'-P\d+[A-Z]*$', '', name)
    return name

# =========================
# STEP 1 — DATA ACQUISITION
# =========================
def acquire_data():
    print("\n" + "=" * 60)
    print("[STEP 1] Data acquisition")
    print("=" * 60)
    print("[STORY] The pipeline integrates three complementary databases:")
    print("  • GENCODE   — genomic loci and gene-level annotation (GTF)")
    print("  • miRBase   — mature miRNA sequences, including 5p/3p arms (FASTA)")
    print("  • miRGeneDB — evolutionarily validated, high-confidence miRNA genes (BED)")
    print("[STORY] Each database captures a different biological layer.")
    print("[STORY] Only by integrating all three can we define a truly robust miRNA set.")
    time.sleep(2)

    download_file(GENCODE_URL,  GENCODE_GTF)
    download_file(MIRBASE_URL,  MIRBASE_FA)
    download_file(MIRGENEDB_URL, MIRGENEDB_FILE)

    print("\n[CHECKPOINT] Validating downloaded files...")
    ok = all([
        check_file(GENCODE_GTF,    "GENCODE GTF"),
        check_file(MIRBASE_FA,     "miRBase FASTA"),
        check_file(MIRGENEDB_FILE, "miRGeneDB BED"),
    ])
    if not ok:
        print("[FATAL] One or more required files are missing. Aborting.")
        sys.exit(1)

# =========================
# STEP 2 — PARSE GENCODE
# =========================
def parse_gencode():
    print("\n" + "=" * 60)
    print("[STEP 2] Parsing GENCODE annotation")
    print("=" * 60)
    print("[STORY] GENCODE encodes miRNAs as genomic loci with gene/transcript structure.")
    print("[STORY] We filter entries where gene_type = 'miRNA' and collapse to unique gene symbols.")
    print("[INFO]  Fields extracted: CHROM, START, END, STRAND, GENE_NAME, GENE_ID")
    time.sleep(2)

    output_path = os.path.join(PROC_DIR, "gencode_mirna.txt")
    mirnas = set()
    example_records = []

    with gzip.open(GENCODE_GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if 'gene_type "miRNA"' not in line:
                continue

            cols = line.strip().split("\t")
            chrom  = cols[0]
            start  = cols[3]
            end    = cols[4]
            strand = cols[6]

            gn_match  = re.search(r'gene_name "([^"]+)"', line)
            gid_match = re.search(r'gene_id "([^"]+)"',   line)

            if not gn_match:
                continue

            gene_name = gn_match.group(1)
            mirnas.add(gene_name.lower())

            if len(example_records) < 5:
                example_records.append({
                    "CHROM": chrom, "START": start, "END": end, "STRAND": strand,
                    "GENE_NAME": gene_name,
                    "GENE_ID": gid_match.group(1) if gid_match else "NA",
                })

    print("[INFO] Example parsed GENCODE records (pre-deduplication):")
    for rec in example_records:
        print(" ", rec)
    print("[INFO] Multiple records per gene are expected — GENCODE annotates transcripts individually.")
    print(f"[DONE] Unique GENCODE miRNA gene symbols: {len(mirnas)}")

    with open(output_path, "w") as fh:
        for m in sorted(mirnas):
            fh.write(m + "\n")

    return output_path, mirnas

# =========================
# STEP 3 — GENERATE BED + DEDUPLICATE
# =========================
def generate_bed_from_gencode():
    print("\n" + "=" * 60)
    print("[STEP 3] Generating deduplicated BED file from GENCODE")
    print("=" * 60)
    print("[STORY] BED format encodes genomic intervals and enables coordinate-based operations.")
    print("[STORY] GENCODE contains multiple GTF entries per gene (one per transcript/feature).")
    print("[PROBLEM] Without deduplication, BEDTOOLS would extract the same locus multiple times,")
    print("          inflating sequence counts and biasing downstream analyses.")
    print("[SOLUTION] We retain only unique (CHROM, START, END, GENE_NAME) combinations.")
    time.sleep(2)

    raw_bed_path    = os.path.join(PROC_DIR, "mirna_coordinates_chr.bed")
    unique_bed_path = os.path.join(PROC_DIR, "mirna_coordinates_unique.bed")

    # --- Build raw BED ---
    bed_entries = []
    with gzip.open(GENCODE_GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if 'gene_type "miRNA"' not in line:
                continue
            cols  = line.strip().split("\t")
            chrom  = cols[0]
            start  = cols[3]
            end    = cols[4]
            strand = cols[6]
            match  = re.search(r'gene_name "([^"]+)"', line)
            if match:
                bed_entries.append((chrom, start, end, match.group(1).upper(), "0", strand))

    with open(raw_bed_path, "w") as fh:
        for b in bed_entries:
            fh.write("\t".join(b) + "\n")

    # --- Deduplicate ---
    seen = set()
    with open(raw_bed_path) as fin, open(unique_bed_path, "w") as fout:
        for line in fin:
            key = tuple(line.strip().split("\t")[:4])  # chr, start, end, name
            if key not in seen:
                seen.add(key)
                fout.write(line)

    print(f"[DONE] Raw BED entries    : {len(bed_entries)}")
    print(f"[DONE] Unique BED entries : {len(seen)}")
    print("[EXPECTATION] Unique loci should approximate annotated miRNA genes in GENCODE (~1800–1900).")
    return unique_bed_path

# =========================
# STEP 4 — PARSE miRBase
# =========================
def parse_mirbase():
    print("\n" + "=" * 60)
    print("[STEP 4] Parsing miRBase mature miRNA sequences")
    print("=" * 60)
    print("[STORY] miRBase catalogs mature miRNA sequences for all species.")
    print("[STORY] We extract only human entries (prefix 'hsa-') from the FASTA headers.")
    print("[STORY] Each entry represents a processed arm (5p or 3p) of a miRNA precursor.")
    time.sleep(2)

    output_path = os.path.join(PROC_DIR, "mirbase_hsa.txt")
    mirnas = set()

    with open(MIRBASE_FA) as fh:
        for line in fh:
            if line.startswith(">") and "hsa" in line:
                name = line.split()[0].lstrip(">").lower()
                mirnas.add(name)

    with open(output_path, "w") as fh:
        for m in sorted(mirnas):
            fh.write(m + "\n")

    print(f"[DONE] Human miRBase mature miRNAs: {len(mirnas)}")
    return output_path, mirnas

# =========================
# STEP 5 — SEQUENCE EXTRACTION (bedtools)
# =========================
def extract_sequences(bed_file):
    print("\n" + "=" * 60)
    print("[STEP 5] Extracting genomic miRNA sequences with BEDTOOLS")
    print("=" * 60)
    print("[STORY] We now move from annotation to sequences.")
    print("[STORY] Using the deduplicated GENCODE loci, BEDTOOLS retrieves the actual DNA")
    print("        sequence of each miRNA gene from the human reference genome.")
    print("[INFO]  Options used: -nameOnly (FASTA headers = gene name), -s (strand-aware)")
    time.sleep(2)

    import shutil
    if not shutil.which("bedtools"):
        print("[ERROR] bedtools is not installed or not in PATH.")
        print("        Linux/Mac: conda install -c bioconda bedtools")
        print("        Windows  : use WSL2 or a conda environment with bioconda")
        sys.exit(1)
    print("[INFO] bedtools detected ✔")

    genome_fa = os.path.join(RAW_DIR, "genome.fa")
    genome_gz = os.path.join(RAW_DIR, "genome.fa.gz")
    output_fa = os.path.join(PROC_DIR, "mirna_sequences.fa")

    if os.path.exists(genome_fa):
        print("[SKIP] genome.fa already exists ✔")
    elif os.path.exists(genome_gz):
        print("[INFO] Decompressing genome.fa.gz ...")
        subprocess.run(["gunzip", genome_gz], check=True)
        print("[DONE] genome.fa ready ✔")
    else:
        print("[ERROR] Reference genome not found. Expected:")
        print(f"        {genome_fa}  OR  {genome_gz}")
        sys.exit(1)

    cmd = ["bedtools", "getfasta",
           "-fi", genome_fa, "-bed", bed_file,
           "-fo", output_fa, "-nameOnly", "-s"]
    try:
        subprocess.run(cmd, check=True)
        print("[DONE] Sequence extraction completed ✔")
    except subprocess.CalledProcessError:
        print("[ERROR] BEDTOOLS execution failed. Check genome file and BED format.")
        sys.exit(1)

    # --- Quick FASTA validation ---
    seq_count, bad_headers = 0, 0
    with open(output_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                seq_count += 1
                h = line.strip().lstrip(">")
                if not h.startswith("MIR") and not h.startswith("LET"):
                    bad_headers += 1

    print(f"[VALIDATION] Sequences extracted : {seq_count}")
    print(f"[VALIDATION] Unexpected headers  : {bad_headers}")
    if bad_headers == 0:
        print("[SUCCESS] All FASTA headers are correctly normalized (MIR/LET format) ✔")
    else:
        print("[WARNING] Some headers are unexpected — review BED name field.")

    return output_fa

# =========================
# STEP 6 — CLEAN + NORMALIZE GENCODE
# =========================
def clean_and_normalize_gencode(raw_file, raw_set):
    print("\n" + "=" * 60)
    print("[STEP 6] Cleaning and normalizing GENCODE identifiers")
    print("=" * 60)
    print("[STORY] The raw GENCODE gene symbol list may contain non-canonical entries")
    print("        (e.g., ENSG IDs or read-through genes).")
    print("[SOLUTION] We retain only entries starting with 'mir' or 'let', then apply")
    print("           the shared normalization function to convert to MIR/LET format.")
    time.sleep(1)

    # --- Clean ---
    clean_path = os.path.join(PROC_DIR, "gencode_clean.txt")
    with open(raw_file) as fin, open(clean_path, "w") as fout:
        for line in fin:
            if line.startswith("mir") or line.startswith("let"):
                fout.write(line)

    clean_count = sum(1 for _ in open(clean_path))
    removed = raw_set - set(l.strip() for l in open(clean_path))
    print(f"[DEBUG] Raw symbols      : {len(raw_set)}")
    print(f"[DEBUG] After filtering  : {clean_count}")
    print(f"[DEBUG] Removed entries  : {len(removed)}")
    print(f"[DEBUG] Removed examples : {list(removed)[:8]}")

    # --- Normalize ---
    norm_path = os.path.join(PROC_DIR, "gencode_norm.txt")
    norm_set  = set()
    with open(clean_path) as fh:
        for line in fh:
            norm_set.add(normalize(line.strip()))

    with open(norm_path, "w") as fh:
        for m in sorted(norm_set):
            fh.write(m + "\n")

    print(f"[DONE] Normalized GENCODE miRNA symbols: {len(norm_set)}")
    print(f"[EXAMPLE] {list(norm_set)[:8]}")
    return norm_path, norm_set

# =========================
# STEP 7 — NORMALIZE miRBase
# =========================
def normalize_mirbase(raw_file):
    print("\n" + "=" * 60)
    print("[STEP 7] Normalizing miRBase identifiers")
    print("=" * 60)
    print("[STORY] miRBase uses a different naming convention from GENCODE.")
    print("[EXAMPLE]")
    print("  miRBase : hsa-miR-21-5p  →  normalized: MIR21")
    print("  miRBase : hsa-let-7g-3p  →  normalized: LET7G")
    print("[STORY] Removing species prefix and arm annotation allows gene-level comparison")
    print("        across databases that would otherwise appear to have low overlap.")
    time.sleep(2)

    norm_path = os.path.join(PROC_DIR, "mirbase_norm.txt")
    norm_set  = set()
    with open(raw_file) as fh:
        for line in fh:
            norm_set.add(normalize(line.strip()))

    with open(norm_path, "w") as fh:
        for m in sorted(norm_set):
            fh.write(m + "\n")

    print(f"[DONE] Normalized miRBase miRNA symbols: {len(norm_set)}")
    print(f"[EXAMPLE] {list(norm_set)[:8]}")
    return norm_path, norm_set

# =========================
# STEP 8 — PARSE miRGeneDB  ← FIX: file now reliably written before compare()
# =========================
def parse_mirgenedb():
    print("\n" + "=" * 60)
    print("[STEP 8] Parsing miRGeneDB — evolutionary validation layer")
    print("=" * 60)
    print("[STORY] Not every annotated miRNA is a bona fide miRNA gene.")
    print("[STORY] miRGeneDB applies strict evolutionary criteria to retain only")
    print("        miRNAs with conserved biogenesis hallmarks across species.")
    print("[STORY] We use this curated set as a high-confidence reference filter.")
    print("[INFO]  Parsing precursor-level entries (_pri / _pre) to avoid arm duplicates.")
    time.sleep(2)

    mirnas = set()
    debug_examples = []   # capture a few transformations for diagnostic output

    with open(MIRGENEDB_FILE) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 4:
                continue

            raw_name = cols[3]

            # Keep only precursor-level entries (_pre preferred; _pri also accepted).
            # This avoids counting 5p/3p arms as separate genes.
            if not (raw_name.endswith("_pre") or raw_name.endswith("_pri")):
                continue

            # Strip structural suffixes BEFORE calling normalize()
            # so the shared function receives a clean stem like "Mir-21" or "Let-7"
            clean = raw_name
            clean = re.sub(r'_(pre|pri)$', '', clean)   # remove _pre / _pri
            clean = re.sub(r'_[35]p\*?$',  '', clean)  # remove _5p / _3p if present

            # normalize()       : strips species prefix, unifies Mir-/miR-/mir-, uppercases
            # collapse_family() : strips -P<n><letters> paralog suffix  (MIR8-P1A → MIR8)
            norm_name = collapse_family(normalize(clean))
            mirnas.add(norm_name)

            if len(debug_examples) < 8:
                debug_examples.append((raw_name, clean, norm_name))

    print("[DEBUG] Normalization trace (raw_name → cleaned → normalized):")
    for raw, cln, nrm in debug_examples:
        print(f"  {raw:<40} → {cln:<28} → {nrm}")

    # ── Write output ──────────────────────────────────────────
    output_path = os.path.join(PROC_DIR, "mirgenedb_norm.txt")
    with open(output_path, "w") as fh:
        for m in sorted(mirnas):
            fh.write(m + "\n")

    print(f"[DONE] miRGeneDB high-confidence miRNAs : {len(mirnas)}")
    print(f"[OUTPUT] Saved to : {output_path}")
    return output_path, mirnas

# =========================
# STEP 9 — CROSS-DATABASE COMPARISON
# =========================
def compare_datasets(file_a, file_b, label_a, label_b, out_filename):
    """
    Generic pairwise comparison between two normalized miRNA sets.
    Writes results to PROC_DIR/<out_filename>.
    Returns (common, only_a, only_b) as sets.
    """
    set_a = set(open(file_a).read().split())
    set_b = set(open(file_b).read().split())

    common  = set_a & set_b
    only_a  = set_a - set_b
    only_b  = set_b - set_a

    print(f"\n  {label_a:<12}: {len(set_a):>5}")
    print(f"  {label_b:<12}: {len(set_b):>5}")
    print(f"  Overlap      : {len(common):>5}")
    print(f"  Only {label_a:<8}: {len(only_a):>5}")
    print(f"  Only {label_b:<8}: {len(only_b):>5}")

    output_path = os.path.join(PROC_DIR, out_filename)
    with open(output_path, "w") as fh:
        fh.write(f"# Comparison: {label_a} vs {label_b}\n\n")
        fh.write(f"COMMON ({len(common)})\n")
        fh.write("\n".join(sorted(common)) + "\n\n")
        fh.write(f"ONLY_{label_a.upper()} ({len(only_a)})\n")
        fh.write("\n".join(sorted(only_a)) + "\n\n")
        fh.write(f"ONLY_{label_b.upper()} ({len(only_b)})\n")
        fh.write("\n".join(sorted(only_b)) + "\n")

    print(f"  [OUTPUT] Saved to: {output_path}")
    return common, only_a, only_b

def run_comparisons(gencode_norm_file, mirbase_norm_file, mirgenedb_norm_file,
                    gencode_norm_set,  mirbase_norm_set,  mirgenedb_set):
    print("\n" + "=" * 60)
    print("[STEP 9] Cross-database comparisons")
    print("=" * 60)
    print("[STORY] Now that all three databases share the same naming convention,")
    print("        we can perform meaningful set comparisons.")
    print("[GOAL]  Identify which miRNAs are shared, database-specific, or universally supported.")
    time.sleep(2)

    print("\n--- 9a. GENCODE vs miRBase ---")
    print("[STORY] Both annotate human miRNAs from different angles (genomic vs sequence).")
    print("[EXPECTATION] Substantial but incomplete overlap — each database has unique entries.")
    compare_datasets(gencode_norm_file, mirbase_norm_file,
                     "GENCODE", "miRBase", "overlap_gencode_vs_mirbase.txt")

    print("\n--- 9b. GENCODE vs miRGeneDB ---")
    print("[STORY] miRGeneDB is a strict evolutionary filter.")
    print("[EXPECTATION] Most GENCODE miRNAs should overlap, but some low-confidence GENCODE")
    print("              entries will not be found in miRGeneDB.")
    compare_datasets(gencode_norm_file, mirgenedb_norm_file,
                     "GENCODE", "miRGeneDB", "overlap_gencode_vs_mirgenedb.txt")

    print("\n--- 9c. miRBase vs miRGeneDB ---")
    print("[STORY] miRBase contains entries with variable confidence levels.")
    print("[EXPECTATION] miRGeneDB validates a subset of miRBase entries.")
    compare_datasets(mirbase_norm_file, mirgenedb_norm_file,
                     "miRBase", "miRGeneDB", "overlap_mirbase_vs_mirgenedb.txt")

    # --- Three-way high-confidence set ---
    print("\n--- 9d. High-confidence set (GENCODE ∩ miRBase ∩ miRGeneDB) ---")
    print("[STORY] The intersection of all three databases represents the most robust miRNA set:")
    print("        annotated genomically (GENCODE), sequenced (miRBase), AND evolutionarily")
    print("        validated (miRGeneDB).")
    high_confidence = gencode_norm_set & mirbase_norm_set & mirgenedb_set
    hc_path = os.path.join(PROC_DIR, "high_confidence_mirnas.txt")
    with open(hc_path, "w") as fh:
        for m in sorted(high_confidence):
            fh.write(m + "\n")
    print(f"  High-confidence miRNAs: {len(high_confidence)}")
    print(f"  [OUTPUT] Saved to: {hc_path}")

    return high_confidence

# =========================
# STEP 10 — VENN DIAGRAM
# =========================
def generate_venn_plot(gencode_set, mirbase_set):
    print("\n" + "=" * 60)
    print("[STEP 10] Generating Venn diagram (GENCODE vs miRBase)")
    print("=" * 60)
    print("[STORY] Visual summary of the two-database overlap after normalization.")
    time.sleep(1)

    output_fig = os.path.join(PROC_DIR, "venn_mirna.png")
    plt.figure(figsize=(6, 6))
    venn2([gencode_set, mirbase_set], set_labels=("GENCODE", "miRBase"))
    plt.title("miRNA overlap: GENCODE vs miRBase (normalized)")
    plt.savefig(output_fig, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[DONE] Venn diagram saved to: {output_fig}")

# =========================
# STEP 11 — QC REPORT
# =========================
def generate_qc_report(bed_file, fasta_file, gencode_set, mirbase_set, mirgenedb_set):
    print("\n" + "=" * 60)
    print("[STEP 11] Generating QC report")
    print("=" * 60)
    time.sleep(1)

    bed_count   = sum(1 for _ in open(bed_file)) if os.path.exists(bed_file) else 0
    fasta_count = 0
    bad_headers = 0

    if os.path.exists(fasta_file):
        with open(fasta_file) as fh:
            for line in fh:
                if line.startswith(">"):
                    fasta_count += 1
                    h = line.strip().lstrip(">")
                    if not h.startswith("MIR") and not h.startswith("LET"):
                        bad_headers += 1

    overlap_gc_mb  = len(gencode_set & mirbase_set)
    overlap_gc_mgd = len(gencode_set & mirgenedb_set)
    overlap_mb_mgd = len(mirbase_set & mirgenedb_set)
    high_conf      = len(gencode_set & mirbase_set & mirgenedb_set)

    report_path = os.path.join(PROC_DIR, "qc_report.txt")
    with open(report_path, "w") as fh:
        fh.write("QC REPORT — Integrative miRNA Annotation Pipeline\n")
        fh.write("=" * 50 + "\n\n")
        fh.write(f"BED unique loci              : {bed_count}\n")
        fh.write(f"FASTA sequences extracted    : {fasta_count}\n")
        fh.write(f"FASTA unexpected headers     : {bad_headers}\n\n")
        fh.write(f"GENCODE miRNA symbols        : {len(gencode_set)}\n")
        fh.write(f"miRBase miRNA symbols        : {len(mirbase_set)}\n")
        fh.write(f"miRGeneDB curated symbols    : {len(mirgenedb_set)}\n\n")
        fh.write(f"Overlap GENCODE ∩ miRBase    : {overlap_gc_mb}\n")
        fh.write(f"Overlap GENCODE ∩ miRGeneDB  : {overlap_gc_mgd}\n")
        fh.write(f"Overlap miRBase ∩ miRGeneDB  : {overlap_mb_mgd}\n")
        fh.write(f"High-confidence (3-way)      : {high_conf}\n\n")
        fh.write("INTERPRETATION\n")
        fh.write("-" * 30 + "\n")
        fh.write(f"FASTA normalization  : {'OK' if bad_headers == 0 else 'ISSUES DETECTED'}\n")
        fh.write(f"Database agreement   : {'HIGH' if overlap_gc_mb > 1000 else 'MODERATE/LOW'}\n")

    print(f"[DONE] QC report saved to: {report_path}")

# =========================
# MAIN
# =========================
def main():
    os.chdir(BASE_DIR)

    print("\n" + "=" * 60)
    print("[PIPELINE START] Integrative miRNA Annotation Pipeline")
    print("=" * 60)
    print("[OBJECTIVE] Define a high-confidence set of human miRNA genes by integrating")
    print("            genomic annotation, sequence databases, and evolutionary curation.")
    print()
    print("[DATABASES]")
    print("  1. GENCODE   — genomic loci (GTF)")
    print("  2. miRBase   — mature miRNA sequences (FASTA)")
    print("  3. miRGeneDB — evolutionarily validated miRNA genes (BED)")
    print()
    print("[WORKFLOW]")
    print("  Acquire → Parse → Normalize → Compare → Validate → Curate → Report")
    time.sleep(3)

    ensure_dirs()
    debug_paths()

    # ── Steps ──────────────────────────────────────────────────
    acquire_data()                                          # Step 1

    gencode_raw_file, raw_set  = parse_gencode()           # Step 2
    bed_file                   = generate_bed_from_gencode()  # Step 3
    mirbase_raw_file, _        = parse_mirbase()           # Step 4
    output_fa                  = extract_sequences(bed_file)  # Step 5

    gencode_norm_file, gencode_norm_set = clean_and_normalize_gencode(gencode_raw_file, raw_set)  # Step 6
    mirbase_norm_file, mirbase_norm_set = normalize_mirbase(mirbase_raw_file)                      # Step 7
    mirgenedb_norm_file, mirgenedb_set  = parse_mirgenedb()                                        # Step 8

    high_confidence = run_comparisons(                     # Step 9
        gencode_norm_file, mirbase_norm_file, mirgenedb_norm_file,
        gencode_norm_set,  mirbase_norm_set,  mirgenedb_set,
    )

    generate_venn_plot(gencode_norm_set, mirbase_norm_set)  # Step 10
    generate_qc_report(                                     # Step 11
        bed_file, output_fa,
        gencode_norm_set, mirbase_norm_set, mirgenedb_set,
    )

    # ── Final summary table ────────────────────────────────────
    print("\n" + "=" * 60)
    print("[FINAL OUTPUT SUMMARY]")
    print("=" * 60)
    rows = [
        ("Step 2 – GENCODE parse",    "gencode_mirna.txt",                  "Raw miRNA gene symbols"),
        ("Step 3 – BED (unique)",     "mirna_coordinates_unique.bed",       "Deduplicated genomic loci"),
        ("Step 4 – miRBase parse",    "mirbase_hsa.txt",                    "Human mature miRNA names"),
        ("Step 5 – Sequences",        "mirna_sequences.fa",                 "Genomic miRNA sequences (FASTA)"),
        ("Step 6 – GENCODE norm",     "gencode_norm.txt",                   "Normalized GENCODE symbols"),
        ("Step 7 – miRBase norm",     "mirbase_norm.txt",                   "Normalized miRBase symbols"),
        ("Step 8 – miRGeneDB norm",   "mirgenedb_norm.txt",                 "Curated & normalized symbols"),
        ("Step 9a – Overlap GC/MB",   "overlap_gencode_vs_mirbase.txt",     "GENCODE ∩ miRBase"),
        ("Step 9b – Overlap GC/MGD",  "overlap_gencode_vs_mirgenedb.txt",   "GENCODE ∩ miRGeneDB"),
        ("Step 9c – Overlap MB/MGD",  "overlap_mirbase_vs_mirgenedb.txt",   "miRBase ∩ miRGeneDB"),
        ("Step 9d – High-confidence", "high_confidence_mirnas.txt",         "3-way intersection"),
        ("Step 10 – Venn diagram",    "venn_mirna.png",                     "GENCODE vs miRBase overlap"),
        ("Step 11 – QC report",       "qc_report.txt",                      "Pipeline quality metrics"),
    ]
    print(f"  {'Step':<28} {'File':<42} {'Description'}")
    print("  " + "-" * 100)
    for step, fname, desc in rows:
        print(f"  {step:<28} {fname:<42} {desc}")

    print("\n" + "=" * 60)
    print("[PIPELINE COMPLETE]")
    print("[RESULTS]")
    print(f"  High-confidence miRNAs (3-way): {len(high_confidence)}")
    print()
    print("[BIOLOGICAL CONCLUSIONS]")
    print("  • GENCODE provides genomic loci but no mature sequences")
    print("  • miRBase provides sequences but includes variable-confidence entries")
    print("  • miRGeneDB filters both to a biologically validated core set")
    print("  • The 3-way intersection is the most robust set for downstream analysis")
    print()
    print("[NEXT STEPS]")
    print("  • Integrate target databases: TargetScan, miRTarBase, TarBase")
    print("  • Add expression data from RNA-seq / small RNA-seq")
    print("  • Extend analysis to additional species")
    print("=" * 60)


if __name__ == "__main__":
    main()