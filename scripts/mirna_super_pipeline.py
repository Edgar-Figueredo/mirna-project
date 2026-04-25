import os
import re
import gzip
import urllib.request
import time
import sys
import subprocess

def ensure_package(pip_name, import_name=None):
    import_name = import_name or pip_name
    try:
        return __import__(import_name)
    except ImportError:
        print(f"[AUTO-INSTALL] Installing missing package: {pip_name}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pip_name])
        return __import__(import_name)


# =========================
# AUTO IMPORTS (CRITICAL)
# =========================
matplotlib = ensure_package("matplotlib")
import matplotlib.pyplot as plt
venn2 = ensure_package("matplotlib-venn", "matplotlib_venn").venn2

# =========================
# CONFIG (ROBUST PATHS)
# =========================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_DIR = os.path.join(BASE_DIR, "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
PROC_DIR = os.path.join(DATA_DIR, "processed")

GENCODE_GTF = os.path.join(RAW_DIR, "gencode.v49.chr.gtf.gz")
MIRBASE_FA = os.path.join(RAW_DIR, "mature.fa.gz")
MIRGENEDB_FILE = os.path.join(RAW_DIR, "mirgenedb.bed")

# URLs 
GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.chr.gtf.gz"
MIRBASE_URL = "https://www.mirbase.org/download/mature.fa"
MIRGENEDB_URL = "https://mirgenedb.org/gff/hsa?all=1&sort=pos"

def debug_paths():
    print("\n[DEBUG PATHS]")
    print("BASE_DIR:", BASE_DIR)
    print("RAW_DIR:", RAW_DIR)
    print("GENCODE_GTF:", GENCODE_GTF)
    print("MIRBASE_FA:", MIRBASE_FA)
    print("MIRGENEDB_FILE:", MIRGENEDB_FILE)
    print("-" * 40)
# =========================
# ENVIRONMENT CHECK
# =========================
REQUIRED = [
    "matplotlib",
    "matplotlib_venn",
]

missing = []

for pkg in REQUIRED:
    try:
        __import__(pkg)
    except ImportError:
        missing.append(pkg)

if missing:
    print("[ERROR] Missing dependencies:")
    for m in missing:
        print(" -", m)
    print("\nRun:")
    print("pip install -r env/requirements.txt")
    exit()

# =========================
# UTILS
# =========================
def ensure_dirs():
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(PROC_DIR, exist_ok=True)
    print("[DEBUG] Working directory:", os.getcwd())

def is_valid_gzip(filepath):
    import os
    if not os.path.exists(filepath):
        return False
    with open(filepath, "rb") as f:
        return f.read(2) == b"\x1f\x8b"

def download_file(url, output):
    if os.path.exists(output):
        print(f"[SKIP] File already exists: {output}")
        return

    print(f"[INFO] Downloading {url}")

    try:
        urllib.request.urlretrieve(url, output)
    except Exception as e:
        print(f"[ERROR] Download failed: {e}")

    print(f"[DONE] File downloaded: {output}")

# =========================
# GENCODE PARSE
# =========================
def parse_gencode():
    print("\n[STEP 3] Parsing GENCODE annotation...")
    print("[STORY] We extract miRNA genes based on gene_type = 'miRNA'.")
    print("[STORY] Important: GENCODE represents miRNAs as genomic loci, not sequences.")
    time.sleep(2)
    print("[INFO] Parsing GENCODE...")
    print("[INFO] Extracting fields: CHROM, START, END, STRAND, GENE_NAME, GENE_ID, TRANSCRIPT_NAME, TRANSCRIPT_ID")
    print("[INFO] Step explanation: Parsing GTF file and extracting miRNA entries based on gene_type = 'miRNA'.")
    time.sleep(1)

    example_records = []

    output = f"{PROC_DIR}/gencode_mirna.txt"
    mirnas = set()

    with gzip.open(GENCODE_GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            if 'gene_type "miRNA"' not in line:
                continue

            chrom = line.split("\t")[0]
            start = line.split("\t")[3]
            end = line.split("\t")[4]
            strand = line.split("\t")[6]

            gene_name_match = re.search(r'gene_name "([^"]+)"', line)
            gene_id_match = re.search(r'gene_id "([^"]+)"', line)
            transcript_name_match = re.search(r'transcript_name "([^"]+)"', line)
            transcript_id_match = re.search(r'transcript_id "([^"]+)"', line)

            if gene_name_match and len(example_records) < 5:
                example_records.append({
                    "CHROM": chrom,
                    "START": start,
                    "END": end,
                    "STRAND": strand,
                    "GENE_NAME": gene_name_match.group(1),
                    "GENE_ID": gene_id_match.group(1) if gene_id_match else "NA",
                    "TRANSCRIPT_NAME": transcript_name_match.group(1) if transcript_name_match else "NA",
                    "TRANSCRIPT_ID": transcript_id_match.group(1) if transcript_id_match else "NA"
                })

            if gene_name_match:
                name = gene_name_match.group(1).lower()
                mirnas.add(name)

    print("[INFO] Example parsed GENCODE entries (raw records before deduplication; note multiple entries can correspond to the same GENE_NAME due to transcripts):")
    for rec in example_records:
        print(rec)

    print("[INFO] Note: Final miRNA list is collapsed at the GENE_NAME level (unique gene symbols).")

    with open(output, "w") as out:
        for m in sorted(mirnas):
            out.write(m + "\n")

    print(f"[DONE] GENCODE miRNA: {len(mirnas)}")
    return output, mirnas

# =========================
# GENERATE BED FROM GENCODE
# =========================
def generate_bed_from_gencode():
    print("\n[STEP] Generating BED file from GENCODE annotation...")

    input_file = f"{PROC_DIR}/gencode_mirna.txt"
    output_file = f"{PROC_DIR}/mirna_coordinates_chr.bed"

    print("[STORY] We convert miRNA gene annotations into BED format.")
    print("[STORY] BED format enables genome coordinate-based operations (e.g., BEDTOOLS extraction).")

    time.sleep(2)

    bed_entries = []

    with gzip.open(GENCODE_GTF, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            if 'gene_type "miRNA"' not in line:
                continue

            cols = line.strip().split("\t")
            chrom = cols[0]
            start = cols[3]
            end = cols[4]
            strand = cols[6]

            match = re.search(r'gene_name "([^"]+)"', line)
            if match:
                gene = match.group(1).upper()

                bed_entries.append((chrom, start, end, gene, "0", strand))

    with open(output_file, "w") as out:
        for b in bed_entries:
            out.write("\t".join(b) + "\n")

    print(f"[DONE] BED file created: {output_file}")
    print(f"[INFO] Total entries: {len(bed_entries)}")

    return output_file
# =========================
# MIRBASE PARSE
# =========================
def parse_mirbase():
    print("\n[STEP 4.1] Parsing miRBase sequences...")
    print("[STORY] Extracting human mature miRNAs (hsa-) from FASTA headers.")
    print("[STORY] miRBase provides mature miRNAs including 5p and 3p arms.")
    time.sleep(2)

    fa_path = MIRBASE_FA

    if not os.path.exists(fa_path):
        print(f"[ERROR] miRBase FASTA not found: {fa_path}")
        print("👉 Ensure download step completed successfully.")
        exit()

    print("[INFO] Using miRBase FASTA (deterministic mode, no gzip handling).")

    output = f"{PROC_DIR}/mirbase_hsa.txt"
    mirnas = set()

    with open(fa_path) as f:
        for line in f:
            if line.startswith(">") and "hsa" in line:
                name = line.split()[0].replace(">", "").lower()
                mirnas.add(name)

    with open(output, "w") as out:
        for m in sorted(mirnas):
            out.write(m + "\n")

    print(f"[DONE] miRBase miRNA: {len(mirnas)}")

    return output
# =========================
# MIRGENEDB PARSE
# =========================

def parse_mirgenedb():
    print("\n[STEP 4.2] Parsing miRGeneDB (bona fide miRNA genes)...")
    print("[STORY] We now validate miRNA genes using evolutionary curated dataset.")
    print("[STORY] miRGeneDB distinguishes true miRNA genes from low-confidence annotations.")
    time.sleep(2)

    mirnas = set()

    with open(MIRGENEDB_FILE) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            if len(cols) < 4:
                continue

            # BED format: chr, start, end, name
            name = cols[3]
            mirnas.add(normalize(name))

    print(f"[DONE] miRGeneDB curated miRNAs: {len(mirnas)}")

    return mirnas
# =========================
# CLEAN GENCODE (REMOVE ENSG)
# =========================
def clean_gencode(input_file, raw_set=None):
    print("[INFO] Filtering GENCODE for valid miRNA gene symbols (removing non-miRNA entries and technical artifacts)...")
    import time
    print("[INFO] Step explanation: Filtering to retain canonical miRNA gene symbols (mir/let families).")
    time.sleep(1)

    output = f"{PROC_DIR}/gencode_clean.txt"

    with open(input_file) as f, open(output, "w") as out:
        for line in f:
            # keep only canonical miRNA gene names (mir / let families)
            if line.startswith("mir") or line.startswith("let"):
                out.write(line)

    if raw_set is not None:
        print(f"[DEBUG] Raw GENCODE miRNA: {len(raw_set)}")
        print(f"[DEBUG] After filtering: {sum(1 for _ in open(output))}")
        filtered_set = set(line.strip() for line in open(output))
        removed = raw_set - filtered_set
        print(f"[DEBUG] Removed during filtering: {len(removed)} entries")
        print(f"[DEBUG] Example removed entries: {list(removed)[:10]}")

    print("[DONE] GENCODE cleaned")
    return output


# =========================
# NORMALIZATION
# =========================
def normalize(name):
    name = name.strip()

    # remove species prefix
    name = name.replace("hsa-", "")

    # unify miR / let notation
    name = name.replace("miR-", "mir")
    name = name.replace("mir-", "mir")
    name = name.replace("let-", "let")

    # remove arm annotation (-5p / -3p)
    name = re.sub(r'-[35]p$', '', name)

    # remove any trailing isoform numbers like -1, -2 ONLY if duplicated (optional, keep for gene-level)
    # name = re.sub(r'-\d+$', '', name)

    # convert to uppercase gene-style nomenclature
    name = name.upper()

    return name


def normalize_file(input_file, output_file):
    print(f"[INFO] Applying normalization to: {input_file}")

    # Identify dataset source for clearer messaging
    if "gencode" in input_file.lower():
        source_label = "GENCODE"
    elif "mirbase" in input_file.lower():
        source_label = "miRBase"
    else:
        source_label = "Unknown dataset"

    print(f"[INFO] Step explanation: Normalizing miRNA names from {source_label} (removing species prefix, harmonizing naming, removing arm info, collapsing to gene-level).")

    time.sleep(1)

    dataset = set()

    with open(input_file) as f:
        for line in f:
            norm = normalize(line.strip())
            dataset.add(norm)

    print("[INFO] Example normalization (first 10 entries):")
    preview = list(dataset)[:10]
    for item in preview:
        print(item)

    with open(output_file, "w") as out:
        for m in sorted(dataset):
            out.write(m + "\n")

    print(f"[DONE] Normalized {source_label} miRNA names: {len(dataset)}")
    return output_file, dataset


# =========================
# COMPARISON
# =========================
def compare(gencode_file, mirbase_file):
    print("\n[STEP 6] Cross-database comparison...")
    print("[STORY] We now compare gene-level miRNA identifiers between GENCODE and miRBase.")
    print("[GOAL] Identify shared annotations and database-specific entries.")
    time.sleep(2)
    print("[INFO] Comparing datasets...")
    print("[INFO] Step explanation: Comparing normalized miRNA sets to identify overlap and database-specific entries.")
    time.sleep(1)

    gencode = set(open(gencode_file).read().split())
    mirbase = set(open(mirbase_file).read().split())

    common = gencode & mirbase
    only_gencode = gencode - mirbase
    only_mirbase = mirbase - gencode

    print("[INFO] Example overlap entries:")
    print(list(common)[:10])

    print("[INFO] Example GENCODE-only entries:")
    print(list(only_gencode)[:10])

    print("[INFO] Example miRBase-only entries:")
    print(list(only_mirbase)[:10])

    print("[INFO] Note: GENCODE count may differ from raw parsed values due to filtering and normalization steps.")
    print("\n=== FINAL RESULTS ===")
    print(f"GENCODE: {len(gencode)}")
    print(f"miRBase: {len(mirbase)}")
    print(f"Overlap: {len(common)}")
    print(f"Only GENCODE: {len(only_gencode)}")
    print(f"Only miRBase: {len(only_mirbase)}")

    output = f"{PROC_DIR}/final_overlap.txt"

    with open(output, "w") as out:
        out.write("COMMON\n")
        out.write("\n".join(sorted(common)) + "\n\n")

        out.write("ONLY_GENCODE\n")
        out.write("\n".join(sorted(only_gencode)) + "\n\n")

        out.write("ONLY_MIRBASE\n")
        out.write("\n".join(sorted(only_mirbase)))

    print(f"[OUTPUT] Saved to {output}")


# =========================
# FIGURE 1 - VENN DIAGRAM
# =========================
def generate_venn_plot(gencode_set, mirbase_set):
    print("\n[FIGURE 1] Generating Venn diagram (GENCODE vs miRBase)...")

    output_fig = f"{PROC_DIR}/venn_mirna.png"

    plt.figure(figsize=(6,6))
    venn2([gencode_set, mirbase_set], set_labels=("GENCODE", "miRBase"))
    plt.title("miRNA overlap: GENCODE vs miRBase")

    plt.savefig(output_fig, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[DONE] Venn diagram saved to: {output_fig}")


# =========================
# QC REPORT
# =========================
def generate_qc_report(bed_file, fasta_file, gencode_set, mirbase_set):
    print("\n[QC REPORT] Generating automatic quality control summary...")
    time.sleep(1)

    report_file = f"{PROC_DIR}/qc_report.txt"

    bed_count = sum(1 for _ in open(bed_file)) if os.path.exists(bed_file) else 0

    fasta_count = 0
    bad_headers = 0

    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                fasta_count += 1
                h = line.strip().replace(">", "")
                if not h.startswith("MIR") and not h.startswith("LET"):
                    bad_headers += 1

    overlap = len(gencode_set & mirbase_set)

    with open(report_file, "w") as out:
        out.write("QC REPORT - miRNA PIPELINE\n")
        out.write("="*40 + "\n")
        out.write(f"BED entries: {bed_count}\n")
        out.write(f"FASTA sequences: {fasta_count}\n")
        out.write(f"FASTA bad headers: {bad_headers}\n")
        out.write(f"GENCODE-miRBase overlap: {overlap}\n")

        out.write("\nINTERPRETATION:\n")
        if bad_headers == 0:
            out.write("- FASTA normalization: OK\n")
        else:
            out.write("- FASTA normalization: ISSUES DETECTED\n")

        if overlap > 1000:
            out.write("- Database agreement: HIGH\n")
        else:
            out.write("- Database agreement: MODERATE/LOW\n")

    print(f"[DONE] QC report saved to: {report_file}")

# =========================
# MAIN PIPELINE
# =========================
def main():
    # Ensure working directory is project root
    os.chdir(BASE_DIR)
    print("[DEBUG] Forced working directory to BASE_DIR:", os.getcwd())
    print("\n" + "="*60)
    print("[PIPELINE START] miRNA cross-database analysis")
    print("="*60)
    print("[STORY] We begin by gathering reference annotations from two major resources:")
    print("[STORY] - GENCODE: comprehensive genome annotation (gene-centric)")
    print("[STORY] - miRBase: curated miRNA sequences (sequence-centric)")
    time.sleep(2)

    ensure_dirs()
    debug_paths()
    
    print("\n[STEP 1] Acquiring GENCODE annotation...")
    print("[INFO] GENCODE GTF contains genomic coordinates, gene/transcript structure, and biotypes (e.g., miRNA).")
    print("[INFO] It provides gene-level annotation, not mature miRNA sequences or targeting data.")
    time.sleep(2)

    # =========================
    # DOWNLOAD SECTION
    # =========================
    download_file(GENCODE_URL, GENCODE_GTF)
    download_file(MIRBASE_URL, MIRBASE_FA)
    download_file(MIRGENEDB_URL, MIRGENEDB_FILE)

    print("\n[STEP 2] Acquiring miRBase data...")
    print("[INFO] miRBase provides mature miRNA sequences (including -5p and -3p arms).")
    print("[INFO] It does NOT provide target genes or affinity data.")
    print("[INFO] It focuses on sequence-level annotation, not genomic structure.")
    time.sleep(2)

    mirbase_gz = f"{RAW_DIR}/mature.fa.gz"

    if os.path.exists(mirbase_gz):
        print("[SKIP] Using existing miRBase gzip file")
    elif os.path.exists(MIRBASE_FA):
        print("[SKIP] Using existing miRBase FASTA file")
    else:
        print("[INFO] miRBase file not found, attempting download (FASTA)...")
        download_file(MIRBASE_URL, MIRBASE_FA)

    # =========================
    # DATA INTEGRITY CHECKPOINT
    # =========================
    print("\n[CHECKPOINT] Validating downloaded files...")

    def check_file(path, min_size_mb=None):
        if not os.path.exists(path):
            print(f"[ERROR] Missing file: {path}")
            return False

        size_mb = os.path.getsize(path) / (1024 * 1024)

        if min_size_mb is not None and size_mb < min_size_mb:
            print(f"[ERROR] File too small ({size_mb:.2f} MB): {path}")
            return False

        print(f"[OK] {path} ({size_mb:.2f} MB)")
        return True


    check_file(GENCODE_GTF)
    check_file(MIRBASE_FA)
    check_file(MIRGENEDB_FILE)
    # ---- PARSE ----
    gencode_raw_file, raw_set = parse_gencode()
    bed_file = generate_bed_from_gencode()
    # Deduplicate BED entries
    print("[INFO] Collapsing BED entries to unique coordinates (removing duplicates)...")
    print("[STORY] GENCODE provides multiple entries per miRNA due to transcripts and feature annotations.")
    print("[PROBLEM] This creates duplicated genomic coordinates for the same miRNA gene.")
    print("[IMPACT] If not corrected, sequence extraction would produce redundant sequences and inflate downstream analyses.")
    print("[SOLUTION] We collapse entries to retain unique genomic loci (CHR, START, END, GENE_NAME).")
    print("[BIOLOGICAL CONTEXT] This step converts transcript-level annotation into gene-level representation.")
    time.sleep(2)

    unique_bed = f"{PROC_DIR}/mirna_coordinates_unique.bed"
    seen = set()
    with open(bed_file) as f, open(unique_bed, "w") as out:
        for line in f:
            key = tuple(line.strip().split("\t")[:4])  # chr, start, end, name
            if key not in seen:
                seen.add(key)
                out.write(line)
    print(f"[DONE] Unique BED entries: {len(seen)}")
    print("[INTERPRETATION] This number represents distinct miRNA genomic loci after removing redundant annotations.")
    print("[EXPECTATION] This value should approximate the number of annotated miRNA genes in GENCODE (~1800-1900).")
    # update bed_file to use the deduplicated version
    bed_file = unique_bed
    mirbase_raw = parse_mirbase()

    # SEQUENCE EXTRACTION STEP
    print("\n[STEP 4.5] Extracting genomic miRNA sequences...")
    print("[STORY] We now move from annotation to sequence-level data.")
    print("[STORY] Using genomic coordinates from GENCODE, we retrieve the underlying DNA sequences.")
    print("[INFO] This step uses BEDTOOLS to map miRNA loci back to the reference genome.")
    print("[INFO] Output: FASTA file containing genomic miRNA sequences.")
    time.sleep(2)

    # ---- CHECK BEDTOOLS ----
    import shutil
    if not shutil.which("bedtools"):
        print("[ERROR] bedtools is not installed or not in PATH.")
        print("👉 Install with: brew install bedtools")
        exit()

    print("[INFO] bedtools detected ✔")

    # ---- CHECK REQUIRED FILES ----
    # bed_file already defined and deduplicated above
    genome_fa = f"{RAW_DIR}/genome.fa"
    genome_gz = f"{RAW_DIR}/genome.fa.gz"

    print("[INFO] Checking required input files...")
    print(f"  BED: {bed_file}")
    print(f"  Genome FA: {genome_fa}")
    print(f"  Genome GZ: {genome_gz}")
    print("[INFO] BED file should now exist from previous generation step.")

    if not os.path.exists(bed_file):
        print(f"[ERROR] BED file not found: {bed_file}")
        print("👉 This should be generated in previous steps.")
        exit()

    # ---- HANDLE GENOME ----

    if os.path.exists(genome_fa):
        print("[SKIP] genome.fa already exists ✔")

    elif os.path.exists(genome_gz):
        print("[INFO] genome.fa not found, decompressing genome.fa.gz...")
        subprocess.run(["gunzip", genome_gz], check=True)
        print("[DONE] genome.fa ready ✔")

    else:
        print("[ERROR] Genome file not found.")
        print("👉 Expected:")
        print(f"   {genome_fa}")
        print("   OR")
        print(f"   {genome_gz}")
        exit()

    # ---- RUN EXTRACTION SCRIPT ----
    print("[STORY] Extracting strand-aware sequences (-s) and preserving miRNA names (-name).")
    time.sleep(1)

    # ---- CHECK POINT OF BEDTOOLS ----

    output_fa = f"{PROC_DIR}/mirna_sequences.fa"

    print(f"[INFO] Output FASTA: {output_fa}")
    print("[INFO] Running BEDTOOLS getfasta...")

    cmd = [
        "bedtools", "getfasta",
        "-fi", genome_fa,
        "-bed", bed_file,
        "-fo", output_fa,
        "-nameOnly",
        "-s"
    ]

    try:
        subprocess.run(cmd, check=True)
        print("[DONE] Sequence extraction completed successfully ✔")
    except subprocess.CalledProcessError:
        print("[ERROR] BEDTOOLS execution failed.")
        print("👉 Check genome file, BED format, and installation.")
        exit()

    # =========================
    # FASTA VALIDATION
    # =========================
    print("\n[FASTA VALIDATION]")
    print("[INFO] Validating extracted miRNA FASTA file...")

    seq_count = 0
    bad_headers = 0

    with open(output_fa) as f:
        for line in f:
            if line.startswith(">"):
                seq_count += 1
                header = line.strip().replace(">", "")
                if not header.startswith("MIR") and not header.startswith("LET"):
                    bad_headers += 1

    print(f"[RESULT] Total sequences: {seq_count}")
    print(f"[RESULT] Headers not normalized (unexpected): {bad_headers}")

    if bad_headers == 0:
        print("[SUCCESS] All FASTA headers are correctly normalized (MIR/LET format)")
    else:
        print("[WARNING] Some FASTA headers are not normalized. Check BED generation step.")

    # ---- POST EXTRACTION SUMMARY ----
    time.sleep(1)
    print("\n[STORY] We now have three complementary biological layers:")
    print("  1. Genomic annotation (GENCODE)")
    print("  2. Mature miRNA sequences (miRBase)")
    print("  3. Genomic DNA sequences (BEDTOOLS extraction)")
    time.sleep(2)
    print("\n[STEP 4.5] Extracting genomic miRNA sequences...")
    print("[STORY] We now move from annotation to sequence-level data.")
    print("[STORY] Using genomic coordinates from GENCODE, we retrieve the underlying DNA sequences.")
    print("[INFO] This step uses BEDTOOLS to map miRNA loci back to the reference genome.")
    print("[INFO] Output: FASTA file containing genomic miRNA sequences.")
    time.sleep(2)

    time.sleep(1)
    print("[STORY] We now have three complementary layers:")
    print("- Genomic annotation (GENCODE)")
    print("- Mature miRNA sequences (miRBase)")
    print("- Genomic DNA sequences (BEDTOOLS extraction)")
    time.sleep(2)

    print("[STORY] At this point:")
    print("- GENCODE → gene-level miRNA annotation")
    print("- miRBase → mature miRNA sequences (5p/3p variants)")
    print("[STORY] Neither database directly provides miRNA target interactions.")
    print("[STORY] For that, databases like TargetScan, miRTarBase, or TarBase are required.")
    time.sleep(3)

    # ---- CLEAN ----
    gencode_clean = clean_gencode(gencode_raw_file, raw_set)

    # ---- NORMALIZE ----
    print("\n[STEP 5] Normalization of miRNA identifiers...")
    print("[PROBLEM] GENCODE and miRBase use different naming conventions.")
    print("[EXAMPLE]")
    print("  GENCODE: MIR21")
    print("  miRBase: hsa-miR-21-5p")
    print("[ISSUE] These refer to the same biological entity but appear different computationally.")
    print("[SOLUTION] We normalize names to a unified gene-level nomenclature (uppercase MIR/LET format).")
    print("[TRANSFORMATION EXAMPLES | Gene-level harmonization]")
    print("  hsa-miR-21-5p  →  MIR21")
    print("  hsa-let-7g-3p →  LET7G")
    print("  MIR1302-2     →  MIR1302-2")
    time.sleep(3)

    gencode_norm_file, gencode_norm_set = normalize_file(
        gencode_clean,
        f"{PROC_DIR}/gencode_norm.txt"
    )

    mirbase_norm_file, mirbase_norm_set = normalize_file(
        mirbase_raw,
        f"{PROC_DIR}/mirbase_norm.txt"
    )

    # ---- COMPARE ----
    compare(gencode_norm_file, mirbase_norm_file)

    # =========================
    # FIGURE GENERATION
    # =========================
    generate_venn_plot(gencode_norm_set, mirbase_norm_set)

    # =========================
    # FINAL OUTPUT SUMMARY TABLE
    # =========================
    print("\n[FINAL SUMMARY TABLE]")

    header = f"{'Process':<25}{'Output File':<50}{'Description'}"
    print(header)
    print("-" * len(header))

    rows = [
        ("GENCODE parsing", gencode_raw_file, "Raw miRNA gene names from GENCODE"),
        ("BED generation", bed_file, "Genomic coordinates of miRNA loci (deduplicated)"),
        ("miRBase parsing", mirbase_raw, "Human mature miRNAs (hsa)"),
        ("Sequence extraction", output_fa, "Genomic miRNA sequences (FASTA)"),
        ("GENCODE cleaned", gencode_clean, "Filtered canonical miRNA genes (mir/let)"),
        ("GENCODE normalized", gencode_norm_file, "Gene-level normalized miRNA names"),
        ("miRBase normalized", mirbase_norm_file, "Normalized miRNA names"),
        ("Final overlap", f"{PROC_DIR}/final_overlap.txt", "Shared and unique miRNAs")
    ]

    for process, file, desc in rows:
        print(f"{process:<25}{file:<50}{desc}")

    # =========================
    # QC REPORT CALL
    # =========================
    generate_qc_report(
        bed_file,
        output_fa,
        gencode_norm_set,
        mirbase_norm_set
    )

    print("\n" + "="*60)
    print("[PIPELINE COMPLETE]")
    print("[SUMMARY]")
    print("- GENCODE: gene-level annotation (no sequences)")
    print("- miRBase: mature miRNA sequences (5p/3p)")
    print("- Neither provides miRNA-target interaction data")
    print("- Additional resources needed: TargetScan, miRTarBase, TarBase")
    print("="*60)


if __name__ == "__main__":
    main()
