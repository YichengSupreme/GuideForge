# CRISPR Target Automation

A modular pipeline for CRISPR gRNA choosing using UCSC Genome Browser and IDT CRISPR tools.

## 🔄 Pipeline Overview

```mermaid
flowchart TD
    A["Input: Genomic Coordinates<br/>chr17:7668402-7668521:+"] --> B["get_ucsc_sequences.py"]

    B --> C{Options}
    C -->|Basic| D["Fetch Upstream/Downstream<br/>Sequences from UCSC"]
    C -->|--scan-pam| E["Scan for SpCas9 PAM sites<br/>NGG patterns"]
    C -->|--qc| F["Apply Quality Control<br/>GC content, homopolymers, etc."]

    D --> G["Upstream_sequences.txt<br/>Downstream_sequences.txt"]
    E --> H["CRISPR_candidates.txt<br/>FASTA format"]
    F --> I["CRISPR_candidates_qc.csv<br/>QC results"]

    G --> J["IDT Analysis"]
    H --> J
    I --> J

    J --> K["idt_batch_crispr.py<br/>Batch API calls"]
    K --> L["CRISPR_candidates_idt.csv<br/>Ranked results with on/off-target scores"]
    L --> M["Final Results<br/>Best CRISPR sites"]

    style A fill:#e1f5fe,stroke:#0277bd,stroke-width:1px
    style B fill:#fff3e0,stroke:#ffb300,stroke-width:1px
    style C fill:#fce4ec
    style K fill:#fff3e0,stroke:#ffb300,stroke-width:1px
    style J fill:#f3e5f5,stroke:#8e24aa,stroke-width:1px
    style M fill:#c8e6c9,stroke:#2e7d32,stroke-width:1px
```

## 🧬 Modular Scripts that can be ran as a pipeline or as individual tasks

#### 1. **`get_ucsc_sequences.py`** - Fetch Sequences Only
**When to use:** You want to get DNA sequences from UCSC but analyze them elsewhere

```bash
# Basic usage 
python get_ucsc_sequences.py targets.txt

# Single coordinate
python get_ucsc_sequences.py chr17:7668402-7668421:+

# Custom distances
python get_ucsc_sequences.py targets.txt --up 200 --down 200

# Different genome (default = hg38)
python get_ucsc_sequences.py targets.txt --genome hg19 

# Scan for PAM sites
python get_ucsc_sequences.py chr17:7668402-7668421:+ --scan-pam
```

**Input:** Genomic coordinates only (e.g., `chr17:7668402-7668421:+`)
**Output:** 
- `Upstream_sequences.txt` - Upstream DNA sequences (FASTA format)
- `Downstream_sequences.txt` - Downstream DNA sequences (FASTA format)
- `CRISPR_candidates.txt` - PAM sites in FASTA format (if `--scan-pam` used)
- `CRISPR_candidates_qc.csv` - QC results for PAM sites (if `--qc` used)

**PAM Scanning:** Use `--scan-pam` to automatically find SpCas9 PAM sites (NGG) in the sequences and output CRISPR candidates in CSV format.

---

#### 2. **`idt_batch_crispr.py`** - Analyze Existing Sequences
**When to use:** You already have FASTA files and want IDT analysis

**Uses:** https://eu.idtdna.com/site/order/designtool/index/CRISPR_SEQUENCE

```bash
# Analyze one file
python idt_batch_crispr.py sequences.txt

# Analyze multiple files
python idt_batch_crispr.py upstream.txt downstream.txt exon.txt
```

**Input:** FASTA-style .txt files
**Output:** `{filename}_idt.csv` with ranked results
- Contains columns: `sequence_name`, `dna_sequence`, `on_target_score`, `off_target_score`, `on_plus_off`
- Sorted by `on_plus_off` (highest scores first)

---

#### 3. **`run_CRISPR_target_automation.py`** - Complete Pipeline
**When to use:** You want the full workflow from chromosome coordinates to ranked results

```bash
# Complete pipeline
python run_CRISPR_target_automation.py targets.txt

# Custom settings
python run_CRISPR_target_automation.py targets.txt --up 200 --down 200 --batch-size 5
```

**Input:** Genomic coordinates only
**Output:** Complete analysis with ranked results
- `Upstream_sequences_idt.csv` - Ranked upstream CRISPR sites
- `Downstream_sequences_idt.csv` - Ranked downstream CRISPR sites
- `CRISPR_candidates_idt.csv` - Ranked PAM sites (if --scan-pam used)

---

## 🎯 Which Script Should I Use?

| Your Situation | Use This Script | Why |
|----------------|-----------------|-----|
| I have coordinates, want full analysis | `run_CRISPR_target_automation.py targets.txt --scan-pam --qc` | Complete pipeline with PAM + QC |
| I have coordinates, want sequences + PAM + QC | `get_ucsc_sequences.py targets.txt --scan-pam --qc` | Find and filter CRISPR targets |
| I have coordinates, want sequences + PAM only | `get_ucsc_sequences.py targets.txt --scan-pam` | Find CRISPR targets (no filtering) |
| I have coordinates, want up/downstream sequences only | `get_ucsc_sequences.py targets.txt` | Just fetch sequences |
| I have FASTA files, want IDT analysis | `idt_batch_crispr.py file1.txt [file2.txt]` | Analyze existing sequences |
| I want to test IDT connection | `idt_batch_crispr.py` | Quick connectivity test |

---

## ⚙️ Configuration

Edit `config.sh` to customize:

```bash
# Sequence distances
UPSTREAM_DISTANCE=100
DOWNSTREAM_DISTANCE=100

# Genome assembly
GENOME_ASSEMBLY="hg38"

# IDT settings
IDT_BATCH_SIZE=10
IDT_TIMEOUT=60
IDT_SESSION_COOKIE="ASP.NET_SessionId=YOUR_SESSION_ID_HERE; ..."

# Retry settings
UCSC_RETRIES=3

# Output verbosity
VERBOSE=1  # Set to 0 for minimal output
```

**Important**: You must update the `IDT_SESSION_COOKIE` in `config.sh` with your valid IDT session cookie.

---

## 📊 Output Files

#### **Basic Sequence Files:**
- `Upstream_sequences.txt` - Upstream DNA sequences (FASTA format)
- `Downstream_sequences.txt` - Downstream DNA sequences (FASTA format)

#### **PAM Analysis Files:**
- `CRISPR_candidates.txt` - PAM sites in FASTA format (if --scan-pam used)
- `CRISPR_candidates_qc.csv` - QC results for PAM sites (if --qc used)

#### **IDT Analysis Files:**
- `CRISPR_candidates_idt.csv` - Ranked PAM sites (if --scan-pam used)

**CRISPR_candidates_idt columns:**
- `sequence_name` - Original identifier
- `dna_sequence` - The actual DNA sequence
- `on_target_score` - IDT on-target analysis
- `off_target_score` - IDT off-target analysis
- `on_plus_off` - Combined score (higher = better)

---

## 🚀 Quick Start

1. **Set up your cookie** in `config.sh`
2. **Create a targets.txt file** with gene coordinates:
   ```
   chr17:7668402-7668421:+
   ```
3. **Run the pipeline:**
   ```bash
   python run_CRISPR_target_automation.py targets.txt --scan-pam --qc 
   ```
4. **Open the results** in Excel and sort by `on_plus_off` column

---

## 📋 Input Format

The targets.txt file supports coordinates format only:

**Coordinates:**
```
chr14:103928378-103928397:+
chr1:1000000-1000020:-
chr17:7668402-7668421:+
```

**Format:** `chr:start-end:strand`

**Note:** Gene names are not supported to avoid ambiguity. Use coordinates instead.

---

## 🔑 How to Get Your IDT Session Cookie

The IDT CRISPR API requires a temporary browser session cookie to authenticate your requests. Follow these simple steps to get yours:

#### Instructions

1. **Open IDT CRISPR Designer**
   - Go to: https://eu.idtdna.com/site/order/designtool/index/CRISPR_SEQUENCE

2. **Open Developer Tools**
   - **Mac**: `Cmd + Option + I`
   - **Windows/Linux**: `Ctrl + Shift + I`

3. **Navigate to Network Tab**
   - Click on the "Network" tab in developer tools

4. **Reload the Page**
   - Press `F5` or `Cmd+R` to refresh
   - This populates the Network panel with requests

5. **Find IDT Request**
   - Look for any request to `eu.idtdna.com`
   - Click on it to inspect

   ![IDT Cookie Location](Cookie.png)
   *The Cookie header is highlighted in red - this is what you need to copy!*

6. **Copy Cookie Value**
   - In the right panel: **Headers** → **Request Headers**
   - Find the line starting with `Cookie:`
   - Copy the entire cookie string

7. **Update config.sh**
   ```bash
   IDT_SESSION_COOKIE="ASP.NET_SessionId=your_session_id; Anon=your_anon_id; ..."
   ```

#### ⚠️ Important Notes

- **Cookies expire every few hours** — if the script stops returning scores, just refresh your cookie
- **Keep your cookie private** — don't share it publicly
- **Example cookie format**:
  ```
  ASP.NET_SessionId=dzzwdv3adomtkqept1zgp0hc; Anon=t6ocoWScxF8=; ARRAffinity=874c5298ae0e2eca12812a980102a414521df46497427d5bbed67654bd42654b
  ```

#### 🔧 Troubleshooting

**❌ Getting "Authentication failed" errors?**
- Your cookie has expired → Get a fresh one using the steps above

**❌ Can't find the Cookie header?**
- Make sure you're looking in **Request Headers** (not Response Headers)!
- Scroll down once inside Headers (You should see General, Response Headers, Request Headers)
- Try reloading the page and looking for requests to `eu.idtdna.com`

## 🔧 Requirements

- Python 3.6+
- Required Python packages: `requests`, `pandas`
- Valid IDT session cookie (update in `config.sh`)
- see requirements.txt

---

## 🛠️ Troubleshooting

**Session expired**: Update the `IDT_SESSION_COOKIE` in `config.sh`

**Network timeouts**: Increase `IDT_TIMEOUT` in `config.sh`

**Verbose output**: Set `VERBOSE=1` in `config.sh`

---

## 🧪 Lab Usage

1. Copy the entire `CRISPR_target_automation` folder
2. Update the IDT session cookie in `config.sh`
3. Create your targets.txt file
4. Run: `run_CRISPR_target_automation.py targets.txt --scan-pam --qc`

## Future Notes
- In future versions, grepping IDT cookie could be automated by fetching an anonymous session cookie via requests.get(), but for stability and reproducibility, the current version uses a user-supplied cookie.