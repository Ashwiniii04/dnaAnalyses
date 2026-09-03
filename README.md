# 🧬 GeneTrace

**Interactive DNA Fingerprinting & Sequence Identification System**

GeneTrace is a single-file Streamlit app that takes a raw DNA sequence through a full
identification pipeline: cleaning & validation, composition statistics, a deterministic
visual fingerprint (Chaos Game Representation), a shareable QR code, and a real remote
BLAST search against NCBI's `nt` nucleotide database.

> **Disclaimer:** GeneTrace is an educational and research prototype. Results should not
> be interpreted as clinical diagnoses, legal forensic evidence, or definitive species
> identification.

## Features

- **DNA Analyzer** — paste or upload a FASTA/TXT sequence; get cleaned sequence,
  base composition, GC%/AT%, and downloadable stats.
- **Chaos Game Representation (CGR)** — a deterministic visual fingerprint of the
  sequence (not decorative — identical sequences always produce identical patterns).
- **QR Code** — encodes the sequence plus metadata (SHA-256 ID, length, GC%, timestamp)
  as scannable JSON.
- **NCBI Identification** — submits a real `blastn` search against NCBI's `nt` database
  via Biopython (`Bio.Blast.NCBIWWW` / `NCBIXML`). No results are ever simulated or
  hardcoded.
- **Compare Sequences** — real global pairwise alignment (Biopython `Align.PairwiseAligner`)
  between two sequences, with identity %, mismatches, gaps, and side-by-side CGR plots.
- **Mystery DNA** — a set of built-in practice sequences for the classroom/self-study.

## Installation

```bash
git clone https://github.com/<your-username>/genetrace.git
cd genetrace
pip install -r requirements.txt
```

## Usage

```bash
streamlit run app.py
```

Then open the local URL Streamlit prints (usually `http://localhost:8501`).

Before running an NCBI BLAST search, enter your email in the sidebar under
**NCBI Settings** — NCBI requires an email address on every BLAST/Entrez request. No
account or API key is needed for BLAST searches.

⏳ Real BLAST searches against NCBI's `nt` database typically take 30 seconds to a few
minutes depending on NCBI server load — this is expected NCBI behavior, not a bug.

## Requirements

- Python 3.9+
- See [`requirements.txt`](requirements.txt) for packages
- Internet access to `blast.ncbi.nlm.nih.gov` for the NCBI Identification page

## Project structure

```
genetrace/
├── app.py            # entire application
├── requirements.txt
├── .gitignore
└── README.md
```

