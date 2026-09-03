"""
================================================================================
GeneTrace — Interactive DNA Fingerprinting & Sequence Identification System
================================================================================

A single-file Streamlit application. Everything lives in this one file:
app.py — no other modules, no database of any kind.

REQUIRED PACKAGES (install with pip):

    pip install streamlit biopython numpy pandas matplotlib qrcode pillow

RUN WITH:

    streamlit run app.py

--------------------------------------------------------------------------------
WORKFLOW
--------------------------------------------------------------------------------
    DNA Sequence
        -> Cleaning & Validation
        -> DNA Statistics
        -> CGR Fingerprint
        -> QR Code
        -> NCBI BLAST
        -> Sequence Matching
        -> Result Interpretation

--------------------------------------------------------------------------------
IMPORTANT NOTES ON NCBI BLAST
--------------------------------------------------------------------------------
This app performs REAL remote BLAST searches against NCBI's "nt" nucleotide
database using Bio.Blast.NCBIWWW.qblast(). It never fabricates, hardcodes, or
simulates results. A real BLAST search against "nt" can take anywhere from
about 30 seconds to several minutes depending on NCBI server load — this is
normal and expected, not a bug.

DISCLAIMER: GeneTrace is an educational and research prototype. Results
should not be interpreted as clinical diagnoses, legal forensic evidence, or
definitive species identification.
================================================================================
"""

import io
import json
import hashlib
import socket
import datetime
import warnings

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import qrcode
from PIL import Image

from Bio import Align
from Bio.Blast import NCBIWWW, NCBIXML

warnings.filterwarnings("ignore", category=UserWarning)

# ==============================================================================
# PAGE CONFIG & DARK THEME
# ==============================================================================
st.set_page_config(
    page_title="GeneTrace",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

CUSTOM_CSS = """
<style>
    .stApp {
        background-color: #0b0f14;
        color: #e6f1ef;
    }
    .block-container {padding-top: 1.5rem; padding-bottom: 2rem;}
    section[data-testid="stSidebar"] {
        background: linear-gradient(180deg, #041b1a 0%, #062b28 100%);
        border-right: 1px solid #113330;
    }
    section[data-testid="stSidebar"] * {
        color: #e6f1ef !important;
    }
    div[data-testid="stMetric"] {
        background: #0f1e1c;
        border: 1px solid #16463f;
        border-radius: 10px;
        padding: 10px 16px;
    }
    div[data-testid="stMetricLabel"] {
        color: #7fd8c9 !important;
    }
    h1, h2, h3, h4 {
        color: #2dd4bf;
    }
    .stButton > button[kind="primary"] {
        background-color: #0f766e;
        border-color: #0f766e;
        color: #f0fdfa;
    }
    .stButton > button[kind="primary"]:hover {
        background-color: #14b8a6;
        border-color: #14b8a6;
    }
    .gt-card {
        background: #0f1e1c;
        border: 1px solid #16463f;
        border-radius: 12px;
        padding: 18px 20px;
        margin-bottom: 14px;
    }
    .gt-best-match {
        background: #08211d;
        border: 2px solid #2dd4bf;
        border-radius: 12px;
        padding: 18px 20px;
    }
    code, .stCodeBlock {
        background-color: #08110f !important;
    }
</style>
"""
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)

# ==============================================================================
# CONSTANTS
# ==============================================================================
VALID_BASES = set("ATGC")

# Chaos Game Representation corners, exactly as specified.
CGR_CORNERS = {
    "A": (0.0, 0.0),
    "T": (1.0, 0.0),
    "G": (0.0, 1.0),
    "C": (1.0, 1.0),
}

MYSTERY_SEQUENCES = {
    "Mystery Sample #1 (short, GC-rich)": (
        "GCGCGCATGCGGCTAGCGCGATCGGCTAGCTAGCGCGCTAGCTAGCGCGATCGATCGGCGCTAGCTAGC"
    ),
    "Mystery Sample #2 (AT-rich)": (
        "ATATATTTAAATATTTAAATATTAAATTTATAAATATTTAAATTTAAATATATTTAAATTTAAATATA"
    ),
    "Mystery Sample #3 (repetitive motif)": (
        "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
    ),
    "Mystery Sample #4 (mixed composition)": (
        "TTGGCACCTCCCTGGCAGAACATTTCTTACACAGTTCTCCACGTAGATCCTGCTCTGGCCTCCCAAAG"
    ),
    "Mystery Sample #5 (palindrome-like)": (
        "GAATTCGGATCCAAGCTTGAATTCGGATCCAAGCTTGAATTCGGATCCAAGCTTGAATTCGGATCC"
    ),
}


# ==============================================================================
# SESSION STATE INITIALIZATION
# ==============================================================================
def init_session_state():
    defaults = {
        "current_sequence": "",
        "current_header": "",
        "current_stats": None,
        "cgr_fig_bytes": None,
        "qr_image_bytes": None,
        "qr_metadata_json": None,
        "blast_results_df": None,
        "blast_best_match": None,
        "blast_raw_status": "No BLAST search performed yet.",
        "analysis_count": 0,
        "comparison_results_text": None,
        "ncbi_email": "",
        "ncbi_api_key": "",
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


# ==============================================================================
# CORE FUNCTIONS: CLEANING, VALIDATION, STATISTICS
# ==============================================================================
def clean_sequence(raw_text: str) -> str:
    """Remove FASTA headers, whitespace, and digits; return uppercase A/T/G/C(/other) string."""
    if not raw_text:
        return ""
    lines = raw_text.strip().splitlines()
    seq_lines = [ln for ln in lines if not ln.strip().startswith(">")]
    joined = "".join(seq_lines)
    no_spaces = "".join(joined.split())
    no_digits = "".join(ch for ch in no_spaces if not ch.isdigit())
    return no_digits.upper()


def validate_sequence(sequence: str):
    """Return (is_valid: bool, message: str). Only A, T, G, C are allowed."""
    if not sequence:
        return False, "Sequence is empty. Please paste or upload a DNA sequence."
    invalid_chars = sorted(set(sequence) - VALID_BASES)
    if invalid_chars:
        return (
            False,
            f"Invalid characters found: {', '.join(invalid_chars)}. "
            f"Only A, T, G, and C are allowed.",
        )
    return True, "Sequence is valid."


def calculate_statistics(sequence: str) -> dict:
    """Compute nucleotide composition statistics for a cleaned DNA sequence."""
    seq = sequence.upper()
    length = len(seq)
    a_count = seq.count("A")
    t_count = seq.count("T")
    g_count = seq.count("G")
    c_count = seq.count("C")

    def pct(n):
        return round((n / length) * 100, 2) if length else 0.0

    gc_percentage = pct(g_count + c_count)
    at_percentage = pct(a_count + t_count)

    return {
        "length": length,
        "a_count": a_count,
        "t_count": t_count,
        "g_count": g_count,
        "c_count": c_count,
        "a_pct": pct(a_count),
        "t_pct": pct(t_count),
        "g_pct": pct(g_count),
        "c_pct": pct(c_count),
        "gc_pct": gc_percentage,
        "at_pct": at_percentage,
    }


def parse_fasta_or_raw(text: str):
    """Parse pasted text that may or may not be FASTA-formatted.
    Returns (header, cleaned_sequence)."""
    if not text:
        return "Unnamed sequence", ""
    stripped = text.strip()
    if stripped.startswith(">"):
        lines = stripped.splitlines()
        header = lines[0][1:].strip() or "Unnamed sequence"
    else:
        header = "Pasted sequence"
    return header, clean_sequence(text)


# ==============================================================================
# CHAOS GAME REPRESENTATION (CGR)
# ==============================================================================
def generate_cgr(sequence: str, point_size: float = 1.2):
    """
    Real Chaos Game Representation.

    Coordinates:
        A = (0, 0)
        T = (1, 0)
        G = (0, 1)
        C = (1, 1)

    Starting at the center (0.5, 0.5), each nucleotide moves the current
    point halfway toward its corresponding corner. This is a deterministic
    mapping — identical sequences always produce identical fingerprints,
    and it is NOT a random or decorative image.
    """
    x, y = 0.5, 0.5
    xs, ys = [], []
    for base in sequence.upper():
        if base not in CGR_CORNERS:
            continue
        cx, cy = CGR_CORNERS[base]
        x, y = (x + cx) / 2.0, (y + cy) / 2.0
        xs.append(x)
        ys.append(y)

    fig, ax = plt.subplots(figsize=(5, 5), facecolor="#0b0f14")
    ax.set_facecolor("#0b0f14")
    ax.scatter(xs, ys, s=point_size, c="#2dd4bf", alpha=0.6, edgecolors="none")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    corner_style = dict(fontsize=12, color="#e6f1ef", fontweight="bold")
    ax.text(0.0, -0.04, "A", ha="center", va="top", transform=ax.transAxes, **corner_style)
    ax.text(1.0, -0.04, "T", ha="center", va="top", transform=ax.transAxes, **corner_style)
    ax.text(0.0, 1.03, "G", ha="center", va="bottom", transform=ax.transAxes, **corner_style)
    ax.text(1.0, 1.03, "C", ha="center", va="bottom", transform=ax.transAxes, **corner_style)
    fig.tight_layout()
    return fig


def fig_to_png_bytes(fig) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=200, bbox_inches="tight", facecolor=fig.get_facecolor())
    buf.seek(0)
    return buf.getvalue()


# ==============================================================================
# QR CODE GENERATION
# ==============================================================================
def generate_genetrace_id(sequence: str) -> str:
    """Deterministic unique ID: SHA-256 hash of the DNA sequence (first 16 hex chars)."""
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:16]


def generate_qr(sequence: str, stats: dict):
    """
    Generate a real QR code encoding the DNA sequence plus metadata as JSON.

    Returns (qr_image_bytes: bytes, metadata_json_str: str).
    """
    genetrace_id = generate_genetrace_id(sequence)
    metadata = {
        "GeneTrace_ID": genetrace_id,
        "DNA_sequence": sequence,
        "length": stats["length"],
        "GC_percentage": stats["gc_pct"],
        "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
    }
    payload = json.dumps(metadata)

    qr = qrcode.QRCode(
        version=None,
        error_correction=qrcode.constants.ERROR_CORRECT_M,
        box_size=8,
        border=4,
    )
    qr.add_data(payload)
    qr.make(fit=True)
    img = qr.make_image(fill_color="#0f766e", back_color="#f0fdfa").convert("RGB")

    buf = io.BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return buf.getvalue(), payload


# ==============================================================================
# NCBI BLAST (REAL — no fake or simulated results)
# ==============================================================================
def run_ncbi_blast(sequence: str, email: str, program: str = "blastn", database: str = "nt"):
    """
    Submit a real BLAST search to NCBI and return the parsed Bio.Blast.Record.
    Raises exceptions on failure — the caller is responsible for catching them
    and showing a friendly message. Never fabricates results.
    """
    if not email or not email.strip():
        raise ValueError("Please enter your NCBI email before searching.")

    from Bio import Entrez
    Entrez.email = email.strip()

    old_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(180)
    try:
        result_handle = NCBIWWW.qblast(program, database, sequence)
        blast_record = NCBIXML.read(result_handle)
        result_handle.close()
        return blast_record
    finally:
        socket.setdefaulttimeout(old_timeout)


def parse_blast_results(blast_record, top_n: int = 10) -> pd.DataFrame:
    """Convert a Bio.Blast.Record into a clean results DataFrame (top N hits)."""
    rows = []
    if not blast_record or not blast_record.alignments:
        return pd.DataFrame(columns=[
            "Rank", "Accession", "Description", "Identity %",
            "Alignment Length", "E-value", "Bit Score",
        ])

    for rank, alignment in enumerate(blast_record.alignments[:top_n], start=1):
        if not alignment.hsps:
            continue
        best_hsp = alignment.hsps[0]
        identity_pct = round((best_hsp.identities / best_hsp.align_length) * 100, 2) if best_hsp.align_length else 0.0
        rows.append({
            "Rank": rank,
            "Accession": getattr(alignment, "accession", None) or alignment.hit_id,
            "Description": alignment.hit_def,
            "Identity %": identity_pct,
            "Alignment Length": best_hsp.align_length,
            "E-value": best_hsp.expect,
            "Bit Score": round(best_hsp.bits, 2),
        })
    return pd.DataFrame(rows)


# ==============================================================================
# PAIRWISE ALIGNMENT (REAL, via Biopython)
# ==============================================================================
def _get_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    return aligner


def perform_alignment(seq1: str, seq2: str):
    """Run a real global pairwise alignment. Returns (aligned1, aligned2, score)."""
    if not seq1 or not seq2:
        return "", "", 0.0
    aligner = _get_aligner()
    alignments = aligner.align(seq1.upper(), seq2.upper())
    best = alignments[0]
    return str(best[0]), str(best[1]), float(best.score)


def calculate_identity(aligned1: str, aligned2: str):
    """
    Given two aligned (gapped) sequences of equal length, return a dict with:
    identity %, alignment length, mismatch count, gap count.
    """
    if not aligned1 or not aligned2 or len(aligned1) != len(aligned2):
        return {"identity_pct": 0.0, "alignment_length": 0, "mismatches": 0, "gaps": 0}

    alignment_length = len(aligned1)
    matches = mismatches = gaps = 0
    for a, b in zip(aligned1, aligned2):
        if a == "-" or b == "-":
            gaps += 1
        elif a == b:
            matches += 1
        else:
            mismatches += 1

    identity_pct = round((matches / alignment_length) * 100, 2) if alignment_length else 0.0
    return {
        "identity_pct": identity_pct,
        "alignment_length": alignment_length,
        "mismatches": mismatches,
        "gaps": gaps,
    }


# ==============================================================================
# SHARED UI HELPERS
# ==============================================================================
def sequence_input_widget(key_prefix: str, label: str = "Paste your DNA sequence here"):
    """Reusable paste-or-upload widget. Returns (header, cleaned_sequence)."""
    method = st.radio(
        "Input method", ["Paste sequence", "Upload FASTA/TXT file"],
        key=f"{key_prefix}_method", horizontal=True,
    )
    header, sequence = "", ""
    if method == "Paste sequence":
        text = st.text_area(
            label, height=160, key=f"{key_prefix}_text",
            placeholder=">Optional FASTA header\nATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGC",
        )
        if text and text.strip():
            header, sequence = parse_fasta_or_raw(text)
    else:
        uploaded = st.file_uploader(
            "Upload FASTA or TXT file", type=["fasta", "fa", "txt"], key=f"{key_prefix}_upload"
        )
        if uploaded is not None:
            content = uploaded.read().decode("utf-8", errors="ignore")
            header, sequence = parse_fasta_or_raw(content)
    return header, sequence


def show_cgr_side_by_side(seq1: str, label1: str, seq2: str, label2: str):
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"**{label1}**")
        fig1 = generate_cgr(seq1)
        st.pyplot(fig1, use_container_width=True)
        plt.close(fig1)
    with col2:
        st.markdown(f"**{label2}**")
        fig2 = generate_cgr(seq2)
        st.pyplot(fig2, use_container_width=True)
        plt.close(fig2)


def render_ncbi_status_badge():
    email_set = bool(st.session_state.get("ncbi_email", "").strip())
    if email_set:
        st.success(f"NCBI email configured: {st.session_state['ncbi_email']}")
    else:
        st.warning("NCBI email not set — configure it in the sidebar before searching NCBI.")


# ==============================================================================
# PAGE: DASHBOARD
# ==============================================================================
def render_dashboard():
    st.title("🧬 GeneTrace")
    st.caption("Interactive DNA Fingerprinting & Sequence Identification System")

    stats = st.session_state.get("current_stats")
    seq_length = stats["length"] if stats else 0
    gc_pct = stats["gc_pct"] if stats else 0.0
    ncbi_status = "Configured" if st.session_state.get("ncbi_email", "").strip() else "Not configured"

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Current Sequence Length", f"{seq_length} bp")
    c2.metric("GC Percentage", f"{gc_pct}%")
    c3.metric("Analyses Performed", st.session_state.get("analysis_count", 0))
    c4.metric("NCBI Status", ncbi_status)

    st.markdown("---")
    st.markdown(
        """
        <div class="gt-card">
        <b>GeneTrace</b> is an educational bioinformatics prototype that takes a raw DNA
        sequence through a full identification pipeline: cleaning &amp; validation,
        composition statistics, a deterministic visual fingerprint (Chaos Game
        Representation), a shareable QR code, and a real remote BLAST search against
        NCBI's nucleotide database.
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.subheader("Workflow")
    st.markdown(
        """
        ```
        DNA Sequence
             ↓
        Cleaning & Validation
             ↓
        DNA Statistics
             ↓
        CGR Fingerprint
             ↓
        QR Code
             ↓
        NCBI BLAST
             ↓
        Sequence Matching
             ↓
        Result Interpretation
        ```
        """
    )

    if stats:
        st.subheader("Current Sequence Snapshot")
        comp_df = pd.DataFrame({
            "Base": ["A", "T", "G", "C"],
            "Count": [stats["a_count"], stats["t_count"], stats["g_count"], stats["c_count"]],
        })
        fig, ax = plt.subplots(figsize=(5, 3), facecolor="#0b0f14")
        ax.set_facecolor("#0b0f14")
        colors = ["#2563eb", "#dc2626", "#16a34a", "#ca8a04"]
        ax.bar(comp_df["Base"], comp_df["Count"], color=colors)
        ax.tick_params(colors="#e6f1ef")
        for spine in ax.spines.values():
            spine.set_color("#16463f")
        st.pyplot(fig, use_container_width=True)
        plt.close(fig)
    else:
        st.info("No sequence analyzed yet. Go to **🧬 DNA Analyzer** to get started.")


# ==============================================================================
# PAGE: DNA ANALYZER
# ==============================================================================
def render_analyzer():
    st.title("🧬 DNA Analyzer")
    st.write("Paste a DNA sequence (or upload a FASTA/TXT file) to clean, validate, and analyze it.")

    header, raw_sequence = sequence_input_widget("analyzer", label="Paste your DNA sequence here")

    if st.button("🔬 Analyze DNA", type="primary", disabled=not raw_sequence):
        sequence = clean_sequence(raw_sequence)
        valid, message = validate_sequence(sequence)
        if not valid:
            st.error(message)
            return

        stats = calculate_statistics(sequence)
        st.session_state["current_sequence"] = sequence
        st.session_state["current_header"] = header or "Pasted sequence"
        st.session_state["current_stats"] = stats
        st.session_state["analysis_count"] += 1
        # Reset downstream artifacts tied to the previous sequence
        st.session_state["cgr_fig_bytes"] = None
        st.session_state["qr_image_bytes"] = None
        st.session_state["qr_metadata_json"] = None

    sequence = st.session_state.get("current_sequence", "")
    stats = st.session_state.get("current_stats")

    if not sequence or not stats:
        st.info("Analyze a sequence above to see statistics, CGR, and QR code.")
        return

    st.success(f"Analysis complete for: {st.session_state.get('current_header', 'Pasted sequence')}")

    st.subheader("Sequence Statistics")
    c1, c2, c3, c4, c5, c6 = st.columns(6)
    c1.metric("Length", f"{stats['length']} bp")
    c2.metric("A count / %", f"{stats['a_count']} / {stats['a_pct']}%")
    c3.metric("T count / %", f"{stats['t_count']} / {stats['t_pct']}%")
    c4.metric("G count / %", f"{stats['g_count']} / {stats['g_pct']}%")
    c5.metric("C count / %", f"{stats['c_count']} / {stats['c_pct']}%")
    c6.metric("GC % / AT %", f"{stats['gc_pct']}% / {stats['at_pct']}%")

    st.subheader("Nucleotide Composition")
    comp_df = pd.DataFrame({
        "Base": ["A", "T", "G", "C"],
        "Count": [stats["a_count"], stats["t_count"], stats["g_count"], stats["c_count"]],
    })
    fig, ax = plt.subplots(figsize=(6, 3.5), facecolor="#0b0f14")
    ax.set_facecolor("#0b0f14")
    colors = ["#2563eb", "#dc2626", "#16a34a", "#ca8a04"]
    ax.bar(comp_df["Base"], comp_df["Count"], color=colors)
    ax.set_ylabel("Count", color="#e6f1ef")
    ax.tick_params(colors="#e6f1ef")
    for spine in ax.spines.values():
        spine.set_color("#16463f")
    st.pyplot(fig, use_container_width=True)
    plt.close(fig)

    st.subheader("Cleaned DNA Sequence")
    st.text_area("Cleaned sequence", sequence, height=100, key="cleaned_seq_display")
    st.download_button(
        "⬇️ Download Clean DNA Sequence (TXT)",
        sequence,
        file_name="genetrace_clean_sequence.txt",
        mime="text/plain",
    )

    metadata = {
        "header": st.session_state.get("current_header", ""),
        "GeneTrace_ID": generate_genetrace_id(sequence),
        "length": stats["length"],
        "gc_percentage": stats["gc_pct"],
        "at_percentage": stats["at_pct"],
        "a_count": stats["a_count"], "t_count": stats["t_count"],
        "g_count": stats["g_count"], "c_count": stats["c_count"],
        "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
    }
    st.download_button(
        "⬇️ Download DNA Metadata (JSON)",
        json.dumps(metadata, indent=2),
        file_name="genetrace_metadata.json",
        mime="application/json",
    )

    # --- Chaos Game Representation ---
    st.subheader("🧬 Chaos Game Representation")
    fig_cgr = generate_cgr(sequence)
    st.pyplot(fig_cgr, use_container_width=True)
    cgr_bytes = fig_to_png_bytes(fig_cgr)
    st.session_state["cgr_fig_bytes"] = cgr_bytes
    plt.close(fig_cgr)
    st.caption(
        "CGR maps each nucleotide to a movement halfway toward a fixed corner of a unit "
        "square (A, T, G, C). Starting from the center, the resulting point cloud is a "
        "deterministic fingerprint: identical sequences always produce identical patterns, "
        "and similar sequence composition/order tends to produce visually related patterns."
    )
    st.download_button(
        "⬇️ Download CGR Image (PNG)",
        cgr_bytes,
        file_name="genetrace_cgr.png",
        mime="image/png",
    )

    # --- QR Code ---
    st.subheader("QR Code")
    qr_bytes, qr_json = generate_qr(sequence, stats)
    st.session_state["qr_image_bytes"] = qr_bytes
    st.session_state["qr_metadata_json"] = qr_json
    col1, col2 = st.columns([1, 2])
    with col1:
        st.image(qr_bytes, caption="GeneTrace QR Code", width=240)
        st.download_button(
            "⬇️ Download QR Code (PNG)",
            qr_bytes,
            file_name="genetrace_qr.png",
            mime="image/png",
        )
    with col2:
        st.write("**JSON data encoded inside the QR code:**")
        st.code(qr_json, language="json")


# ==============================================================================
# PAGE: NCBI IDENTIFICATION (REAL BLAST)
# ==============================================================================
def render_ncbi_identification():
    st.title("🔎 NCBI Identification")
    st.write(
        "Submit your DNA sequence to NCBI's real BLAST service (`blastn` against the `nt` "
        "nucleotide database) and view real matches. No results are ever simulated or "
        "hardcoded."
    )

    render_ncbi_status_badge()

    sequence = st.session_state.get("current_sequence", "")
    if not sequence:
        st.info("No analyzed sequence found. You can analyze one below, or go to the DNA Analyzer page.")
        header, raw_sequence = sequence_input_widget("ncbi_direct", label="Paste your DNA sequence here")
        if raw_sequence:
            candidate = clean_sequence(raw_sequence)
            valid, msg = validate_sequence(candidate)
            if valid:
                sequence = candidate
            else:
                st.error(msg)
    else:
        current_header = st.session_state.get("current_header", "Unknown sequence")
        cap_col, btn_col = st.columns([4, 1])
        with cap_col:
            st.caption(
                f"Using currently analyzed sequence: **{current_header}** ({len(sequence)} bp). "
                "You can re-analyze a different sequence on the DNA Analyzer page."
            )
        with btn_col:
            if st.button("🔄 Clear sequence", key="ncbi_clear_sequence"):
                st.session_state["current_sequence"] = ""
                st.session_state["current_header"] = ""
                st.session_state["current_stats"] = None
                st.session_state["blast_results_df"] = None
                st.session_state["blast_best_match"] = None
                st.session_state["blast_raw_status"] = "No BLAST search performed yet."
                st.rerun()

    st.info(
        "⏳ A real BLAST search against NCBI's `nt` database typically takes anywhere from "
        "about 30 seconds to a few minutes, depending on NCBI server load."
    )

    search_disabled = not sequence
    if st.button("🔎 Search NCBI BLAST", type="primary", disabled=search_disabled):
        email = st.session_state.get("ncbi_email", "").strip()
        if not email:
            st.error("Please enter your NCBI email before searching.")
        else:
            valid, msg = validate_sequence(sequence)
            if not valid:
                st.error(msg)
            else:
                try:
                    with st.spinner("Submitting sequence to NCBI BLAST — this can take a while..."):
                        blast_record = run_ncbi_blast(sequence, email)
                    results_df = parse_blast_results(blast_record, top_n=10)
                    st.session_state["blast_results_df"] = results_df
                    if results_df.empty:
                        st.session_state["blast_best_match"] = None
                        st.session_state["blast_raw_status"] = "NCBI BLAST completed but returned no matches."
                    else:
                        st.session_state["blast_best_match"] = results_df.iloc[0].to_dict()
                        st.session_state["blast_raw_status"] = "NCBI BLAST completed successfully."
                except ValueError as e:
                    st.error(str(e))
                except socket.timeout:
                    st.error(
                        "NCBI BLAST could not be reached. Please check your internet "
                        "connection and try again."
                    )
                except (OSError, ConnectionError):
                    st.error(
                        "NCBI BLAST could not be reached. Please check your internet "
                        "connection and try again."
                    )
                except Exception as e:
                    st.error(
                        "NCBI BLAST could not be reached. Please check your internet "
                        f"connection and try again. (Details: {e})"
                    )

    results_df = st.session_state.get("blast_results_df")
    best_match = st.session_state.get("blast_best_match")

    if results_df is not None:
        if results_df.empty:
            st.warning("NCBI BLAST completed but no matches were found for this sequence.")
        else:
            st.subheader("🏆 Best NCBI Match")
            st.markdown(
                f"""
                <div class="gt-best-match">
                <b>Description / Organism:</b> {best_match['Description']}<br>
                <b>Accession:</b> {best_match['Accession']}<br>
                <b>Identity %:</b> {best_match['Identity %']}%<br>
                <b>Alignment Length:</b> {best_match['Alignment Length']}<br>
                <b>E-value:</b> {best_match['E-value']}<br>
                <b>Bit Score:</b> {best_match['Bit Score']}
                </div>
                """,
                unsafe_allow_html=True,
            )

            st.subheader("Top NCBI BLAST Matches")
            st.dataframe(results_df, use_container_width=True)

            csv_bytes = results_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                "⬇️ Download NCBI Results (CSV)",
                csv_bytes,
                file_name="genetrace_ncbi_results.csv",
                mime="text/csv",
            )

    st.caption(st.session_state.get("blast_raw_status", ""))
    st.info(
        "A BLAST match reflects sequence similarity to entries in NCBI's public database. "
        "It is not a definitive clinical diagnosis, legal forensic evidence, or definitive "
        "species identification."
    )


# ==============================================================================
# PAGE: COMPARE SEQUENCES
# ==============================================================================
def render_comparison():
    st.title("⚖️ Compare Sequences")
    st.write("Enter two DNA sequences to clean, validate, align, and compare them.")

    col1, col2 = st.columns(2)
    with col1:
        raw1 = st.text_area("Sequence 1", height=140, key="cmp_raw1",
                             placeholder=">Optional header\nATGCGTACG...")
    with col2:
        raw2 = st.text_area("Sequence 2", height=140, key="cmp_raw2",
                             placeholder=">Optional header\nATGCGTACG...")

    if st.button("⚖️ Compare Sequences", type="primary", disabled=not (raw1 and raw2)):
        header1, seq1 = parse_fasta_or_raw(raw1)
        header2, seq2 = parse_fasta_or_raw(raw2)

        valid1, msg1 = validate_sequence(seq1)
        valid2, msg2 = validate_sequence(seq2)
        if not valid1:
            st.error(f"Sequence 1: {msg1}")
            return
        if not valid2:
            st.error(f"Sequence 2: {msg2}")
            return

        stats1 = calculate_statistics(seq1)
        stats2 = calculate_statistics(seq2)
        aligned1, aligned2, score = perform_alignment(seq1, seq2)
        identity_info = calculate_identity(aligned1, aligned2)

        st.session_state["cmp_seq1"] = seq1
        st.session_state["cmp_seq2"] = seq2
        st.session_state["cmp_header1"] = header1
        st.session_state["cmp_header2"] = header2
        st.session_state["cmp_stats1"] = stats1
        st.session_state["cmp_stats2"] = stats2
        st.session_state["cmp_aligned1"] = aligned1
        st.session_state["cmp_aligned2"] = aligned2
        st.session_state["cmp_score"] = score
        st.session_state["cmp_identity_info"] = identity_info

        report_lines = [
            "GeneTrace — Sequence Comparison Report",
            f"Generated: {datetime.datetime.now().isoformat(timespec='seconds')}",
            "",
            f"Sequence 1 ({header1}): length={stats1['length']} bp, GC%={stats1['gc_pct']}",
            f"Sequence 2 ({header2}): length={stats2['length']} bp, GC%={stats2['gc_pct']}",
            "",
            f"Alignment score: {score}",
            f"Identity %: {identity_info['identity_pct']}",
            f"Alignment length: {identity_info['alignment_length']}",
            f"Mismatches: {identity_info['mismatches']}",
            f"Gaps: {identity_info['gaps']}",
            "",
            "Aligned Sequence 1:",
            aligned1,
            "Aligned Sequence 2:",
            aligned2,
        ]
        st.session_state["comparison_results_text"] = "\n".join(report_lines)

    if "cmp_seq1" in st.session_state and "cmp_seq2" in st.session_state:
        seq1, seq2 = st.session_state["cmp_seq1"], st.session_state["cmp_seq2"]
        stats1, stats2 = st.session_state["cmp_stats1"], st.session_state["cmp_stats2"]
        aligned1, aligned2 = st.session_state["cmp_aligned1"], st.session_state["cmp_aligned2"]
        score = st.session_state["cmp_score"]
        identity_info = st.session_state["cmp_identity_info"]

        st.subheader("Composition")
        c1, c2 = st.columns(2)
        with c1:
            st.markdown(f"**{st.session_state['cmp_header1']}**")
            st.write(f"Length: {stats1['length']} bp")
            st.write(f"A/T/G/C: {stats1['a_count']}/{stats1['t_count']}/{stats1['g_count']}/{stats1['c_count']}")
            st.write(f"GC%: {stats1['gc_pct']}%  •  AT%: {stats1['at_pct']}%")
        with c2:
            st.markdown(f"**{st.session_state['cmp_header2']}**")
            st.write(f"Length: {stats2['length']} bp")
            st.write(f"A/T/G/C: {stats2['a_count']}/{stats2['t_count']}/{stats2['g_count']}/{stats2['c_count']}")
            st.write(f"GC%: {stats2['gc_pct']}%  •  AT%: {stats2['at_pct']}%")

        st.subheader("Global Pairwise Alignment")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Alignment Score", score)
        c2.metric("Identity %", f"{identity_info['identity_pct']}%")
        c3.metric("Mismatches", identity_info["mismatches"])
        c4.metric("Gaps", identity_info["gaps"])
        st.caption(f"Alignment length: {identity_info['alignment_length']}")

        with st.expander("View full alignment", expanded=True):
            st.text(aligned1)
            match_line = "".join(
                "|" if a == b and a != "-" else " " for a, b in zip(aligned1, aligned2)
            )
            st.text(match_line)
            st.text(aligned2)

        st.subheader("Visual Fingerprint Comparison")
        show_cgr_side_by_side(seq1, st.session_state["cmp_header1"], seq2, st.session_state["cmp_header2"])

        if st.session_state.get("comparison_results_text"):
            st.download_button(
                "⬇️ Download Comparison Results (TXT)",
                st.session_state["comparison_results_text"],
                file_name="genetrace_comparison.txt",
                mime="text/plain",
            )


# ==============================================================================
# PAGE: MYSTERY DNA
# ==============================================================================
def render_mystery():
    st.title("🕵️ Mystery DNA")
    st.info(
        "**Educational challenge only.** Mystery DNA results are for learning purposes and "
        "are not forensic proof of anything."
    )

    choice = st.selectbox("Choose a mystery sample", list(MYSTERY_SEQUENCES.keys()), key="mystery_choice")
    mystery_seq = MYSTERY_SEQUENCES[choice]

    st.code(mystery_seq, language=None)
    st.caption(f"Length: {len(mystery_seq)} bp")

    if st.button("🔬 Analyze Mystery Sample", type="primary"):
        valid, msg = validate_sequence(mystery_seq)
        if not valid:
            st.error(msg)
            return
        stats = calculate_statistics(mystery_seq)
        st.session_state["current_sequence"] = mystery_seq
        st.session_state["current_header"] = choice
        st.session_state["current_stats"] = stats
        st.session_state["analysis_count"] += 1
        st.success(f"Loaded '{choice}' into the analyzer. Statistics computed below.")

        c1, c2, c3 = st.columns(3)
        c1.metric("Length", f"{stats['length']} bp")
        c2.metric("GC %", f"{stats['gc_pct']}%")
        c3.metric("AT %", f"{stats['at_pct']}%")

        st.subheader("🧬 CGR Fingerprint")
        fig = generate_cgr(mystery_seq)
        st.pyplot(fig, use_container_width=True)
        cgr_bytes = fig_to_png_bytes(fig)
        plt.close(fig)
        st.download_button(
            "⬇️ Download CGR Image (PNG)", cgr_bytes,
            file_name="mystery_cgr.png", mime="image/png", key="mystery_cgr_dl",
        )

        st.subheader("QR Code")
        qr_bytes, qr_json = generate_qr(mystery_seq, stats)
        col1, col2 = st.columns([1, 2])
        with col1:
            st.image(qr_bytes, width=220)
            st.download_button(
                "⬇️ Download QR Code (PNG)", qr_bytes,
                file_name="mystery_qr.png", mime="image/png", key="mystery_qr_dl",
            )
        with col2:
            st.code(qr_json, language="json")

        st.info(
            "To search NCBI BLAST for this mystery sample, head to **🔎 NCBI Identification** — "
            "it will automatically use this sequence since it's now the current analyzed sequence."
        )


# ==============================================================================
# PAGE: ABOUT
# ==============================================================================
def render_about():
    st.title("ℹ️ About GeneTrace")
    st.markdown(
        """
        <div class="gt-card">

        **GeneTrace** is a single-file educational bioinformatics prototype built with
        Streamlit and Biopython. It walks a DNA sequence through a complete
        identification pipeline.

        #### DNA Sequence Analysis
        Raw pasted or uploaded text is cleaned (headers, whitespace, and digits removed),
        validated to contain only A, T, G, C, and reduced to composition statistics
        (base counts, GC%, AT%).

        #### Chaos Game Representation (CGR)
        A deterministic algorithm that maps each nucleotide to a step halfway toward a
        fixed corner of a unit square (A, T, G, C corners). The resulting point cloud is
        a visual fingerprint that depends entirely on sequence composition and order —
        not randomness or decoration.

        #### QR Codes
        A QR code is generated containing the actual DNA sequence plus metadata (a
        SHA-256-derived GeneTrace ID, sequence length, GC percentage, and timestamp) as
        JSON, so the sequence and its key statistics can be shared or scanned.

        #### BLAST & NCBI
        BLAST (Basic Local Alignment Search Tool) compares a query sequence against
        NCBI's public sequence databases to find similar sequences. GeneTrace submits
        **real** `blastn` searches against the `nt` nucleotide database via Biopython's
        `Bio.Blast.NCBIWWW` and parses **real** results via `Bio.Blast.NCBIXML` — nothing
        is hardcoded or simulated.

        #### Sequence Identity, E-value, and Alignment
        - **Identity %** — the percentage of aligned positions that match exactly.
        - **E-value** — the expected number of alignments with a similar score that would
          occur by chance in a database this size; lower is more significant.
        - **Alignment** — the position-by-position pairing of two sequences (with gaps)
          that maximizes similarity under a scoring scheme.

        #### Limitations
        GeneTrace is a teaching and research tool. A BLAST match or comparison result is
        a similarity-based estimate, not proof of identity, ancestry, or diagnosis.

        </div>

        <div class="gt-card" style="border-color:#f59e0b;">
        <b>Disclaimer:</b> GeneTrace is an educational and research prototype. Results
        should not be interpreted as clinical diagnoses, legal forensic evidence, or
        definitive species identification.
        </div>
        """,
        unsafe_allow_html=True,
    )


# ==============================================================================
# MAIN APP
# ==============================================================================
def main():
    init_session_state()

    with st.sidebar:
        st.markdown("## 🧬 GeneTrace")
        st.caption("DNA Fingerprinting & Identification")
        page = st.radio(
            "Navigate",
            [
                "🏠 Dashboard",
                "🧬 DNA Analyzer",
                "🔎 NCBI Identification",
                "⚖️ Compare Sequences",
                "🕵️ Mystery DNA",
                "ℹ️ About",
            ],
            label_visibility="collapsed",
        )

        st.markdown("---")
        st.markdown("### NCBI Settings")
        st.text_input(
            "NCBI Email", key="ncbi_email", placeholder="you@example.com",
            help="Required by NCBI for every BLAST/Entrez request.",
        )
        st.text_input(
            "NCBI API Key (optional)", key="ncbi_api_key", type="password",
            help="Optional. Not required for qblast, but useful for other Entrez requests.",
        )

    pages = {
        "🏠 Dashboard": render_dashboard,
        "🧬 DNA Analyzer": render_analyzer,
        "🔎 NCBI Identification": render_ncbi_identification,
        "⚖️ Compare Sequences": render_comparison,
        "🕵️ Mystery DNA": render_mystery,
        "ℹ️ About": render_about,
    }
    pages[page]()


if __name__ == "__main__":
    main()
