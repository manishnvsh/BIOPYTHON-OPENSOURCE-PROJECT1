# == Imports and Streamlit setup ==
import streamlit as st; from Bio import SeqIO; from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio import pairwise2; from Bio.pairwise2 import format_alignment
import io; import pandas as pd; import matplotlib.pyplot as plt; import base64; import textwrap; import traceback

# == Page configuration and main title ==
st.set_page_config(page_title="Biopython Toolkit", layout="wide")
st.title("Biopython Toolkit — Streamlit App")
st.markdown("""A polished toolbox for common sequence tasks using **Biopython**.
Upload FASTA files, inspect sequences, compute GC content, translate, find ORFs, run pairwise alignment, search motifs, and visualize codon usage and GC distribution.
""")
st.sidebar.header("Quick actions"); st.sidebar.markdown("Upload a multi-FASTA and select tools from the main panel."); st.sidebar.markdown("---")

# == Helper: Parse FASTA files and create download links ==
def parse_fasta_bytes(uploaded_file):
    try:
        if uploaded_file is None: return []
        content = uploaded_file.read().decode('utf-8'); handle = io.StringIO(content)
        return list(SeqIO.parse(handle, "fasta"))
    except Exception as e: st.error("Failed to parse FASTA. Make sure it's a valid FASTA file. Error: " + str(e)); return []

def download_link(data: bytes, fname: str, text: str):
    b64 = base64.b64encode(data).decode(); return f'<a href="data:file/txt;base64,{b64}" download="{fname}">{text}</a>'

uploaded = st.file_uploader("Upload FASTA file (multi-FASTA allowed)", type=["fa","fasta","txt"])

# == Example FASTA display for quick test ==
with st.expander("Example FASTA"):
    st.code(textwrap.dedent(""">seq1_example
ATGCGTACGTTAGCGTAGCTAGCTAGCGTAGCTAGCTGACTGATCGATCGATCGTAGCTAG
>seq2_example
ATGAAATTTGGGCCCTTTAAACCCGGGATGCTAGCTAGCTAA"""))
records = parse_fasta_bytes(uploaded) if uploaded else []

# == Toolbox control panel ==
st.header("Tools"); col1, col2 = st.columns(2)
with col1:
    show_upload = st.checkbox("Show uploaded sequences", True)
    show_summary = st.checkbox("Sequence summary (length, GC%)", True)
    show_revtrans = st.checkbox("Reverse complement & translation", True)
with col2:
    show_orf = st.checkbox("Six-frame translation & ORF finder", True)
    show_align = st.checkbox("Pairwise alignment (global/local)", True)
    show_motif = st.checkbox("Motif search", True)
    show_codon = st.checkbox("Codon usage & plots", True)

# == Uploaded sequence block ==
if show_upload:
    st.subheader("Uploaded sequences")
    if not records: st.info("No FASTA uploaded yet. Paste or upload a FASTA to begin. Use the example above to test.")
    else:
        for rec in records:
            with st.expander(f"{rec.id} | {len(rec.seq)} bp"):
                st.write(f"Description: {rec.description}")
                st.code(str(rec.seq[:1000]) + ("..." if len(rec.seq)>1000 else ""))
                st.markdown(download_link(f">{rec.id}\n{str(rec.seq)}\n".encode('utf-8'), f"{rec.id}.fasta", "Download this sequence (FASTA)"), unsafe_allow_html=True)

# == Sequence summary block: length, GC ==
if show_summary:
    st.subheader("Sequence summary")
    if not records: st.info("Upload FASTA to compute summary.")
    else:
        rows = []
        for rec in records: rows.append({"ID": rec.id, "Length": len(rec.seq), "GC%": round(gc_fraction(rec.seq)*100,3)})
        df = pd.DataFrame(rows); st.dataframe(df); st.markdown("Download summary:")
        st.markdown(download_link(df.to_csv(index=False).encode('utf-8'), "summary.csv", "Download CSV"), unsafe_allow_html=True)

# == Reverse Complement & Translation UI ==
if show_revtrans:
    st.subheader("Reverse Complement & Translation")
    seq_choice = st.selectbox("Choose sequence", [r.id for r in records] if records else [], key="rev_choice")
    table = st.checkbox("Show translation table details", False)
    if seq_choice:
        rec = next(r for r in records if r.id==seq_choice); seq = rec.seq
        try:
            st.markdown("**Reverse Complement**"); st.code(str(seq.reverse_complement()))
            st.markdown("**Translation (standard table)**"); prot = seq.translate(to_stop=False)
            st.code(str(prot)); st.markdown("**Protein molecular weight (Da)**")
            try: st.write(round(molecular_weight(prot),3))
            except Exception: st.write("Unable to compute molecular weight for translation.")
            if table: st.write("Translation info: (showing first 200 aa)"); st.code(str(prot[:200]))
        except Exception as e: st.error("Error computing translation or reverse complement: " + str(e))

# == Six-frame ORF finder ==
if show_orf:
    st.subheader("Six-frame translation and ORF finder")
    seq_choice = st.selectbox("Choose sequence for ORF", [r.id for r in records] if records else [], key="orf_choice2")
    min_orf_len = st.number_input("Minimum ORF length (aa)", 10, 10000, 30)
    search_start_meth = st.selectbox("ORF search method", ["Start with M and end with Stop", "Any open region (no internal stop)"], key="orf_method")
    if seq_choice:
        rec = next(r for r in records if r.id==seq_choice); seq = rec.seq; frames = []
        def find_orfs_simple(s, strand, frame):
            prot = s.translate(to_stop=False); orfs = []
            if search_start_meth=="Start with M and end with Stop":
                for i, aa in enumerate(prot):
                    if aa == "M":
                        for j in range(i+1, len(prot)):
                            if prot[j] == "*":
                                length = j - i + 1
                                if length >= min_orf_len: 
                                    orfs.append({"strand": strand, "frame": frame, "aa_start": i, "aa_end": j, "length_aa": length,
                                                 "protein": str(prot[i:j+1]), "nuc_start": (i*3)+frame, "nuc_end": (j*3)+frame+3})
                                break
            else:
                start = None
                for i, aa in enumerate(prot):
                    if aa != "*" and start is None: start = i
                    if aa == "*" and start is not None:
                        length = i - start
                        if length >= min_orf_len:
                            orfs.append({"strand": strand, "frame": frame, "aa_start": start, "aa_end": i-1, "length_aa": length,
                                "protein": str(prot[start:i]), "nuc_start": (start*3)+frame, "nuc_end": (i*3)+frame})
                        start = None
            return orfs
        seq_str = seq
        for strand, s in [(+1, seq_str), (-1, seq_str.reverse_complement())]:
            for frame in range(3): s_frame = s[frame:]; frames += find_orfs_simple(s_frame, strand, frame)
        if frames:
            df_orfs = pd.DataFrame(frames).sort_values("length_aa", ascending=False)
            st.dataframe(df_orfs)
            st.markdown(download_link(df_orfs.to_csv(index=False).encode('utf-8'), "orfs.csv", "Download ORFs CSV"), unsafe_allow_html=True)
        else: st.write("No ORFs found meeting the length threshold. Try lowering the threshold.")

# == Pairwise alignment tool. Fixed indentation error ==
if show_align:
    st.subheader("Pairwise alignment")
    seqs = [r.id for r in records]
    if len(seqs) < 2: st.info("Upload FASTA with at least two sequences to align.")
    else:
        s1 = st.selectbox("Sequence 1", seqs, key="s1"); s2 = st.selectbox("Sequence 2", seqs, key="s2")
        method = st.radio("Method", ["global", "local"], key="align_method")
        rec1 = next(r for r in records if r.id==s1); rec2 = next(r for r in records if r.id==s2)
        if st.button("Run alignment", key="align_btn"):
            try:
                alns = pairwise2.align.globalxx(str(rec1.seq), str(rec2.seq), one_alignment_only=True) if method=="global" else pairwise2.align.localxx(str(rec1.seq), str(rec2.seq), one_alignment_only=True)
                if alns:
                    aln = alns[0]
                    st.text(format_alignment(*aln))
                    aln_text = format_alignment(*aln)
                    st.markdown(download_link(aln_text.encode('utf-8'), f"alignment_{s1}_vs_{s2}.txt", "Download alignment"), unsafe_allow_html=True)
                else: st.write("No alignment found.")
            except Exception as e: st.error("Alignment failed: " + str(e)); st.text(traceback.format_exc())

# == Motif (restriction/regex) search ==
if show_motif:
    st.subheader("Motif / Restriction search")
    motif = st.text_input("Enter motif (DNA regex, e.g., 'GAATTC' for EcoRI)")
    seq_choice = st.selectbox("Choose sequence", [r.id for r in records] if records else [], key="motif_seq")
    if motif and seq_choice:
        rec = next(r for r in records if r.id==seq_choice); seq = str(rec.seq); import re
        try:
            matches = [(m.start()+1, m.group(0)) for m in re.finditer(motif, seq)]
            st.write(f"Occurrences: {len(matches)}")
            if matches:
                dfm = pd.DataFrame(matches, columns=["start","motif"])
                st.dataframe(dfm)
                st.markdown(download_link(dfm.to_csv(index=False).encode('utf-8'), "motif_hits.csv", "Download motif hits"), unsafe_allow_html=True)
        except re.error as e: st.error("Invalid regex provided: " + str(e))

# == Codon usage and GC distribution plotting ==
if show_codon:
    st.subheader("Codon usage & GC distribution plots")
    seq_choice = st.selectbox("Choose sequence", [r.id for r in records] if records else [], key="codon_seq")
    if seq_choice:
        rec = next(r for r in records if r.id==seq_choice); seq = str(rec.seq)
        codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]; codon_counts = {}
        for c in codons: codon_counts[c] = codon_counts.get(c,0)+1 if len(c)==3 else codon_counts.get(c,0)
        df_codon = pd.DataFrame(list(codon_counts.items()), columns=["Codon","Count"]).sort_values("Count", ascending=False)
        st.dataframe(df_codon)
        st.markdown(download_link(df_codon.to_csv(index=False).encode('utf-8'), "codon_usage.csv", "Download codon table"), unsafe_allow_html=True)
        fig, ax = plt.subplots(figsize=(8,4)); ax.bar(df_codon['Codon'], df_codon['Count'])
        ax.set_xlabel("Codon"); ax.set_ylabel("Count"); ax.set_title("Codon usage (frame 0)")
        plt.xticks(rotation=90); st.pyplot(fig)
        window = st.slider("GC sliding window (bp)", 10, 200, 50); seq_upper = seq.upper(); gc_values, positions = [], []
        for i in range(0, max(1, len(seq_upper)-window+1), window//2):
            win = seq_upper[i:i+window]; positions.append(i+1); gc_values.append(gc_fraction(win)*100)
        fig2, ax2 = plt.subplots(figsize=(8,3)); ax2.plot(positions, gc_values)
        ax2.set_xlabel("Position (bp)"); ax2.set_ylabel("GC%"); ax2.set_title("Sliding-window GC%"); st.pyplot(fig2)

# == Sidebar information and developer export block ==
st.sidebar.markdown("---"); st.sidebar.markdown("Developed with Biopython & Streamlit. No external network calls are made.")
st.sidebar.markdown("You can download the full app package from the main menu (Developer export).")

# == Developer export tool ==
import zipfile, os; from io import BytesIO
def make_zip_bytes(files):
    buf = BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for fpath, arcname in files:
            zf.writestr(arcname, open(fpath, "rb").read())
    return buf.getvalue()

if st.sidebar.button("Create app zip (for download)"):
    try:
        files = [
            ("/mnt/data/streamlit_biopython_app.py", "streamlit_biopython_app.py"),
            ("/mnt/data/requirements.txt", "requirements.txt"),
            ("/mnt/data/README_streamlit_app.txt", "README_streamlit_app.txt")
        ]
        zip_bytes = make_zip_bytes(files)
        st.markdown(download_link(zip_bytes, "biopython_streamlit_app.zip", "Download app ZIP"), unsafe_allow_html=True)
    except Exception as e: st.error("Failed to create zip: " + str(e))

# End of application
