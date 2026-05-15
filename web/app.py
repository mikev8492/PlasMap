#!/usr/bin/env python3
"""
app.py — Streamlit web interface for PlasMap
=============================================

Wraps the PlasMap CLI tool (src/main.py) with a browser-based UI so users can:
  - Upload a plasmid sequence file (FASTA or GenBank)
  - Select restriction enzymes via a multiselect widget
  - Choose single-stranded or double-stranded linear map output
  - View the circular and linear maps inline
  - Download the maps and CSV results

All analysis logic lives in src/motif_id_lib/ and is unchanged.
This file only handles the UI layer.

Usage
-----
    # From the repo root:
    streamlit run web/app.py

Dependencies (add to environment.yml or pip install):
    streamlit >= 1.35
    All existing PlasMap dependencies remain required.
"""

import io
import sys
import csv
import tempfile
from pathlib import Path

import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------------------------------------
# Path setup — allow imports from src/ without installing the package
# ---------------------------------------------------------------------------
SRC_DIR = Path(__file__).parent.parent / "src"
sys.path.insert(0, str(SRC_DIR))

from motif_id_lib.input import Sequence, Enzymes
from motif_id_lib.motif_locator import Motifs
from motif_id_lib.output import PlasmidMap

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Path to the bundled restriction enzyme database CSV.
DB_FILE = SRC_DIR / "database" / "enzymes.csv"

#: File extensions PlasMap accepts as input.
GENBANK_EXTS = {"gb", "gbk"}
FASTA_EXTS   = {"fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn"}
VALID_EXTS   = GENBANK_EXTS | FASTA_EXTS

#: Default enzyme panel shown when the app first loads.
DEFAULT_ENZYMES = [
    "EcoRI",  "HindIII", "BamHI",  "XhoI",   "NotI",
    "SalI",   "PstI",    "KpnI",   "XbaI",   "EcoRV",
    "SmaI",   "NdeI",    "SacI",   "SpeI",   "BglII",
    "ApaI",   "SphI",    "MluI",   "ClaI",   "HaeIII",
    "Eco91I", "Eco24I",
]

#: Demo plasmid files bundled with the repository.
DEMO_FILES = {
    "pUC19 (FASTA)"         : SRC_DIR.parent / "inputs/test/pUC19.fasta",
    "pUC19 (GenBank)"       : SRC_DIR.parent / "inputs/test/pUC19.gb",
    "pCMV-GLuc (GenBank)"   : SRC_DIR.parent / "inputs/test/pCMV-GLuc.gb",
    "M13mp18 (GenBank)"     : SRC_DIR.parent / "inputs/test/M13mp18.gb",
    "pBeloBAC11 (GenBank)"  : SRC_DIR.parent / "inputs/test/pBeloBAC11.gb",
}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
@st.cache_data
def load_all_enzyme_names(db_path: Path) -> list[str]:
    """
    Read the enzyme database CSV and return a sorted list of all enzyme names.
    @st.cache_data: Decorator that prevents re-running the function when it is called with the same parameters. It prevents re-reading the enzymes database everytime the app state changes.
    
    Parameters
    ----------
    db_path : Path
        Path to src/database/enzymes.csv.

    Returns
    -------
    list[str]
        All enzyme names present in the database, alphabetically sorted.
    """
    names = []
    with open(db_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            names.append(row["enzyme"])
    return sorted(names)


def save_uploaded_file(uploaded_file) -> Path:
    """
    Write a Streamlit UploadedFile to a NamedTemporaryFile on disk and return
    its path.  The caller is responsible for deleting the file when done.

    PlasMap's Sequence class needs a real file path (it opens the file itself),
    so we cannot pass the in-memory buffer directly.

    Parameters
    ----------
    uploaded_file : streamlit.runtime.uploaded_file_manager.UploadedFile
        The object returned by st.file_uploader().

    Returns
    -------
    Path
        Absolute path to the temporary file.  The suffix preserves the
        original file extension so the format-detection logic in input.py works.
    """
    suffix = Path(uploaded_file.name).suffix  # e.g. ".fasta" or ".gb"
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    tmp.write(uploaded_file.getvalue())
    tmp.flush()
    tmp.close()
    return Path(tmp.name)


def run_analysis(
    sequence_path: Path,
    selected_enzymes: list[str],
    single_stranded: bool,
) -> tuple[dict, list, plt.Figure, plt.Figure]:
    """
    Core analysis pipeline — mirrors the logic in src/main.py but returns
    matplotlib Figure objects instead of saving PNGs to disk.

    Parameters
    ----------
    sequence_path : Path
        Path to the plasmid sequence file (FASTA or GenBank).
    selected_enzymes : list[str]
        Enzyme names to search for (must be present in the database).
    single_stranded : bool
        If True, produce a single-stranded linear map;
        otherwise produce a double-stranded linear map.

    Returns
    -------
    results : dict
        Raw motif-search results keyed by enzyme name.
        Each value is [motif, cut_notation, count, [positions]].
    plasmid : list[str]
        [header, sequence] parsed from the input file.
    fig_circular : matplotlib.figure.Figure
        Rendered circular plasmid map.
    fig_linear : matplotlib.figure.Figure
        Rendered linear (single- or double-stranded) plasmid map.

    Raises
    ------
    ValueError
        If no enzymes are found in the database, or the sequence is empty.
    """
    # --- 1. Load and parse the sequence ----------------------------------
    seq_obj = Sequence(str(sequence_path))
    seq_obj.load_sequence()
    plasmid = seq_obj.sequence  # [header, nucleotide_string]

    if not plasmid or not plasmid[1]:
        raise ValueError("Could not parse a nucleotide sequence from the uploaded file.")

    # --- 2. Filter enzymes from the database -----------------------------
    re_list = Enzymes(str(DB_FILE))
    # Pass interface=False to skip the TUI; enzymes come from the widget instead.
    re_list.interface(selected_enzymes, interface=False)
    re_list.filter_enzymes()
    enzymes = re_list.filtered

    if not enzymes:
        raise ValueError(
            "None of the selected enzymes were found in the database. "
            "Check the enzyme names against src/database/enzymes.csv."
        )

    # --- 3. Locate motifs in the plasmid ---------------------------------
    motif_loc = Motifs(enzymes)
    motif_loc.array_set(plasmid[1])
    results = motif_loc.get_motif_results()

    # --- 4. Render maps (output_path=None → do not save to disk) ---------
    plasmid_map = PlasmidMap(
        results          = results,
        plasmid_sequence = plasmid[1],
        title            = plasmid[0],
    )

    # Circular map — always generated
    fig_circular = plasmid_map.annotate_circular(output_path=None)

    # Linear map — single- or double-stranded depending on checkbox
    if single_stranded:
        fig_linear = plasmid_map.annotate_linear(output_path=None)
    else:
        fig_linear = plasmid_map.annotate_double_stranded(output_path=None)

    return results, plasmid, fig_circular, fig_linear


def fig_to_png_bytes(fig: plt.Figure, dpi: int = 150) -> bytes:
    """
    Render a matplotlib Figure to a PNG byte string suitable for
    st.download_button() or st.image().

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        A fully rendered figure.
    dpi : int
        Resolution of the exported PNG.

    Returns
    -------
    bytes
        Raw PNG bytes.
    """
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
    buf.seek(0)
    return buf.read()


def results_to_csv_bytes(results: dict) -> bytes:
    """
    Serialise the motif-search results dictionary to a CSV byte string for
    use with st.download_button().

    Parameters
    ----------
    results : dict
        Output of Motifs.get_motif_results():
        { enzyme: [motif, cut_site, count, [positions]] }

    Returns
    -------
    bytes
        UTF-8 encoded CSV bytes.
    """
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow(["enzyme", "motif", "cut_site", "observed_count", "start_positions"])
    for enzyme, (motif, cut_site, count, positions) in results.items():
        writer.writerow([enzyme, motif, cut_site, count, positions])
    return buf.getvalue().encode("utf-8")


def results_to_dataframe(results: dict) -> pd.DataFrame:
    """
    Convert the motif-search results dictionary to a pandas DataFrame for
    display with st.dataframe().

    Parameters
    ----------
    results : dict
        Output of Motifs.get_motif_results().

    Returns
    -------
    pd.DataFrame
        Columns: enzyme, motif, cut_site, observed_count, start_positions.
        Rows are sorted by observed_count descending.
    """
    rows = []
    for enzyme, (motif, cut_site, count, positions) in results.items():
        rows.append({
            "enzyme"          : enzyme,
            "motif"           : motif,
            "cut_site"        : cut_site,
            "observed_count"  : count,
            "start_positions" : ", ".join(str(p) for p in positions) if positions else "—",
        })
    df = pd.DataFrame(rows)
    return df.sort_values("observed_count", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Page config  (must be the first Streamlit call)
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title = "PlasMap — Plasmid Map Generator",
    page_icon  = "🧬",
    layout     = "wide",
)


# ---------------------------------------------------------------------------
# Sidebar — inputs & settings
# ---------------------------------------------------------------------------
with st.sidebar:
    st.title("🧬 PlasMap")
    st.caption("Annotated plasmid map generator")
    st.divider()

    # ---- Sequence input --------------------------------------------------
    st.subheader("1 · Sequence input")

    input_mode = st.radio(
        "Input mode",
        options   = ["Upload file", "Use demo plasmid"],
        horizontal = True,
        help       = "Upload your own FASTA/GenBank file or try one of the bundled demo plasmids.",
    )

    sequence_path: Path | None = None  # resolved below

    if input_mode == "Upload file":
        # Accept FASTA and GenBank extensions
        uploaded = st.file_uploader(
            "Upload plasmid sequence",
            type   = list(VALID_EXTS),
            help   = "Accepted formats: FASTA (.fasta .fa .fas …) or GenBank (.gb .gbk)",
        )
        if uploaded:
            # Write to a temp file so Sequence() can open it by path
            sequence_path = save_uploaded_file(uploaded)

    else:  # demo plasmid
        demo_label = st.selectbox(
            "Choose demo plasmid",
            options = list(DEMO_FILES.keys()),
        )
        sequence_path = DEMO_FILES[demo_label]

    # ---- Enzyme selection ------------------------------------------------
    st.divider()
    st.subheader("2 · Enzyme selection")

    # Load every enzyme name from the database for the full picker
    all_enzymes = load_all_enzyme_names(DB_FILE)

    # Seed the multiselect key with the default enzyme list on first load.
    # Buttons mutate this key directly before the widget renders, so the
    # widget reflects the new value on the same rerun without needing `default=`.
    if "multiselect_enzymes" not in st.session_state:
        st.session_state["multiselect_enzymes"] = DEFAULT_ENZYMES

    col_a, col_b = st.columns(2)
    with col_a:
        if st.button("Default 22", width="stretch"):
            st.session_state["multiselect_enzymes"] = DEFAULT_ENZYMES
    with col_b:
        if st.button("Clear all", width="stretch"):
            st.session_state["multiselect_enzymes"] = []

    selected_enzymes: list[str] = st.multiselect(
        "Restriction enzymes",
        options  = all_enzymes,
        key      = "multiselect_enzymes",
        help     = "Select one or more enzymes. The default panel is the top-20 common Type II enzymes.",
    )

    # ---- Map options -----------------------------------------------------
    st.divider()
    st.subheader("3 · Map options")

    single_stranded: bool = st.checkbox(
        "Single-stranded linear map",
        value = False,
        help  = (
            "When checked, the linear map shows only the forward (sense) strand. "
            "When unchecked, both strands with base-pairing ticks are shown."
        ),
    )

    # ---- Run button ------------------------------------------------------
    st.divider()
    run_button = st.button(
        "▶  Generate maps",
        type              = "primary",
        width="stretch",
        disabled          = (sequence_path is None or len(selected_enzymes) == 0),
    )

    if sequence_path is None:
        st.info("Upload a sequence file or choose a demo plasmid to continue.")
    if len(selected_enzymes) == 0:
        st.warning("Select at least one enzyme.")


# ---------------------------------------------------------------------------
# Main panel — results
# ---------------------------------------------------------------------------
st.title("PlasMap — Plasmid Map Generator")
st.markdown(
    "Upload a **FASTA** or **GenBank** plasmid file, choose your restriction enzymes, "
    "then click **Generate maps** to visualise cut sites."
)

if run_button and sequence_path is not None and selected_enzymes:
    # Run the full analysis pipeline and cache results in session_state so
    # they survive Streamlit reruns (e.g. when the user clicks a download button).
    with st.spinner("Analysing sequence and rendering maps…"):
        try:
            results, plasmid, fig_circular, fig_linear = run_analysis(
                sequence_path    = sequence_path,
                selected_enzymes = selected_enzymes,
                single_stranded  = single_stranded,
            )
            # Store in session state for download buttons to reference
            st.session_state["results"]       = results
            st.session_state["plasmid"]       = plasmid
            st.session_state["fig_circular"]  = fig_circular
            st.session_state["fig_linear"]    = fig_linear
            st.session_state["single_strand"] = single_stranded

        except (ValueError, IOError) as err:
            st.error(f"Analysis failed: {err}")
            st.stop()

# Render the results whenever they are present in session state
# (persists across reruns so download buttons don't clear the maps)
if "results" in st.session_state:
    results      = st.session_state["results"]
    plasmid      = st.session_state["plasmid"]
    fig_circular = st.session_state["fig_circular"]
    fig_linear   = st.session_state["fig_linear"]

    # ---- Summary metrics ------------------------------------------------
    st.divider()
    plasmid_name  = plasmid[0].split(",")[0]   # trim long GenBank definition lines
    plasmid_len   = len(plasmid[1])
    cutters       = sum(1 for v in results.values() if v[2] > 0)
    non_cutters   = sum(1 for v in results.values() if v[2] == 0)

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Plasmid",          plasmid_name)
    col2.metric("Length (bp)",      f"{plasmid_len:,}")
    col3.metric("Cutting enzymes",  cutters)
    col4.metric("Non-cutters",      non_cutters)

    # ---- Maps -----------------------------------------------------------
    st.divider()
    tab_circ, tab_lin, tab_table = st.tabs(["🔵 Circular map", "📏 Linear map", "📋 Results table"])

    with tab_circ:
        st.subheader("Circular plasmid map")
        # Convert figure → PNG bytes for display (avoids pyplot deprecation warning)
        circ_bytes_display  = fig_to_png_bytes(fig_circular, dpi=100)  # screen
        circ_bytes_download = fig_to_png_bytes(fig_circular, dpi=150)  # print
        st.image(circ_bytes_display, width="stretch")

        st.download_button(
            label    = "⬇  Download circular map (PNG)",
            data     = circ_bytes_download,
            file_name = "circular_map.png",
            mime     = "image/png",
        )

    with tab_lin:
        map_type = "single-stranded" if st.session_state.get("single_strand") else "double-stranded"
        st.subheader(f"Linear plasmid map ({map_type})")
        lin_bytes_display  = fig_to_png_bytes(fig_linear, dpi=100)  # screen
        lin_bytes_download = fig_to_png_bytes(fig_linear, dpi=150)  # print
        st.image(lin_bytes_display, width="stretch")

        st.download_button(
            label     = "⬇  Download linear map (PNG)",
            data      = lin_bytes_download,
            file_name = f"{map_type.replace('-','_')}_linear_map.png",
            mime      = "image/png",
        )

    with tab_table:
        st.subheader("Restriction enzyme results")
        df = results_to_dataframe(results)

        # Highlight rows with at least one cut site.
        # Use a muted teal with explicit dark text so the label stays readable
        # in both light and dark Streamlit themes.  The light green (#d4edda)
        # used previously blended into dark-mode backgrounds and hid the text.
        def _highlight_cutters(row):
            """
            Apply a theme-safe highlight to enzymes that cut the plasmid.

            A mid-saturation teal background (#1a6b55) with white text is
            visible on both light and dark Streamlit themes.  Non-cutters
            receive no styling, so they inherit the theme default.
            """
            if row["observed_count"] > 0:
                style = "background-color: #1a6b55; color: #ffffff;"
            else:
                style = ""
            return [style] * len(row)

        st.dataframe(
            df.style.apply(_highlight_cutters, axis=1),
            width="stretch",
            height              = 450,
        )

        # CSV download
        csv_bytes = results_to_csv_bytes(results)
        st.download_button(
            label     = "⬇  Download results (CSV)",
            data      = csv_bytes,
            file_name = "plasmid_results.csv",
            mime      = "text/csv",
        )

else:
    # ---- Placeholder shown before first run ----------------------------
    st.info(
        "Configure the sequence and enzyme options in the sidebar, "
        "then click **▶ Generate maps** to run the analysis."
    )

    # Show a brief feature summary so users know what to expect
    col_f1, col_f2, col_f3 = st.columns(3)
    with col_f1:
        st.markdown("**🔵 Circular map**")
        st.caption("Radial view of the plasmid with labelled cut sites and position ruler.")
    with col_f2:
        st.markdown("**📏 Linear map**")
        st.caption("Base-level sequence view with highlighted recognition motifs and cut annotations.")
    with col_f3:
        st.markdown("**📋 Results table**")
        st.caption("Enzyme name, motif, cut notation, count, and all start positions as a downloadable CSV.")