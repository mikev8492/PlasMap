# PlasMap — Streamlit Web Interface

> A browser-based UI for the plasmid annotation tool [PlasMap](https://github.com/mikev8492/PlasMap).


---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Launching the App](#launching-the-app)
- [UI Walkthrough](#ui-walkthrough)
  - [Sidebar](#sidebar)
    - [1 · Sequence Input](#1--sequence-input)
    - [2 · Enzyme Selection](#2--enzyme-selection)
    - [3 · Map Options](#3--map-options)
    - [Generate Maps Button](#generate-maps-button)
  - [Main Panel](#main-panel)
    - [Summary Metrics](#summary-metrics)
    - [Circular Map Tab](#-circular-map-tab)
    - [Linear Map Tab](#-linear-map-tab)
    - [Results Table Tab](#-results-table-tab)
- [CLI vs. Streamlit — Feature Comparison](#cli-vs-streamlit--feature-comparison)
- [File Structure](#file-structure)
- [Notes for Developers](#notes-for-developers)

---

## Overview

`app.py` wraps the PlasMap analysis pipeline with a [Streamlit](https://streamlit.io) web interface. It replaces the `argparse` CLI and the optional terminal UI (`-i` flag) with interactive widgets that run in any browser.

All core logic in `src/motif_id_lib/` is **unchanged**. The Streamlit layer only handles input collection, calling `run_analysis()`, and rendering results — it does not write any files to disk during normal use.

---

## Requirements

All existing PlasMap dependencies are still required, with these additional depenencies:

```yaml
# environment.yml 
- streamlit>=1.35
- pandas
```
---

## Installation

```bash
# 1. Clone the repository 
git clone https://github.com/mikev8492/PlasMap_demo.git
cd PlasMap_demo

# 2. Create and activate the conda environment
conda env create -f environment.yml
conda activate plasmap_demo

```

---

## Launching the App

Run this from the **root** (the same directory that contains `app.py`):

```bash
streamlit run app.py
```

Streamlit will print a local URL (typically `http://localhost:8501`) and open the app in your default browser automatically.

---

## UI Walkthrough

The interface is split into two areas: a **sidebar** on the left for inputs and settings, and a **main panel** on the right for results.

---

### Sidebar

The sidebar walks through three numbered sections. The **▶ Generate maps** button at the bottom activates once all required inputs are set.

---

#### 1 · Sequence Input

Controls how the plasmid sequence is loaded.

**Input mode** (radio buttons)

| Option | Behaviour |
|---|---|
| **Upload file** | Shows a file-upload widget. Drag-and-drop or click to browse. |
| **Use demo plasmid** | Shows a dropdown of five plasmids bundled with the repository. |

**Accepted file formats** (upload mode)

| Format | Extensions |
|---|---|
| FASTA | `.fasta` `.fa` `.fas` `.fna` `.ffn` `.faa` `.mpfa` `.frn` |
| GenBank | `.gb` `.gbk` |

**Bundled demo plasmids** (demo mode)

| Label | Description |
|---|---|
| pUC19 (FASTA) | Standard cloning vector, 2,686 bp, FASTA format |
| pUC19 (GenBank) | Same vector in annotated GenBank format |
| pCMV-GLuc (GenBank) | Mammalian expression vector |
| M13mp18 (GenBank) | Single-stranded bacteriophage vector |
| pBeloBAC11 (GenBank) | Bacterial artificial chromosome |

---

#### 2 · Enzyme Selection

Controls which restriction enzymes are searched against the plasmid sequence.

**Restriction enzymes** (multiselect)

A searchable dropdown containing every enzyme in `src/database/enzymes.csv`. Type to filter by name. Selected enzymes appear as removable tags above the dropdown.

The app loads pre-selected with a default panel of **22 common Type II restriction enzymes**:

```
EcoRI   HindIII  BamHI   XhoI    NotI    SalI    PstI    KpnI
XbaI    EcoRV    SmaI    NdeI    SacI    SpeI    BglII   ApaI
SphI    MluI     ClaI    HaeIII  Eco91I  Eco24I
```

**Quick-select buttons**

| Button | Action |
|---|---|
| **Default 22** | Resets the selection to the 22-enzyme default panel above. |
| **Clear all** | Removes all selected enzymes. |

> **Note:** Enzyme names are case-sensitive and must match the names in `src/database/enzymes.csv` exactly (e.g. `EcoRI`, not `ecori` or `ECORI`).

---

#### 3 · Map Options

**Single-stranded linear map** (checkbox)

| State | Linear map output |
|---|---|
| ☐ Unchecked (default) | **Double-stranded** — both the forward (5'→3') and reverse complementary (3'→5') strands, connected by base-pairing tick marks, with cut annotations on each strand. |
| ☑ Checked | **Single-stranded** — forward strand only, with top-strand cut site annotations. |

This option does not affect the circular map, which is always generated.

---

#### Generate Maps Button

```
▶  Generate maps
```

- Appears **disabled** (greyed out) until both a sequence source and at least one enzyme are selected.
- Displays a spinner (`Analysing sequence and rendering maps…`) while the pipeline runs.
- On success, the results panel updates immediately.
- On failure (unreadable file, no enzyme matches, etc.), an error message is shown in the main panel.

---

### Main Panel

Before running any analysis, the main panel shows a brief feature summary. After clicking **▶ Generate maps**, results appear in four areas.

---

#### Summary Metrics

Four metric cards displayed in a row immediately after a successful run:

| Card | Content |
|---|---|
| **Plasmid** | Sequence identifier or definition line from the file header. |
| **Length (bp)** | Total nucleotide length of the loaded sequence, comma-formatted. |
| **Cutting enzymes** | Number of selected enzymes with at least one recognition site found. |
| **Non-cutters** | Number of selected enzymes with zero matches in this sequence. |

---

#### 🔵 Circular Map Tab

A radial plasmid map rendered at 150 DPI, displayed full-width.

**Map elements:**

- **Backbone circle** — represents the full plasmid sequence.
- **Position ruler** — tick marks and bp labels at 500 bp intervals around the backbone.
- **Cut-site ticks** — each tick is coloured by enzyme and points radially outward from the backbone at the cut position.
- **Leader lines** — dashed lines connecting each tick to its label.
- **Labels** — enzyme name and sequence position `(bp)`, stacked radially outward to avoid overlap when multiple cut sites fall close together.
- **Legend** (lower-left) — colour swatch and enzyme name for every cutting enzyme.
- **Centre text** — plasmid name, total length in bp, and unique cutter count.
- **Non-cutter note** (bottom) — comma-separated list of enzymes with zero cuts, if any.

**Download button:** `⬇ Download circular map (PNG)` — saves `circular_map.png` at 150 DPI.

---

#### 📏 Linear Map Tab

A base-level sequence viewer rendered at 150 DPI, displayed full-width. The tab subtitle shows whether the map is single-stranded or double-stranded based on the sidebar checkbox.

**Map elements (both modes):**

- Sequence is wrapped into lines of 80 bases per row.
- Each nucleotide is coloured by base identity:

  | Base | Colour |
  |---|---|
  | A | Green |
  | T | Red |
  | G | Yellow |
  | C | Blue |

- Recognition motif spans are highlighted with the enzyme's assigned palette colour.
- Top-strand cut sites are annotated with a vertical tick and a labelled arrow pointing upward.
- Position numbers are shown at the left margin of each row.

**Additional elements in double-stranded mode:**

- The reverse complement strand is rendered directly below the forward strand.
- Short vertical tick marks between paired bases indicate Watson–Crick hydrogen bonding.
- Bottom-strand cut site annotations point downward below the complement strand.
- 5'/3' polarity labels appear at each end of both strands.

**Download button:** `⬇ Download linear map (PNG)` — saves `double_stranded_linear_map.png` or `single_stranded_linear_map.png` at 150 DPI.

---

#### 📋 Results Table Tab

An interactive sortable table with one row per enzyme, sorted by cut count descending.

**Columns:**

| Column | Description |
|---|---|
| `enzyme` | Restriction enzyme name. |
| `motif` | Recognition sequence (IUPAC notation). |
| `cut_site` | Cut notation using `^` (top strand) and `_` (bottom strand). |
| `observed_count` | Number of times the recognition motif was found in the sequence. |
| `start_positions` | Comma-separated list of 0-based start positions of each match. Shown as `—` for non-cutters. |

Rows for enzymes with at least one cut site are highlighted in green. Rows for non-cutters have a plain white background.

Click any column header to sort. The table is scrollable with a fixed height of 450 px.

**Download button:** `⬇ Download results (CSV)` — saves `plasmid_results.csv`.

**CSV columns:** `enzyme, motif, cut_site, observed_count, start_positions`

---

## CLI vs. Streamlit — Feature Comparison

| Feature | CLI (`main.py`) | Streamlit (`app.py`) |
|---|---|---|
| Sequence input | File path via `-s` flag | File upload widget or demo dropdown |
| Enzyme selection | `-e` flag or default list | Multiselect widget with full DB access |
| Interactive enzyme picker | `-i` TUI (terminal checkboxes) | Multiselect + quick-select buttons |
| Single-stranded map | `-ss` flag | Checkbox in sidebar |
| Circular map | Always generated, saved to `results/` | Always generated, view in browser + download |
| Linear map | Saved to `results/` | View in browser + download |
| CSV output | Saved to path via `-c` flag | Download button in Results table tab |
| Output destination | `results/` directory on disk | In-browser display + on-demand download |

---

## File Structure

```
└── 📁PlasMap_demo
    ├── app.py                  ← Streamlit UI 
    └── 📁src
        └── 📁motif_id_lib      ← Unchanged core library
            ├── input.py
            ├── motif_locator.py
            ├── output.py
            └── csv_output.py
        └── 📁database
            └── enzymes.csv
        └── main.py             ← Original CLI entry point
    └── 📁inputs
        └── 📁test              ← Demo plasmid files used by the demo dropdown
            ├── pUC19.fasta
            ├── pUC19.gb
            ├── pCMV-GLuc.gb
            ├── M13mp18.gb
            └── pBeloBAC11.gb
    ├── environment.yml
    └── README.md               
```

---

## Notes for Developers

**Session state** — Analysis results (`results`, `plasmid`, `fig_circular`, `fig_linear`) are stored in `st.session_state` after a successful run. This ensures the maps and table remain visible when a user clicks a download button, which triggers a Streamlit rerun.

**No disk writes during analysis** — `PlasmidMap.annotate_circular()` and the linear map methods accept `output_path=None`, which skips the `fig.savefig()` call and returns the figure directly. PNGs are only materialised in memory when a download button is clicked.

**Temporary file handling** — Streamlit's `UploadedFile` object is an in-memory buffer. Because `Sequence()` in `input.py` opens files by path, uploaded files are written to a `NamedTemporaryFile` with the correct extension (so format detection in `load_sequence()` works) and the path is passed downstream.

**TUI replacement** — The `-i` terminal interface (built with `inquirer`) is not used by the Streamlit app. `Enzymes.interface()` is called with `interface=False`, which bypasses the `inquirer.prompt()` call entirely and assigns the multiselect choices directly to `self.usr_list`.

--- 

### AI assistance:

This project was developed with the help of [Claude Sonnet 4.6](https://claude.ai) (`claude-sonnet-4-6`) by [Anthropic](https://anthropic.com).

Claude assisted with:
- Code architecture and implementation
- Docstring and documentation writing
- Debugging and code review

All generated code was reviewed and tested by the author.

## Author:

Michael Villarreal - mvillar6@charlotte.edu | mikev8492@gmail.com