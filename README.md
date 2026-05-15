# PlasMap - Plasmid annotation tool

>Click here for live demo: [PlasMap]()

## Overview:
**PlasMap** is a python tool that takes a plasmid sequence file as input and generates annotated sequence maps with Restriction Enzyme cut sites.

### Features:
- **CSV Results File** - Results table containing enzyme list with motif, cut type, and binding location. 
- **Circular map**  -  Circular plasmid map displaying enzyme cut locations to visualize plasmid topology and identify candidate enzymes for cloning experiments.  
- **Double-stranded linear map** - Sequence-level map of the entire plasmid containing highlighted recognition motifs and cut marks on forward and complimentary strands. This allows visualization of enzyme overlap and overhang.  
- **Single-stranded linear map** - Optional map containing the same sequence-level map, but with only the forward strand visualized. 
- **Terminal User Interface** - Optional mode that displays enzyme list in the terminal for the user to select from. 
- IUPAC ambiguity code support in recognition sequences
- GenBank and FASTA format compatible


## Project File Structure:
```
└── 📁PlasMap
    └── 📁inputs
        └── 📁test
            ├── pUC19.fasta
            ├── ...
    └── 📁src
        └── 📁database
            ├── enzymes.csv
        └── 📁motif_id_lib
            ├── __init__.py
            ├── csv_output.py
            ├── input.py
            ├── motif_locator.py
            ├── output.py
        ├── main.py         <---- CLI entry point
    └── 📁web
        ├── app.py          <---- Web UI entry point
        ├── UI-README.md
    ├── .gitignore
    ├── dependencies.txt
    ├── environment.yml
    ├── LICENSE
    ├── pseudocode.md
    └── README.md
``` 

1. `inputs`: Folder to upload plasmid sequences for analysis.
    - Contains `test` folder with example plasmid files.
2. `database` : Contains a CSV file of restriction enzymes and recognition sequence motifs. 
3. `src`: Contains the logic for the program.
    - Loads the user input sequence and restriction enzymes database. 
    - Parses the plasmid sequence with each enzyme motif to identify cut locations.
    - Produces output containing an annotated map of the sequence with the enzyme cut sites in the `results` folder.
4. `web`: Contains web UI to demo the program:
    - Uses Streamlit package to run and host 
    - Graphical UI uses same modules as main program
    - Converts user CLI arguments to widgets for input
    - Generates the same output files for download. 

---

## Installation

>Ensure conda environment manager is installed first. 

Create project environment:
```bash
conda env create -f environment.yml
```
Conda will automatically create an environment named `plasmap` with all the specified packages and versions.

---

## Usage:

### 1. Activate the environment:
```bash
conda activate plasmap
```
### 2. Run command:

#### Example 1: Default
- Double stranded map
- Default pUC19 plasmid - FASTA format
- Default enzyme list (all 20). 
```bash
python src/main.py
```

#### Example 2: User plasmid
- Double stranded map
- Default enzyme list (all 20)
- User plasmid pCMV-GLuc (mammalian vector) - GenBank format
```bash
python src/main.py -s inputs/test/pCMV-GLuc.gb
```

#### Example 3: User enzyme list
- Double stranded map
- User provided list of enzymes
- Default pUC19 plasmid
```bash
python src/main.py -e EcoRI HindIII BamHI
```
>NOTE: Refer to `database/enzymes.csv` for full list of Type II restriction enzymes.

#### Example 4: Single stranded
- Single stranded map
- Single stranded M13mp18 bacteriophage vector
- Default enzyme list (all 20)
- 
```bash
python src/main.py -s inputs/test/M13mp18.gb -ss
```

#### BONUS Example: TUI
- Terminal User interface for enzyme selection

- Double stranded map
- Default pUC19 plasmid
```bash
python src/main.py -i
```
>Use keyboard arrows, and spacebar to select enzymes. Press enter to submit.
---

### Command-Line Arguments:
| Argument                | Description                                  | Default         |
| ----------------------- | -------------------------------------------- | --------------- |
| `-s`. `--sequence_filepath`| Plasmid sequence filepath. Accepts FASTA and GenBank formats. | inputs/test/pUC19.fasta         |
| `-c`. `--csv_output`| Filepath to export CSV results. | results/plasmid_results.csv |
| `-ss`, `--single_stranded`        | Output linear map as single-stranded.  | False|
| `-e`, `--enzymes`       | User list of restriction enzyme names     | (list of top 20)|
| `-i`, `--interface`           | Optional User interface that displays enzyme list for selection. Use arrows and space bar to make selection.               | False|

## Output:
All saved to `results` folder:
1. CSV file containing enzyme cut results: 
    - enzyme
    - motif
    - cut_site
    - observed_count
    - start_positions
2. Circular annotated plasmid map `.png`
3. Linear annotated plasmid map (Double stranded or Single stranded) `.png`

---
## License: 
**GNU General Public License Version 3**

The GNU GPL is a license that ensures code is open-source. GNU GPL allows others to utilize, modify, or distribute code. If other users modify the code, then these users are expected to share their changes to the code under a GNU GPL to maintain the open-source integrity of the code.

--- 
## References:

### Python Standard Library
**`re` - Regular Expressions**
Python Software Foundation. (2024). *re - Regular expression operations*. Python 3 Documentation.
https://docs.python.org/3/library/re.html
Used to expand IUPAC ambiguity codes into regex character classes and locate recognition motif positions within the plasmid sequence via `re.finditer()`.
 
**`math`**
Python Software Foundation. (2024). *math - Mathematical functions*. Python 3 Documentation.
https://docs.python.org/3/library/math.html
Used for `math.ceil()` (line-count calculations), `math.floor()` and `math.log10()` (ruler tick interval rounding), and `math.pi` (circular angle arithmetic).
 
**`random`**
Python Software Foundation. (2024). *random - Generate pseudo-random numbers*. Python 3 Documentation.
https://docs.python.org/3/library/random.html
Used to shuffle the enzyme colour palette with a fixed seed for reproducible random colour assignment.

### Third-Party Libraries
**NumPy**

NumPy Developers. (n.d.). NumPy reference: Routines. https://numpy.org/doc/stable/reference/routines.html

NumPy Developers. (n.d.). numpy.lib.stride_tricks.sliding_window_view. https://numpy.org/doc/stable/reference/generated/numpy.lib.stride_tricks.sliding_window_view.html

**Matplotlib**

Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9(3), 90–95.
https://doi.org/10.1109/MCSE.2007.55
https://matplotlib.org/


### Bioinformatics Concepts & Standards

**IUPAC Nucleotide Ambiguity Codes**
Nomenclature Committee of the International Union of Biochemistry (NC-IUB). (1985). Nomenclature for incompletely specified bases in nucleic acid sequences. *European Journal of Biochemistry*, 150(1), 1–5.
https://doi.org/10.1111/j.1432-1033.1985.tb08977.x
The `IUPAC` dictionary maps ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) to their corresponding regex character classes for motif matching.
 
**Restriction Enzyme Cut Notation**
Rebase - The Restriction Enzyme Database. Roberts, R. J., Vincze, T., Posfai, J., & Macelis, D. (2023). REBASE: a database for DNA restriction and modification: enzymes, genes and genomes. *Nucleic Acids Research*, 51(D1), D629–D630.
https://doi.org/10.1093/nar/gkac975
https://rebase.neb.com/
The `^` (top-strand cut) and `_` (bottom-strand cut) notation parsed by `_top_cut_offset()` and `_bot_cut_offset()` follows the REBASE convention for describing staggered and blunt restriction enzyme cleavage sites.
 
**DNA Complementarity**
Watson, J. D., & Crick, F. H. C. (1953). Molecular structure of nucleic acids: A structure for deoxyribose nucleic acid. *Nature*, 171, 737–738.
https://doi.org/10.1038/171737a0
The `COMPLEMENT` dictionary (A↔T, G↔C) used to generate the bottom strand in `DoubleStrandedMap` is based on Watson–Crick base-pairing rules.

### Visualization Design

**Circular Plasmid Map Style**
SnapGene Viewer. GSL Biotech LLC. (2024). *SnapGene - Plasmid map visualisation*.
https://www.snapgene.com/
NEB Cutter. Vincze, T., Posfai, J., & Roberts, R. J. (2003). NEBcutter: A program to cleave DNA with restriction enzymes. *Nucleic Acids Research*, 31(13), 3688–3691.
https://doi.org/10.1093/nar/gkg526
The radial stacking algorithm in `_compute_label_radii()` - pushing overlapping labels outward in discrete radial steps - was designed to reproduce the label layout style used by these tools.
 
**Colour Palette Design**
Crameri, F., Shephard, G. E., & Heron, P. J. (2020). The misuse of colour in science communication. *Nature Communications*, 11, 5444.
https://doi.org/10.1038/s41467-020-19160-7
Informed the decision to use perceptually distinct, high-contrast colours for enzyme annotation rather than sequential or single-hue palettes.

### AI assistance:

This project was developed with the help of [Claude Sonnet 4.6](https://claude.ai) (`claude-sonnet-4-6`) by [Anthropic](https://anthropic.com).

Claude assisted with:
- Code architecture and implementation
- Docstring and documentation writing
- Debugging and code review

All generated code was reviewed and tested by the author.

## Author:

Michael Villarreal - mvillar6@charlotte.edu | mikev8492@gmail.com