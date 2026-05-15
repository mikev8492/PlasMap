# PlasMap – Pseudocode Overview

## Project Summary

PlasMap is a plasmid visualization and restriction analysis workflow that:

1. Loads a plasmid DNA sequence
2. Loads a restriction enzyme list
3. Scans the sequence for enzyme recognition sites
4. Detects forward and reverse-complement matches
5. Annotates restriction cut locations
6. Generates circular and linear plasmid maps
7. Produces visual and/or exportable outputs


## High-Level Workflow

```text
START
│
├── Read plasmid sequence file
├── Validate DNA sequence
├── Read restriction enzyme database/list
├── For each enzyme:
│     ├── Get recognition sequence
│     ├── Search forward strand
│     ├── Search reverse complement
│     ├── Store cut site positions
│     └── Store annotation metadata
│
├── Generate plasmid feature annotations
├── Create circular plasmid map
├── Create linear plasmid map
├── Render labels and enzyme positions
├── Export figures/output files
│
└── END
```

# Pseudocode:

### Main Program 

```python
FUNCTION main():

    # -----------------------------------
    # Load Inputs
    # -----------------------------------
    plasmid_sequence = load_fasta(input_fasta)

    enzyme_database = load_enzyme_list(enzyme_file)


    # -----------------------------------
    # Validate Sequence
    # -----------------------------------
    IF not is_valid_dna(plasmid_sequence):
        print("Invalid DNA sequence")
        EXIT PROGRAM


    # -----------------------------------
    # Analyze Restriction Sites
    # -----------------------------------
    restriction_results = []

    FOR each enzyme IN enzyme_database:

        recognition_site = enzyme.sequence

        forward_matches = find_matches(
            plasmid_sequence,
            recognition_site
        )

        reverse_site = reverse_complement(recognition_site)

        reverse_matches = find_matches(
            plasmid_sequence,
            reverse_site
        )

        all_matches = combine_matches(
            forward_matches,
            reverse_matches
        )

        enzyme_record = {
            "name": enzyme.name,
            "site": recognition_site,
            "positions": all_matches,
            "cut_pattern": enzyme.cut_pattern
        }

        append restriction_results with enzyme_record


    # -----------------------------------
    # Generate Feature Annotations
    # -----------------------------------
    annotations = create_annotations(
        restriction_results
    )


    # -----------------------------------
    # Generate Circular Map
    # -----------------------------------
    circular_map = generate_circular_map(
        sequence=plasmid_sequence,
        annotations=annotations
    )


    # -----------------------------------
    # Generate Linear Map
    # -----------------------------------
    linear_map = generate_linear_map(
        sequence=plasmid_sequence,
        annotations=annotations
    )


    # -----------------------------------
    # Export Results
    # -----------------------------------
    save_figure(circular_map, "circular_map.svg")

    save_figure(linear_map, "linear_map.svg")

    export_site_table(restriction_results)


    print("Plasmid analysis complete")

END FUNCTION
```

---

### FASTA Loading Logic

```python
FUNCTION load_fasta(file_path):

    sequence = ""

    OPEN file_path AS fasta_file

    FOR each line IN fasta_file:

        IF line starts with ">":
            CONTINUE

        sequence += strip_whitespace(line)

    RETURN uppercase(sequence)
```

---

### DNA Validation Logic

```python
FUNCTION is_valid_dna(sequence):

    valid_bases = ["A", "T", "C", "G", "N"]

    FOR each base IN sequence:

        IF base NOT IN valid_bases:
            RETURN False

    RETURN True
```

---

### Restriction Site Search Logic

```python
FUNCTION find_matches(sequence, motif):

    match_positions = []

    motif_length = length(motif)


    FOR i FROM 0 TO length(sequence) - motif_length:

        current_window = sequence[i : i + motif_length]

        IF current_window == motif:
            append match_positions with i


    RETURN match_positions
```

---

### Reverse Complement Logic

```python
FUNCTION reverse_complement(sequence):

    complement_table = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    reversed_sequence = reverse(sequence)

    reverse_complement_sequence = ""


    FOR each base IN reversed_sequence:

        reverse_complement_sequence += complement_table[base]


    RETURN reverse_complement_sequence
```

---

### Circular Map Rendering Logic

```python
FUNCTION generate_circular_map(sequence, annotations):

    plasmid_length = length(sequence)

    initialize_canvas()

    draw_outer_circle()


    FOR each annotation IN annotations:

        angle = convert_position_to_angle(
            annotation.position,
            plasmid_length
        )

        draw_tick_mark(angle)

        draw_label(
            text=annotation.name,
            angle=angle
        )


    RETURN completed_canvas
```

---

### Linear Map Rendering Logic

```python
FUNCTION generate_linear_map(sequence, annotations):

    plasmid_length = length(sequence)

    initialize_canvas()

    draw_horizontal_backbone()


    FOR each annotation IN annotations:

        x_coordinate = scale_position(
            annotation.position,
            plasmid_length
        )

        draw_vertical_marker(x_coordinate)

        draw_label(
            text=annotation.name,
            x=x_coordinate
        )


    RETURN completed_canvas
```

---

## Data Structures

### Enzyme Record

```python
enzyme_record = {
    "name": "EcoRI",
    "site": "GAATTC",
    "positions": [120, 842, 2101],
    "cut_pattern": "G^AATTC"
}
```
### Search Results
```
{enzyme name: [motif, cuts, count of cuts, cut positions]}
```
```python
results = {
    'Eco24I': ['GRGCYC', 'G_RGCY^C', 1, [277]],
    'Eco91I': ['GGTNACC', 'G^GTNAC_C', 0, []], 
    'SpeI': ['ACTAGT', 'A^CTAG_T', 0, []], 
    'EcoRI': ['GAATTC', 'G^AATT_C', 1, [283]], 
    'PstI': ['CTGCAG', 'C_TGCA^G', 1, [244]], 
    'HindIII': ['AAGCTT', 'A^AGCT_T', 1, [232]], 
    'NotI': ['GCGGCCGC', 'GC^GGCC_GC', 0, []], 
    'ClaI': ['ATCGAT', 'AT^CG_AT', 0, []], 
    'ApaI': ['GGGCCC', 'G_GGCC^C', 0, []], 
    'BglII': ['AGATCT', 'A^GATC_T', 0, []], 
    'HaeIII': ['GGCC', 'GG^_CC', 6, [35, 292, 394, 692, 1279, 1546]], 
    'XbaI': ['TCTAGA', 'T^CTAG_A', 1, [256]], 
    'EcoRV': ['GATATC', 'GAT^_ATC', 0, []], 
    'SacI': ['GAGCTC', 'G_AGCT^C', 1, [277]], 
    'SalI': ['GTCGAC', 'G^TCGA_C', 1, [250]], 
    'KpnI': ['GGTACC', 'G_GTAC^C', 1, [271]], 
    'SphI': ['GCATGC', 'G_CATG^C', 1, [238]], 
    'SmaI': ['CCCGGG', 'CCC^_GGG', 1, [267]], 
    'NdeI': ['CATATG', 'CA^TA_TG', 1, [496]], 
    'BamHI': ['GGATCC', 'G^GATC_C', 1, [262]], 
    'XhoI': ['CTCGAG', 'C^TCGA_G', 0, []], 
    'MluI': ['ACGCGT', 'A^CGCG_T', 0, []]
}

```
