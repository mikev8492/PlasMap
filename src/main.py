#!/usr/bin/env python3

import sys, argparse, traceback, time
from pathlib import Path
from Bio.Restriction import CommOnly
from motif_id_lib.input import Sequence, Enzymes
from motif_id_lib.motif_locator import Motifs
from motif_id_lib.output import PlasmidMap
from motif_id_lib.csv_output import CreateCSV

"""
This file contains the main functions that will perform the logic for the program. 
"""
# --- Constants -----------------------------
# sys.stdout colors:
HIGHLIGHT = "\033[36m" # teal
RESET = "\033[0m"

DB_FILE = "src/database/enzymes.csv"

DEFAULT_ENZYMES = [
            "EcoRI",
            "HindIII",
            "BamHI",
            "XhoI",
            "NotI",
            "SalI",
            "PstI",
            "KpnI",
            "XbaI",
            "EcoRV",
            "SmaI",
            "NdeI",
            "SacI",
            "SpeI",
            "BglII",
            "ApaI",
            "SphI",
            "MluI",
            "ClaI",
            "HaeIII",
            "Eco91I",
            "Eco24I"
        ]

def create_parser() -> argparse.Namespace:
    '''
    Purpose:
    --------
        Defines CLI arguments for the program
    '''
    motif_parser = argparse.ArgumentParser(description="REcut: Plasmid Sequence Cutting Tool - Generates an annotated plasmid map with restriction enzyme cut sites.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    motif_parser.add_argument(
        "-s", "--sequence_filepath",
        type=str, 
        default="inputs/test/pUC19.fasta", 
        help="Path to the input plasmid sequence file. Accepted formats: FASTA (fasta, .fas, .fa, .fna, .ffn, .faa, .mpfa, .frn) or GenBank (.gb, .gbk)."
    )

    motif_parser.add_argument(
        "-i", "--interface",
        action="store_true",
        help="Enable TUI (terminal user interface) to allow manual selection of Restriction enzymes to search with."
    )

    motif_parser.add_argument(
        "-ss", "--single_stranded",
        action="store_true",
        help="Output linear plasmid map as single-stranded."
    )

    motif_parser.add_argument(
        "-e", "--enzymes",
        nargs="+",
        default=DEFAULT_ENZYMES,
        help="List of restriction enzyme names to map to plasmid sequence. Refer to src/database/enzymes.csv for exact enzyme name format."
    )

    motif_parser.add_argument(
        "-c", "--csv_output",
        type=str,
        default="results/plasmid_results.csv",
        help="File path to CSV output file."
    )

    return motif_parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    '''
    Purpose:
    --------
        Validate the command-line arguments

    Args:
    -----
        args (argparse.Namespace): Parsed command-line arguments

    Returns:
    -------
        bool: True if arguments are valid, False otherwise

    Raises:
    -------
        ValueError: If any argument is invalid
    '''
    # Check filetype
    genbank_extensions = ('gb', 'gbk')
    fasta_extensions = ('fasta', 'fas', 'fa', 'fna', 'ffn', 'faa', 'mpfa', 'frn')

    input_extension = args.sequence_filepath.split(".")[-1]
    
    if input_extension not in genbank_extensions and input_extension not in fasta_extensions:
        raise ValueError("Invalid file type (-s, --sequence_filepath). Please provide a FASTA or GenBank file.")

    # Check plasmid filepath exists
    filepath = Path(args.sequence_filepath)
    if not filepath.is_file():
        raise IOError("Filepath (-s, --sequence_filepath) does not exist")
    
    csvpath = Path(args.csv_output)
    if not csvpath.parent.exists():
        raise IOError(f"Directory {args.csv_output} does not exist.")


def create_db():
    """
    Purpose:
    --------
        Create restriction enzyme database file (CSV). Uses "CommOnly" (common type II enzymes) list of REBASE data (Bio.Restriction).
    """
    # Check if database exists:
    filepath = Path(DB_FILE)
    if not filepath.is_file():
        sys.stdout.write(f"Creating database {filepath}\n")

        with open(filepath, "w") as file:
            file.write("enzyme,motif,cutInfo\n")
            for enzyme in CommOnly:
                file.write(f"{enzyme},{enzyme.site},{enzyme.elucidate()}\n")

def create_results():
    """
    Purpose:
    --------
        Create results folder for user.
    """
    folder_path = Path("results")
    folder_path.mkdir(parents=True, exist_ok=True)

def main():

    try:
        program_start = time.perf_counter()
        create_db()
        create_results()
        args = create_parser()
        validate_arguments(args)

        # =======================
        # INPUT.py
        # =======================
        # load and parse sequence 
        user_file = Sequence(args.sequence_filepath)
        user_file.load_sequence()

        # plasmid: list [header, seq]
        plasmid = user_file.sequence

        # load and filter enzymes
        re_list = Enzymes(DB_FILE)
        re_list.interface(args.enzymes, args.interface) # Interface optional
        re_list.filter_enzymes()

        # Enzymes dict: 
        #   enzyme: [motif, cutInfo]
        enzymes = re_list.filtered

        # --- OUTPUT: Plasmid sequence, Enzyme list -------------
        sys.stdout.write(f"\nPlasmid:\n\t{plasmid[0]}\n")
        sys.stdout.write("\nSelected Enzyme list:\n")
        for enzyme, info in enzymes.items():
            sys.stdout.write(f"\t{enzyme}: {info}\n")


        # =======================
        # MOTIF_LOCATOR.py
        # =======================
        # Find the motif locations in the plasmid sequence. Return counts, and locations for each enzyme.
        # load enzymes and plasmid sequence into motif locator class.
        motif_locations = Motifs(enzymes)
        motif_locations.array_set(plasmid[1])

        # Search for motif locations in plasmid sequence. Store results in a dictionary.
        results = motif_locations.get_motif_results()

        # =======================
        # CSV_OUTPUT.py
        # =======================
        # Generate a csv file from the motif searching results dictionary.
        csv_output_file = CreateCSV(results, args.csv_output)
        csv_output_file.create_csv_output()
        sys.stdout.write(f"\nRESULTS:\n{'-'*8}")
        sys.stdout.write(f"\n\t1. CSV output written to {HIGHLIGHT}{args.csv_output}{RESET}\n")

        # =======================
        # OUTPUT.py
        # =======================
        plasmid_map = PlasmidMap(results = results, plasmid_sequence = plasmid[1], title = plasmid[0])
        plasmid_map.annotate_circular()
        if args.single_stranded:
            plasmid_map.annotate_linear()
        else:
            plasmid_map.annotate_double_stranded()
        sys.stdout.write("\n-> DONE!\n\n")
        elapsed = time.perf_counter() - program_start
        sys.stdout.write(f"[Execution time: {elapsed:.4f} sec.]\n\n")
        
    # --- Error handling --------------
    except ValueError as err:
        sys.stderr.write(f"Error: {err}\n")
        traceback.print_exc()
        return 1
    except IOError as err:
        sys.stderr.write(f"Error: {err}\n")
        traceback.print_exc()
        return 1
    
if __name__=="__main__":
    sys.exit(main())
