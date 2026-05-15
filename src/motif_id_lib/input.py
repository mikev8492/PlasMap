
import inquirer, re, csv
from termcolor import colored

class Sequence:
    """
    Purpose:
    --------
        Handles reading and parsing of biological sequence files in either
        GenBank (.gb) or FASTA (.fa) format, extracting the header metadata
        and raw nucleotide sequence from each.

    Attributes:
    -----------
        filename (str): Path to the sequence file to be parsed.
        sequence (list): Stores the parsed [header, sequence] after loading.
    """

    def __init__(self,filename) -> None:
        self.filename = filename
        self.sequence = []

    def genBank_parse(self) -> None:
        """
        Purpose:
        --------
            Parses a file in GenBank format, extracting the sequence identifier
            and definition from the header fields, and the nucleotide sequence
            from the ORIGIN section. Appends the header and sequence to
            self.sequence.

        Returns:
        --------
            None. Mutates self.sequence by appending the parsed header (str)
            and sequence (str).
        """
        with open(self.filename) as file:
            header = ""
            seq = ""
            for line in file:
                # Combine header lines to match FASTA header
                match1 = re.match(r'^DEFINITION\s(.*)$', line)
                if match1: header = match1.group(1)
                match2 = re.match(r'^VERSION\s*(.*)$', line)
                if match2: header = match2.group(1) + header 

                # Combine Sequence lines
                if "ORIGIN" in line:
                    lines = file.read()
                    seq_matches = re.findall(r'[a-z]+', lines)
                    seq = "".join(seq_matches).upper()
            self.sequence += header, seq

    def fasta_parse(self) -> None:
        """
        Purpose:
        --------
            Parses a file in FASTA format, reading the first line as the
            header (stripping the leading '>') and concatenating all
            subsequent lines as the nucleotide sequence. Appends the header
            and sequence to self.sequence.

        Returns:
        --------
            None. Mutates self.sequence by appending the parsed header (str)
            and sequence (str).
        """
        with open(self.filename) as file:
            header = file.readline().strip()[1:]
            seq = ""
            for line in file:
                seq += line.strip()
            self.sequence += header, seq

    def load_sequence(self) -> None:
        """
        Purpose:
        --------
            Detects the format of the sequence file based on its extension
            and delegates parsing to the appropriate method (genBank_parse
            for .gb files, fasta_parse for .fa files).

        Returns:
        --------
            None. Indirectly mutates self.sequence via the called parse method. 
        """
        seq_file = self.filename
        if ".gb" in seq_file:
            self.genBank_parse()
            
        elif ".fa" in seq_file:
            self.fasta_parse()
        

class Enzymes:
    """
    Purpose:
    --------
        Manages the full lifecycle of enzyme selection — presenting a
        terminal UI for the user to choose enzymes, then reading the
        CSV database and filtering it down to only the selected entries.

    Attributes:
    -----------
        db (str): Path to the CSV enzyme database file.
        usr_list (list): Enzyme names selected by the user.
        filtered (dict): Maps selected enzyme names to their [motif, cutInfo].
    """

    def __init__(self, database) -> None:
        self.db = database
        self.usr_list = []
        self.filtered = {}


    def app_header(self) -> None:
        """
        Purpose:
            Prints a styled application header to the terminal using colored
            text to visually introduce the TUI (terminal user interface).
        """
        print(colored("\nREcut: Plasmid Sequence Cutting Tool\n", 'cyan', on_color='on_dark_grey', attrs=[ 'blink']))

    def interface(self, enzymes, interface) -> None:
        """
        Purpose:
        --------
            Presents an interactive checkbox UI for the user to select
            enzymes from a provided list. If the interface flag is False,
            bypasses the UI and uses the provided list directly, enabling
            programmatic or testing use.

        Arguments:
        ----------
            enzymes (list[str]): 
                - List of enzyme names to display as options.
            interface (bool):    
                - If True, renders the interactive TUI checkbox.
                - If False, assigns the enzymes list directly to self.usr_list without user interaction.

        Returns:
        --------
            None. Mutates self.usr_list with the selected or provided enzymes.
        """""
        if interface:
            self.app_header()
            usr_enzymes = [
                inquirer.Checkbox('enzymes',
                                    message="Choose your enzymes below (use the arrows and space bar to select):",
                                    choices= enzymes
                                    ),
                ]
            answers = inquirer.prompt(usr_enzymes)
            self.usr_list = answers['enzymes']
        else:
            self.usr_list = enzymes

    def filter_enzymes(self) -> None:
        """
        Purpose:
        --------
            Reads the enzyme CSV database and retains only the rows whose
            enzyme name appears in self.usr_list. Stores the matching
            entries in self.filtered as a dictionary keyed by enzyme name.

        Arguments:
        ----------
            None. Relies on self.db (database path) and self.usr_list
            (user-selected enzyme names) set during initialization and
            interface selection.

        Returns:
        --------
            None. Mutates self.filtered, mapping each selected enzyme
            name (str) to a list of [motif (str), cutInfo (str)].

        Example:
        --------
            self.filtered = {
                "EcoRI": ["GAATTC", "G^GTNAC_C"],
                "BamHI": ["GGATCC", "G^GATC_C"]
            }
        """
        with open(self.db, "r") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if row["enzyme"] in self.usr_list:
                    enz = row["enzyme"]
                    motif = row["motif"]
                    cut = row["cutInfo"]
                    self.filtered[enz] = [motif, cut]
