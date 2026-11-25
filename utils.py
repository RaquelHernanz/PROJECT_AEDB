#Install biopython for the project and script (pip install biopython)
import Bio
from Bio import SeqIO

#-------------------------------------------- CONSTANTS -------------------------------------------

GAP = "-"

#-------------------------------------------- SUBMATRIX FUNCTIONS -------------------------------------------

def make_sub_matrix(match: int, mismatch: int, alphabet: str) -> dict:
    sub_matrix = {}
    for char1 in alphabet:
        sub_matrix[char1] = {}
        for char2 in alphabet:
            if char1 == char2:
                sub_matrix[char1][char2] = match
            else:
                sub_matrix[char1][char2] = mismatch
    return sub_matrix

def read_sub_matrix(path: str) -> dict:
    sub_matrix = {}
    with open(path, "r") as f:
        lines = f.readlines()
    alphabet = lines[0].split()
    for line in lines[1:]:
        fields = line.split()
        letter = fields[0]
        scores = fields[1:]
        sub_matrix[letter] = {}
        for i in range(len(alphabet)):
            sub_matrix[letter][alphabet[i]] = int(scores[i])
    return sub_matrix


#-------------------------------------------- FASTA FILE -------------------------------------------

def read_fasta(filename): 
    """Reading the FASTA file. Returns a list with all sequences."""

    # SeqIO.parse lê qualquer extensão (.fasta, .txt, .fa) diretamente
    records = list(SeqIO.parse(filename, "fasta"))

    #Biopython stores sequences as Seq objects. When converting to str
    sequences = [str(record.seq).upper() for record in records]

    # Check if the FASTA file contains sequences.
    if not sequences:
        raise ValueError("No sequences found in the FASTA file.")

    return sequences 

#-------------------------------------------- PRINT MSA -------------------------------------------

def print_msa(msa: list[str]) -> None:
    """
    Print an MSA, one sequence per line.

    Parameters
    msa: list[str]. Alignment to print, one string per sequence.
    """
    for sequence in msa:
        print(sequence)

#-------------------------------------------- GAP VERIFICATION IN MSA COLUMNS -------------------------------------------

def is_all_gaps(msa: list[str], column: int) -> bool:
    """
    Check whether a column of an MSA is composed solely of gaps.

    Parameters
    msa: list[str]. Rectangular alignment.
    column: int. Index of the column to check (0-based).

    Returns
    bool: True if every sequence has a gap at `column`, False otherwise.
    """
    for sequence in msa:
        if sequence[column] != GAP:
            return False
    return True


#-------------------------------------------- CLEAN ALIGNMENT -------------------------------------------

def clean_alignment(msa: list[str]) -> list[str]:
    """
    Remove the columns composed solely of gaps from an MSA.

    Parameters
    msa: list[str]. Rectangular alignment, possibly with gap-only columns.

    Returns
    list[str]. New alignment without gap-only columns.

    Example
    We have a column made of -
    ATC-GG---TT
    A-C-CGCA-TC
    A-C---CGATC

    ATCGG---TT
    A-CCGCA-TC
    A-C--CGATC

    """
    if not msa:
        return []

    kept_columns: list[int] = []
    for column in range(len(msa[0])):
        #We check which index columns are not full of gaps with is_all_gaps()
        if not is_all_gaps(msa, column):
            kept_columns.append(column)

    cleaned: list[str] = []
    #Ad we add the nucleotides to the cleaned sequence, without the gap columns
    for sequence in msa:
        new_sequence = ""
        for column in kept_columns:
            new_sequence += sequence[column]
        cleaned.append(new_sequence)
    return cleaned


#-------------------------------------------- MSA VALIDATION -------------------------------------------

def is_valid_msa(msa: list[str], original_sequences: list[str]) -> bool:
    """
    Check that an MSA is valid for a given set of sequences.

    An MSA is valid if:
    - it has one row per original sequence;
    - all rows have the same length;
    - there are no gap-only columns;
    - removing the gaps from each row gives back the original sequence
      (same order).

    Use it after initialization, after each mutation and, from the
    tests, on the crossover offspring.

    Parameters
    msa : list[str]. Alignment to check.
    original_sequences : list[str]. Unaligned sequences the MSA should contain.

    Returns
    bool
        True if the MSA is valid, False otherwise.
    """
    if len(msa) != len(original_sequences):
        return False
    if not msa:
        return True

    # All rows with the same length
    length = len(msa[0])
    for row in msa:
        if len(row) != length:
            return False

    # No gap-only columns
    for column in range(length):
        only_gaps = True
        for row in msa:
            if row[column] != GAP:
                only_gaps = False
                break
        if only_gaps:
            return False

    # Same residues, in the same order
    for i in range(len(msa)):
        if msa[i].replace(GAP, "") != original_sequences[i]:
            return False

    return True

if __name__ == "__main__":
    FILE = "cytochromes.fa"   # <- cambia por vuestro archivo

    seqs = read_fasta(FILE)
    print(len(seqs))                                   # nº of expected seqs
    print([len(s) for s in seqs])                      # expected lenght
    print(all(s == s.upper() for s in seqs))           # True
    print(all(GAP not in s for s in seqs))             # True (without gaps)
    print(all(set(s) <= set("ACGT") for s in seqs))    # True if ADN (for proteínas, FALSE)

    """
    13
    [105, 105, 105, 105, 105, 105, 104, 104, 104, 104, 105, 104, 104]
    True
    True
    False
    """