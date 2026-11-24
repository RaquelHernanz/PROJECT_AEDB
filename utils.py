
### CONSTANTS
GAP = "-"

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

def print_msa(msa: list[str]) -> None:
    """
    Print an MSA, one sequence per line.

    Parameters
    msa: list[str]. Alignment to print, one string per sequence.
    """
    for sequence in msa:
        print(sequence)

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
