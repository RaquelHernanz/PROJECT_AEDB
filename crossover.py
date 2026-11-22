"""
Crossover operator for multiple sequence alignments (MSA).

This module contains the functions needed to perform the crossover
operator: two parent MSAs are split at a chosen column (the second 
parent at the equivalent residues) and the halves are recombined into two offspring.

Note: a MSA is represented as a list of strings, one per sequence, where
gaps are indicated by '-'. All strings of a valid MSA have the same
length and there are no columns made only of gaps.
"""

### CONSTANTS
GAP = "-"


### AUXILIARY FUNCTIONS

def print_msa(msa: list[str]) -> None:
    """
    Print an MSA, one sequence per line.

    Parameters
    msa: list[str]. Alignment to print, one string per sequence.
    """
    for sequence in msa:
        print(sequence)


def count_residues(sequence: str) -> int:
    """
    Count the residues (non-gap characters) in an aligned sequence.

    Parameters
    sequence: str. A single line of an MSA, possibly containing gaps.

    Returns
    int. Number of characters different from the gap symbol.

    Example
    >>> count_residues("--AC--")
    2
    """
    residue_count = 0
    for symbol in sequence:
        if symbol != GAP:
            residue_count += 1
    return residue_count


def index_at_residue(sequence: str, n_residues: int) -> int:
    """
    Find the split index that leaves exactly `n_residues` on the left.

    The returned index is the position right after the `n_residues`-th
    residue, so any gaps immediately after it are assigned to the
    right part.

    Parameters
    sequence: str. A single line of an MSA, possibly containing gaps.
    n_residues: int. Number of residues the left part must contain.

    Returns
    int. Index such that ``sequence[:index]`` contains exactly
        `n_residues` residues. Returns 0 if `n_residues` is 0
        (empty left part; leading gaps go to the right part).
    
    Example
        for "--AC--GT" and n_residues = 2, the second residue (C)
        is at index 3, so the function returns 4:
        sequence[:4] = "--AC" and sequence[4:] = "--GT".

    Raises ValueError
        If `n_residues` is negative or larger than the number of
        residues in `sequence`.
    """
    if n_residues < 0:
        raise ValueError("The number of residues must be non-negative.")
    if n_residues == 0:
        return 0

    residue_count = 0
    for position in range(len(sequence)):
        if sequence[position] != GAP:
            residue_count += 1
            if residue_count == n_residues:
                return position + 1

    raise ValueError(f"Sequence has fewer than {n_residues} residues.")


### SPLITTING

def split_at(msa: list[str], split_indexes: list[int]) -> tuple[list[str], list[str], list[int]]:
    """
    Split each sequence of an MSA at its own index.

    Parameters
    msa: list[str]. Alignment to split, one string per sequence.
    split_indexes: list[int]. Split position for each sequence: ``msa[i]`` is cut at ``split_indexes[i]``.

    Returns
    left: list[str]. Left part of each sequence (``msa[i][:split_indexes[i]]``).
    right: list[str]. Right part of each sequence (``msa[i][split_indexes[i]:]``).
    residue_counts: list[int]. Number of residues in each line of the left part.

    Example
    split_at(["ATC-G-G-TT"], [3]) -> (["ATC"], ["-G-G-TT"], [3])
    """
    left: list[str] = []
    right: list[str] = []
    residue_counts: list[int] = []
    for i in range(len(msa)):
        sequence = msa[i]
        split_index = split_indexes[i]
        left_part = sequence[:split_index]
        left.append(left_part)
        right.append(sequence[split_index:])
        residue_counts.append(count_residues(left_part))
    return left, right, residue_counts


def residues_to_indexes(msa: list[str], residue_counts: list[int]) -> list[int]:
    """
    Compute, for each sequence, the split index that keeps a given
    number of residues on the left. Used to split the second parent at the same residues where the
    first parent was split.

    Parameters
    msa: list[str]. Alignment to split (second parent).
    residue_counts: list[int]. Number of residues that must stay on the left part of each
        sequence (as returned by `split_at` on the first parent).

    Returns
    list[int]. Split index for each sequence of `msa`.

    Example
    residues_to_indexes(["ATCGG---TT", "AC-CGCA-TC", "A---CCGATC"], [3, 2, 2])
    >> [3, 2, 5]
    """
    split_indexes: list[int] = []
    for i in range(len(msa)):
        split_indexes.append(index_at_residue(msa[i], residue_counts[i]))
    return split_indexes


### PADDING, MERGING AND CLEANING

def pad_alignment(msa: list[str], side: str) -> list[str]:
    """
    Pad the sequences of a (jagged) MSA part with gaps so that all of
    them have the length of the longest one.

    The left part of a split MSA must be padded on the right ("r") and
    the right part on the left ("l"), so that the residues stay next
    to the split point. This may create gap-only columns, which are
    removed later by `clean_alignment`.

    Parameters
    msa: list[str]. Part of an MSA, possibly with sequences of different lengths.
    side: str. "r" to add the gaps at the end, "l" to add them at the start.

    Returns
    list[str]. New list where every sequence has the same length.

    Raises ValueError
        If `side` is not "l" or "r".
    """
    if side not in ("l", "r"):
        raise ValueError('side must be "l" or "r".')

    max_length = 0
    for sequence in msa:
        if len(sequence) > max_length:
            max_length = len(sequence)

    padded: list[str] = []
    for sequence in msa:
        padding = GAP * (max_length - len(sequence))
        if side == "r":
            padded.append(sequence + padding)
        else:
            padded.append(padding + sequence)
    return padded


def merge(left: list[str], right: list[str]) -> list[str]:
    """
    Join, sequence by sequence, a left part and a right part coming
    from two different MSAs.

    Both parts must have the same number of sequences and each part
    must already be rectangular (see `pad_alignment`).

    Parameters
    left : list[str]. Left part (one string per sequence).
    right : list[str]. Right part (one string per sequence).

    Returns
    list[str]. Merged alignment: ``left[i] + right[i]`` for each sequence.
    """
    merged: list[str] = []
    for i in range(len(left)):
        merged.append(left[i] + right[i])
    return merged


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
    """
    if not msa:
        return []

    kept_columns: list[int] = []
    for column in range(len(msa[0])):
        if not is_all_gaps(msa, column):
            kept_columns.append(column)

    cleaned: list[str] = []
    for sequence in msa:
        new_sequence = ""
        for column in kept_columns:
            new_sequence += sequence[column]
        cleaned.append(new_sequence)
    return cleaned


### CROSSOVER FUNCTION

def generate_offspring(parent1: list[str], parent2: list[str], split_column: int) -> tuple[list[str], list[str]]:
    """
    Apply the one-point crossover operator to two parent MSAs.

    1. `parent1` is split vertically at `split_column`.
    2. `parent2` is split at the same residues of each sequence, which
       usually gives a jagged split.
    3. The parts of `parent2` are padded with gaps to make them
       rectangular.
    4. The halves are swapped, merged and cleaned of gap-only columns.

    Parameters
    parent1: list[str]. First parent MSA (split by column).
    parent2: list[str]. Second parent MSA, with the same sequences as `parent1`.
    split_column : int. Column where `parent1` is split. Values 0 and len(parent1[0])
        produce offspring identical to the parents.

    Returns
    offspring1: list[str]. Left part of `parent1` + right part of `parent2`.
    offspring2: list[str]. Left part of `parent2` + right part of `parent1`.

    Raises ValueError
        If the parents have a different number of sequences or if
        `split_column` is outside [0, len(parent1[0])].
    """
    if len(parent1) != len(parent2):
        raise ValueError("Parents must have the same number of sequences.")
    if not 0 <= split_column <= len(parent1[0]):
        raise ValueError(f"split_column must be between 0 and {len(parent1[0])}.")

    # 1. Split the first parent at a single column
    left1, right1, residue_counts = split_at(parent1, [split_column] * len(parent1))
    # 2. Split the second parent at the same residues
    split_indexes = residues_to_indexes(parent2, residue_counts)
    left2, right2, _ = split_at(parent2, split_indexes)
    # 3. Pad the jagged parts of the second parent
    left2_padded = pad_alignment(left2, "r")
    right2_padded = pad_alignment(right2, "l")
    # 4. Merge and remove gap-only columns
    offspring1 = clean_alignment(merge(left1, right2_padded))
    offspring2 = clean_alignment(merge(left2_padded, right1))
    return offspring1, offspring2


### EXAMPLE TESTS (only run when the file is executed directly)

if __name__ == "__main__":
    align1 = ["ATC-G-G-TT", "A-CCGC-ATC", "A-CCG--ATC"]
    align2 = ["ATCGG---TT", "AC-CGCA-TC", "A---CCGATC"]
    print_msa(align1)
    print()
    print_msa(align2)
    print()

    print(count_residues("--AC--"))          # 2
    print(index_at_residue("--AC--GT", 2))   # 4

    left1, right1, residue_counts = split_at(align1, [3] * len(align1))
    print_msa(left1)
    print_msa(right1)
    print(residue_counts)                    # [3, 2, 2]

    split_indexes = residues_to_indexes(align2, residue_counts)
    print(split_indexes)                     # [3, 2, 5]

    left2, right2, _ = split_at(align2, split_indexes)
    left2_padded = pad_alignment(left2, "r")
    print_msa(left2_padded)
    right2_padded = pad_alignment(right2, "l")
    print_msa(right2_padded)

    left1_right2 = merge(left1, right2_padded)
    print_msa(left1_right2)
    print(is_all_gaps(left1_right2, 2))      # False
    print(is_all_gaps(left1_right2, 3))      # True

    left1_right2_clean = clean_alignment(left1_right2)
    print_msa(left1_right2_clean)
    print()

    offspring1, offspring2 = generate_offspring(align1, align2, 3)
    print_msa(offspring1)
    print()
    print_msa(offspring2)
