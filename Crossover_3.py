"""
Crossover operator for multiple sequence alignments (MSA).

This module contains the functions needed to perform the crossover
operator: two parent MSAs are split at a chosen column (the second 
parent at the equivalent residues) and the halves are recombined into two offspring.

Note: a MSA is represented as a list of strings, one per sequence, where
gaps are indicated by '-'. All strings of a valid MSA have the same
length and there are no columns made only of gaps.
"""


#-------------------------------------------- IMPORTS -------------------------------------------
import random 
from Scoring_2 import sum_of_pairs
import utils
from utils import GAP



#-------------------------------------------- AUXILIAR FUNCTIONS -------------------------------------------

def count_residues(sequence: str) -> int:
    """
    Count the residues (non-gap characters) in an aligned sequence.

    Parameters
    - sequence: str. A single line of an MSA, possibly containing gaps.

    Returns
    - int. Number of characters different from the gap symbol.

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
    - sequence: str. A single line of an MSA, possibly containing gaps.
    - n_residues: int. Number of residues the left part must contain.

    Returns
    - int. Index such that ``sequence[:index]`` contains exactly
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


#-------------------------------------------- SPLITTING -------------------------------------------

def split_at(msa: list[str], split_indexes: list[int]) -> tuple[list[str], list[str], list[int]]:
    """
    Split each MSA sequence at its own index (previously obtained from index_at_residue). 
    The idea is to obtain the parent sequences (left and right) from spliting each MSA sequence, 
    and how many residues are in total on the left (l) parent sequences 

    Parameters
    - msa: list[str]. Alignment to split, one string per sequence.
    - split_indexes: list[int]. Split position for each sequence: ``msa[i]`` is cut at ``split_indexes[i]``.

    Returns
    - left: list[str]. Left part of each sequence (``msa[i][:split_indexes[i]]``).
    - right: list[str]. Right part of each sequence (``msa[i][split_indexes[i]:]``).
    - residue_counts: list[int]. Number of residues in each line of the left part.

    Example
    split_at(["ATC-G-G-TT"], [3]) -> (["ATC"], ["-G-G-TT"], [3])

    left        right     residues
    ATC        -G-G-TT    3
    A-C        CGC-ATC    2
    A-C        CG--ATC    2

    residue_counts = [3, 2, 2]

    """
    left: list[str] = []
    right: list[str] = []
    residue_counts: list[int] = []

    #Verification: Ensures the data passed to the function is consistent before processing
    if len(msa) != len(split_indexes):
        raise ValueError("msa and split_indexes must have the same length.")


    for i in range(len(msa)):
        sequence = msa[i]

        #Choose the splitting index previously obtained
        split_index = split_indexes[i]

        #Left: From Seq[0] to Seq[index-1], like ATC
        left_part = sequence[:split_index]
        left.append(left_part)

        #From Seq[index] to Seq[n-1], like -G-G-TT
        right.append(sequence[split_index:])
        residue_counts.append(count_residues(left_part)) #It would be 3

    return left, right, residue_counts


def residues_to_indexes(msa: list[str], residue_counts: list[int]) -> list[int]:
    """
    Compute for each sequence the function index_at_residue, obtaining the split index 
    that keeps a given number of residues on the left. Sometimes the cuts can be irregular
    or jagged (the index positions obtained are not equal between sequences)

    Parameters
    - msa: list[str]. Alignment to split (second parent).
    - residue_counts: list[int]. Number of residues that must stay on the left part of each
        sequence (as returned by `split_at` on the first parent).

    Returns: list[int]. Split index for each sequence of `msa`.

    Example
    residues_to_indexes(["ATCGG---TT", "AC-CGCA-TC", "A---CCGATC"], [3, 2, 2])
    >> [3, 2, 5]

    ATCGG---TT  → 3 residues → index 3   (ATC | GG---TT)
    AC-CGCA-TC  → 2 residues → index 2   (AC  | -CGCA-TC)
    A---CCGATC  → 2 residues → index 5   (A---C | CGATC)
    """

    #Step 1: Ensures the data passed to the function is consistent before processing
    if len(msa) != len(residue_counts):
        raise ValueError("msa and residue_counts must have the same length.")


    #Step2: Finds the exact index where each sequence from the second parent must be cut 
    # to maintain the exact same number of residues that remained in the left half of the first parent.
    split_indexes: list[int] = []
    for i in range(len(msa)):
        split_indexes.append(index_at_residue(msa[i], residue_counts[i]))
    return split_indexes



#---------------------------- PADDING, MERGING and CLEANING -------------------------------------------

def pad_alignment(msa: list[str], side: str) -> list[str]:
    """
    Pad the parent sequences (in other words, it makes their length equal), since some 
    of the them maybe jagged. To do this the function must add gaps so that all have length 
    of the longest one. 

    The left part of a split MSA must be padded on the right ("r") and
    the right part on the left ("l"), so that the residues stay next
    to the split point. This may create gap-only columns, which are
    removed later by `clean_alignment`.

    Parameters
    - msa: list[str]. Part of an MSA, possibly with sequences of different lengths.
    - side: str. "r" to add the gaps at the end, "l" to add them at the start.

    Returns: list[str]. New list where every sequence has the same length.
             - Raises ValueError If `side` is not "l" or "r".

    Example:

    (ATC | GG---TT)
    (AC  | -CGCA-TC)
    (A---C | CGATC)

    left2 → "r"      right2 → "l"
    ATC--            -GG---TT
    AC---            -CGCA-TC
    A---C            ---CGATC

    """
    if side not in ("l", "r"):
        raise ValueError('side must be "l" or "r".')

    #What is the maximun length of the right/left sequences
    max_length = 0
    for sequence in msa:
        if len(sequence) > max_length:
            max_length = len(sequence) 
    
    #We add the gaps to the sequences
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
    must already be rectangular (see the pad_alignment function).

    Parameters:
    - left : list[str]. Left part (one string per sequence).
    - right : list[str]. Right part (one string per sequence).

    Returns: list[str]. Merged alignment: ``left[i] + right[i]`` for each sequence.

    Example
    ATC--      -GG---TT
    AC---      -CGCA-TC
    A---C      ---CGATC

    ATC-GG---TT
    A-C-CGCA-TC
    A-C---CGATC

    """
    #Step 1: Ensures the data passed to the function is consistent before processing
    if len(left) != len(right):
        raise ValueError("Left and right parts must have the same number of sequences.")


    #Step2: Joins the left half of one alignment with the right half of the other, 
    #       sequence by sequence, creating the resulting new alignment (still contains the gap sequences).

    merged: list[str] = []
    #Iterate through the list of parent sequences
    for i in range(len(left)):
        merged.append(left[i] + right[i])
    return merged


def generate_offspring(parent1: list[str], parent2: list[str], split_column: int) -> tuple[list[str], list[str]]:
    """
    Apply the one-point crossover operator to two parent MSAs.

    1. `parent1` is split vertically at `split_column`.
    2. `parent2` is split at the same residues of each sequence, which
       usually gives a jagged split.
    3. The parts of `parent2` are padded with gaps to make them
       rectangular.
    4. The halves are swapped, merged and cleaned of gap-only columns.

    Parameters:
    - parent1: list[str]. First parent MSA (split by column).
    - parent2: list[str]. Second parent MSA, with the same sequences as `parent1`.
    - split_column : int. Column where `parent1` is split. Values 0 and len(parent1[0])
        produce offspring identical to the parents.

    Returns:
    - offspring1: list[str]. Left part of `parent1` + right part of `parent2`.
    - offspring2: list[str]. Left part of `parent2` + right part of `parent1`.

    Raises ValueError: If the parents have a different number of sequences or if
                       `split_column` is outside [0, len(parent1[0])].

    Example:
    offspring1 = left(p1) + right(p2)  offspring2 = left(p2) + right(p1)
    ATCGG---TT                         ATC--G-G-TT
    A-CCGCA-TC                         AC--CGC-ATC
    A-C--CGATC                         A--CCG--ATC
    """

    ############## Verifications ##################
    # 1. Check that parents are not empty
    if not parent1 or not parent2:
        raise ValueError("Parents cannot be empty.")

    # 2. Check that parents have the same number of sequences
    if len(parent1) != len(parent2):
        raise ValueError("Parents must have the same number of sequences.")

    # 3. Check that both parents are valid rectangular MSAs- all sequences must have the same length.
    parent1_length = len(parent1[0])
    parent2_length = len(parent2[0])


    for seq in parent1:
        if len(seq) != parent1_length:
            raise ValueError("Parent1 is not a valid rectangular MSA.")
    for seq in parent2:
        if len(seq) != parent2_length:
            raise ValueError("Parent2 is not a valid rectangular MSA.")

    """5. Checks whether both parents actually contain the exact same original sequences, 
    in the same order, ignoring gaps:

    Ex. Two alignments are compatible, because, when removing the gaps, the sequences are identical
        parent1 = ["ATC-G", "A-CCG"]
        parent2 = ["AT-CG", "AC-CG"]

    """
    for i in range(len(parent1)):
        seq1 = parent1[i].replace(GAP, "")
        seq2 = parent2[i].replace(GAP, "")
        
        if seq1 != seq2:
            raise ValueError("Parents must contain the same ungapped sequences.")


    #6.Checks the variable type
    if type(split_column) != int:
        raise ValueError("split_column must be an integer.")

    if not 0 <= split_column <= len(parent1[0]):
        raise ValueError(f"split_column must be between 0 and {len(parent1[0])}.")
    

    ############# Appling the function ##########

    # 1. Split the first parent sequences at a single column in common
    left1, right1, residue_counts = split_at(parent1, [split_column] * len(parent1))
    #We need the exact number of residues for the parents2 left and right split

    # 2. Split the second parent with the same residues number, but first obtain the indexes for that split. 
    split_indexes = residues_to_indexes(parent2, residue_counts)
    left2, right2, _ = split_at(parent2, split_indexes)

    # 3. Pad the jagged parts of the second parent
    left2_padded = pad_alignment(left2, "r")
    right2_padded = pad_alignment(right2, "l")

    # 4. Merge and remove gap-only columns
    offspring1 = utils.clean_alignment(merge(left1, right2_padded))
    offspring2 = utils.clean_alignment(merge(left2_padded, right1))
    return offspring1, offspring2



#-------------------------------------------- CROSSOVER -------------------------------------------

def crossover(parent1: list[str], parent2: list[str], sub_matrix: dict[str, dict[str, int]], gap: int) -> list[str]:
  """
    Apply one-point crossover at a random column and return the highest-scoring offspring.
    """
  # 1. Choose a random split column (avoiding the ends 0 and max length)
  msa_length = len(parent1[0])
  if msa_length > 1:
    split_column = random.randint(1, msa_length - 1)
  else:
    split_column = 0

  # 2. Generate the two offspring
  offspring1, offspring2 = generate_offspring(parent1, parent2, split_column)

  # 3. Evaluate both offspring and return only the best one (highest score)
  score1 = sum_of_pairs(offspring1, sub_matrix, gap)
  score2 = sum_of_pairs(offspring2, sub_matrix, gap)

  if score1 >= score2:
    return offspring1
  else: 
      return offspring2

#-------------------------------------------- TESTS -------------------------------------------

if __name__ == "__main__":
    align1 = ["ATC-G-G-TT", "A-CCGC-ATC", "A-CCG--ATC"]
    align2 = ["ATCGG---TT", "AC-CGCA-TC", "A---CCGATC"]
    utils.print_msa(align1)
    print()
    utils.print_msa(align2)

    print(count_residues("--AC--"))          # 2
    print(index_at_residue("--AC--GT", 2))   # 4

    left1, right1, residue_counts = split_at(align1, [3] * len(align1))
    utils.print_msa(left1)
    utils.print_msa(right1)
    print(residue_counts)                    # [3, 2, 2]

    split_indexes = residues_to_indexes(align2, residue_counts)
    print(split_indexes)                     # [3, 2, 5]

    left2, right2, _ = split_at(align2, split_indexes)
    left2_padded = pad_alignment(left2, "r")
    utils.print_msa(left2_padded)
    right2_padded = pad_alignment(right2, "l")
    utils.print_msa(right2_padded)

    left1_right2 = merge(left1, right2_padded)
    utils.print_msa(left1_right2)
    print(utils.is_all_gaps(left1_right2, 2))      # False
    print(utils.is_all_gaps(left1_right2, 3))      # True

    left1_right2_clean = utils.clean_alignment(left1_right2)
    utils.print_msa(left1_right2_clean)
    print()

    offspring1, offspring2 = generate_offspring(align1, align2, 3)
    utils.print_msa(offspring1)
    print()
    utils.print_msa(offspring2)

