"""
Genetic algorithm for multiple sequence alignment (MSA) - Phase 1.

Building-block functions for the genetic algorithm (based on SAGA,
Notredame & Higgins, 1996): data reading, population initialization,
scoring, mutation operators and validation.

The crossover operator functions are in their own module (crossover.py) and is
NOT duplicated here: only the shared helpers are imported from it.

An MSA is represented as a list of strings, one per sequence, where
gaps are indicated by '-'. A population is a list of MSAs.
Substitution matrices are nested dictionaries: matrix["A"]["C"] is the
score of aligning A with C.

Functions marked as SKELETON still have to be implemented: their
docstring describes what to do.
"""

import random
import crossover as mc
from crossover import GAP, clean_alignment, print_msa
import msa as msa
import psa as psa


### 3. SCORING

def pair_score(char1: str, char2: str, sub_matrix: dict[str, dict[str, int]], gap: int) -> int:
    """
    Score two characters of the same MSA column.

    (Previously `score_column` in the lab code from class. The gap-gap case was
    added: in an MSA two sequences can both have a gap in the same
    column, and that pair must not be penalized.)

    Parameters
    char1: str. Character of the first sequence (residue or GAP).
    char2: str. Character of the second sequence (residue or GAP).
    sub_matrix : dict[str, dict[str, int]]. Substitution matrix.
    gap : int. Gap penalty (negative value).

    Returns
    int. 0 if both are gaps, `gap` if only one is a gap, otherwise the
        score from `sub_matrix`.
    """
    if char1 == GAP and char2 == GAP:
        return 0
    if char1 == GAP or char2 == GAP:
        return gap
    return sub_matrix[char1][char2]


def score_pairwise(seq1: str, seq2: str, sub_matrix: dict[str, dict[str, int]], gap: int) -> int:
    """
    Score two aligned sequences (two rows of an MSA). It was previously `score_alignment` in the lab code
    from class.

    Parameters
    seq1: str. First aligned sequence (with gaps).
    seq2: str. Second aligned sequence (same length).
    sub_matrix : dict[str, dict[str, int]]. Substitution matrix.
    gap: int. Gap penalty (negative value).

    Returns
    int. Sum of `pair_score` over all columns.
    """
    total = 0
    for i in range(len(seq1)):
        total += pair_score(seq1[i], seq2[i], sub_matrix, gap)
    return total


def sum_of_pairs(msa: list[str], sub_matrix: dict[str, dict[str, int]], gap: int) -> int:
    """
    SKELETON - Compute the sum-of-pairs score of an MSA.

    To implement
    - Add `score_pairwise(msa[i], msa[j], sub_matrix, gap)` for every
      pair of sequences with i < j (each pair counted once).

    Parameters
    msa : list[str]. Alignment to score.
    sub_matrix : dict[str, dict[str, int]]. Substitution matrix.
    gap : int. Gap penalty (negative value).

    Returns
    int. Total score (higher is better).
    """
    raise NotImplementedError


# Optional: sum_of_pairs_affine(msa, sub_matrix, gap_open, gap_extend)
# Same as sum_of_pairs, but distinguishing between opening and
# extending a gap.


### 4. POPULATION INITIALIZATION

def random_msa(sequences: list[str], max_offset: int, rng: random.Random) -> list[str]:
    """
    SKELETON - Create a random MSA from unaligned sequences (SAGA
    approach).

    To implement
    - Add a random number of gaps (between 0 and `max_offset`) at the
      beginning of each sequence.
    - Pad the end of each sequence with gaps until all of them have
      the length of the longest one.
    - Remove gap-only columns with `clean_alignment`.

    Parameters
    sequences: list[str]. Unaligned sequences (no gaps).
    max_offset: int. Maximum number of leading gaps per sequence.
    rng: random.Random. Random number generator (seeded, for reproducibility).

    Returns
    list[str]. A valid MSA of `sequences`.
    """
    raise NotImplementedError


def initialize_population(sequences: list[str], population_size: int,
                          max_offset: int, seed: int) -> list[list[str]]:
    """
    SKELETON - Generate an initial population of random MSAs.

    To implement
    - Create a `random.Random(seed)` generator.
    - Call `random_msa` `population_size` times with that generator.

    Parameters
    sequences: list[str]. Unaligned sequences (no gaps).
    population_size: int. Number of MSAs in the population.
    max_offset: int. Maximum number of leading gaps per sequence.
    seed: int. Seed for the random generator.

    Returns
    list[list[str]]. Population (list of MSAs).
    """
    raise NotImplementedError


### 5. MUTATION OPERATORS

def find_gap_blocks(sequence: str) -> list[tuple[int, int]]:
    """
    SKELETON - Find the blocks of consecutive gaps in an aligned
    sequence.

    To implement
    - Scan the sequence and record where each run of GAP starts and
      ends.
    - Decide (and document) whether `end` is inclusive or exclusive.

    Parameters
    sequence: str. A single line of an MSA.

    Returns
    list[tuple[int, int]]. List of (start, end) positions of each gap block.
    """
    raise NotImplementedError


def insert_gap_block(msa: list[str], rng: random.Random) -> list[str]:
    """
    SKELETON - Mutation 1: insert a gap block of random size at a
    random position of a single sequence.

    To implement
    - Choose a sequence, a position and a block size at random.
    - Insert the block in that sequence.
    - Add the same number of gaps to the other sequences (at the end
      or at a random position) so the MSA stays rectangular.
    - Finish with `clean_alignment`.

    Parameters
    msa : list[str]. Alignment to mutate (must not be modified).
    rng : random.Random. Random number generator.

    Returns
    list[str]. New, valid MSA.
    """
    raise NotImplementedError


def delete_gap_block(msa: list[str], rng: random.Random) -> list[str]:
    """
    SKELETON - Mutation 2: delete a random gap block.

    To implement
    - Choose a sequence that has gaps (use `find_gap_blocks`).
    - Remove one of its gap blocks.
    - Restore the length: add gaps at the end of that sequence (or
      remove columns from the others).
    - Finish with `clean_alignment`.
    - If no sequence has gaps, return an unchanged copy.

    Parameters
    msa : list[str]. Alignment to mutate (must not be modified).
    rng : random.Random. Random number generator.

    Returns
    list[str]. New, valid MSA.
    """
    raise NotImplementedError


def move_gap_block(msa: list[str], rng: random.Random) -> list[str]:
    """
    SKELETON - Mutation 7: move a gap block one position in either
    direction.

    To implement
    - Choose a sequence with gaps and one of its gap blocks.
    - Swap the block with the neighbouring residue on the left or on
      the right (the length does not change).
    - Handle the edges: the block cannot move outside the sequence.
    - Finish with `clean_alignment`.
    - If no sequence has gaps, return an unchanged copy.

    Parameters
    msa : list[str]. Alignment to mutate (must not be modified).
    rng : random.Random. Random number generator.

    Returns
    list[str]. New, valid MSA.
    """
    raise NotImplementedError


def mutate(msa: list[str], rng: random.Random) -> list[str]:
    """
    SKELETON - Apply one mutation operator chosen at random.

    To implement
    - Keep a list of the implemented operators.
    - Choose one with `rng` and apply it to `msa`.

    Parameters
    msa : list[str]. Alignment to mutate.
    rng : random.Random. Random number generator.

    Returns
    list[str]. New, valid MSA.
    """
    raise NotImplementedError


# Rules for every mutation operator:
# - return a new MSA, never modify the original;
# - never change the order of the residues;
# - keep all rows with the same length;
# - finish with clean_alignment;
# - return an unchanged copy if the chosen sequence has no gaps.


### 6. VALIDATION (for testing)

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


### TESTS (only run when the file is executed directly)

if __name__ == "__main__":
    dna_matrix = psa.make_sub_matrix(2, -2, "ACGT")

    # pair_score: gap-gap must be 0
    print(pair_score("A", "A", dna_matrix, -3))   # 2
    print(pair_score("A", "-", dna_matrix, -3))   # -3
    print(pair_score("-", "-", dna_matrix, -3))   # 0

    # score_pairwise: same example as the lab -> 4
    print(score_pairwise("-CAGTGCATG-ACATA", "TCAG-GC-TCTACAGA", dna_matrix, -3))

    # is_valid_msa
    original = ["ATCGGTT", "ACCGCATC", "ACCGATC"]
    print(is_valid_msa(["ATC-G-G-TT", "A-CCGC-ATC", "A-CCG--ATC"], original))  # True
    print(is_valid_msa(["ATC-GG-TT", "A-CCGC-ATC", "A-CCG--ATC"], original))   # False (lengths)
    print(is_valid_msa(["-ATCGGTT", "-ACCGCATC", "-ACCGATC"], original))      # False (lengths)
    print(is_valid_msa(["-ATCGGTT-", "-ACCGCATC", "-ACCGATC-"], original))    # False (gap-only column)

    # Uncomment when the files are in the same folder:
    # identifiers, sequences = read_fasta("cytochromes.fa")
    # blosum62 = read_sub_matrix("blosum62.mat")
