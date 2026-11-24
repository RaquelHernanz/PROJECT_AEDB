import random as r
from crossover import GAP, clean_alignment, print_msa
import utils

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


def insert_gap_block(msa: list[str], rng: r.Random) -> list[str]:
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


def delete_gap_block(msa: list[str], rng: r.Random) -> list[str]:
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


def move_gap_block(msa: list[str], rng: r.Random) -> list[str]:
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


def mutate(msa: list[str], rng: r.Random) -> list[str]:
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
