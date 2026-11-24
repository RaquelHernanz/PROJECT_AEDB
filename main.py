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

import random as r
import initialization as init
import scoring as sc
import mutations as mt
import crossover as cso
from crossover import GAP, clean_alignment, print_msa
import utils

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
    dna_matrix = utils.make_sub_matrix(2, -2, "ACGT")

    # is_valid_msa
    original = ["ATCGGTT", "ACCGCATC", "ACCGATC"]
    print(is_valid_msa(["ATC-G-G-TT", "A-CCGC-ATC", "A-CCG--ATC"], original))  # True
    print(is_valid_msa(["ATC-GG-TT", "A-CCGC-ATC", "A-CCG--ATC"], original))   # False (lengths)
    print(is_valid_msa(["-ATCGGTT", "-ACCGCATC", "-ACCGATC"], original))      # False (lengths)
    print(is_valid_msa(["-ATCGGTT-", "-ACCGCATC", "-ACCGATC-"], original))    # False (gap-only column)

    # Uncomment when the files are in the same folder:
    # identifiers, sequences = read_fasta("cytochromes.fa")
    # blosum62 = read_sub_matrix("blosum62.mat")
