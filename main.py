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
import Population_initialization_1 as init
import Scoring_2 as sc
import Mutations_4 as mt
import Crossover_3 as cso
import utils
from utils import GAP

### TESTS (only run when the file is executed directly)

if __name__ == "__main__":
    dna_matrix = utils.make_sub_matrix(2, -2, "ACGT")

    # is_valid_msa
    original = ["ATCGGTT", "ACCGCATC", "ACCGATC"]
    print(utils.is_valid_msa(["ATC-G-G-TT", "A-CCGC-ATC", "A-CCG--ATC"], original))  # True
    print(utils.is_valid_msa(["ATC-GG-TT", "A-CCGC-ATC", "A-CCG--ATC"], original))   # False (lengths)
    print(utils.is_valid_msa(["-ATCGGTT", "-ACCGCATC", "-ACCGATC"], original))      # False (lengths)
    print(utils.is_valid_msa(["-ATCGGTT-", "-ACCGCATC", "-ACCGATC-"], original))    # False (gap-only column)

    # Uncomment when the files are in the same folder:
    # identifiers, sequences = read_fasta("cytochromes.fa")
    # blosum62 = read_sub_matrix("blosum62.mat")
