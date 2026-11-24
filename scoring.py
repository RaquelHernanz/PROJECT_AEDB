import random as r
from crossover import GAP, clean_alignment, print_msa
import utils

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

### TESTS (only run when the file is executed directly)

if __name__ == "__main__":
    dna_matrix = utils.make_sub_matrix(2, -2, "ACGT")

    # pair_score: gap-gap must be 0
    print(pair_score("A", "A", dna_matrix, -3))   # 2
    print(pair_score("A", "-", dna_matrix, -3))   # -3
    print(pair_score("-", "-", dna_matrix, -3))   # 0

    # score_pairwise: same example as the lab -> 4
    print(score_pairwise("-CAGTGCATG-ACATA", "TCAG-GC-TCTACAGA", dna_matrix, -3))