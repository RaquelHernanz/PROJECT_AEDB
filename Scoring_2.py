import utils
from utils import GAP

#-------------------------------------------- PAIR SCORE -------------------------------------------

def pair_score(char1: str, char2: str, sub_matrix: dict[str, dict[str, int]], gap: int) -> int:
    """
    Score two characters of the same MSA column.
    (Previously `score_column` in the lab code from class. The gap-gap case was added: in 
    an MSA two sequences can both have a gap in the same column, and that pair must not be penalized.)

    Parameters: 
    - char1: str. Character of the first sequence (residue or GAP).
    - char2: str. Character of the second sequence (residue or GAP).
    - sub_matrix : dict[str, dict[str, int]]. Substitution matrix.
    - gap : int. Gap penalty (negative value).

    Returns: int. 0 if both are gaps, `gap` if only one is a gap, otherwise the score from `sub_matrix`.
    """

    if char1 == GAP and char2 == GAP:
        return 0
    if char1 == GAP or char2 == GAP:
        return gap
    return sub_matrix[char1][char2]

#-------------------------------------------- SCORE PAIRWISE -------------------------------------------

def score_pairwise(seq1: str, seq2: str, sub_matrix: dict[str, dict[str, int]], gap: int) -> int:
    """
    Score two aligned sequences (two rows of an MSA). It was previously 
    `score_alignment` in the lab code from class.

    Parameters
    - seq1: str. First aligned sequence (with gaps).
    - seq2: str. Second aligned sequence (same length).
    - sub_matrix : dict[str, dict[str, int]]. Substitution matrix.
    - gap: int. Gap penalty (negative value).

    Returns: int. Sum of `pair_score` over all columns.
    """
    #Check if the sequences have the same length
    if len(seq1) != len(seq2):
        raise ValueError("Aligned sequences must have the same length.")
    
    total = 0
    for i in range(len(seq1)):
        total += pair_score(seq1[i], seq2[i], sub_matrix, gap)


    return total

#-------------------------------------------- SUM OF PAIRS -------------------------------------------

def sum_of_pairs(msa: list[str], sub_matrix: dict[str, dict[str, int]], gap: int) -> int:
    """
    Compute the sum-of-pairs score of an MSA.

    Parameters
    - msa : list[str]. Alignment to score.
    - sub_matrix : dict[str, dict[str, int]]. Substitution matrix.
    - gap : int. Gap penalty (negative value).

    Returns: 
    - int. Total score (higher is better).

     To implement
        - Add `score_pairwise(msa[i], msa[j], sub_matrix, gap)` for every
          pair of sequences with i < j (each pair counted once).
    """
    #Confirm that the MSA has at least two sequences
    if len(msa) < 2:
        raise ValueError("MSA must contain at least two sequences.")

    #Confirm that all sequences have the same length
    for seq in msa:
        if len(seq) != len(msa[0]):
            raise ValueError("All aligned sequences must have the same length.")

    total = 0 
    for i in range(len(msa)):
        for j in range(i + 1, len(msa)):
            total += score_pairwise(msa[i], msa[j], sub_matrix, gap)

    return total

if __name__ == "__main__":
    
    dna_matrix = utils.make_sub_matrix(2, -2, "ACGT")
    msa = ["ACGT", "A-GT", "AC-T"] #4

    # pair_score
    print(pair_score("A", "A", dna_matrix, -3))     # 2
    print(pair_score("A", "C", dna_matrix, -3))     # -2
    print(pair_score("A", "-", dna_matrix, -3))     # -3
    print(pair_score("-", "-", dna_matrix, -3))     # 0

    # score_pairwise
    print(score_pairwise("ACGT", "ACGT", dna_matrix, -3))   # 8
    print(score_pairwise("AC-T", "A-GT", dna_matrix, -3))   # -2  (2 -3 -3 +2)

    # sum_of_pairs
    msa = ["ACGT", "A-GT", "AC-T"]
    print(sum_of_pairs(msa, dna_matrix, -3))                # 4   (3 + 3 - 2)
    print(sum_of_pairs(["A", "-", "-"], dna_matrix, -3))    # -6  (gap-gap cost 0)
