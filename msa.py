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

dna_matrix = make_sub_matrix(2, -2, "ACGT")
blosum62 = read_sub_matrix("blosum62.mat")

def score_column(char1: str, char2: str, sub_matrix: dict, gap: int) -> int:
    if char1 == "-" or char2 == "-":
        return gap
    else:
        return sub_matrix[char1][char2]

def score_alignment(seq1: str, seq2: str, sub_matrix: dict, gap: int) -> int:
    score = 0
    for i in range(len(seq1)):
        score += score_column(seq1[i], seq2[i], sub_matrix, gap)
    return score

DIAG, UP, LEFT = 1,2,3
def best_move(diag: int, up: int, left: int) -> int:
    if diag > up:
        if diag > left:
            return DIAG
        else:
            return LEFT
    else:
        if up > left:
            return UP
        else:
            return LEFT

def print_matrix(matrix: list[list[int]]) -> None:
    for l in matrix:
        print(l)

def needleman_wunsch(seq1: str, seq2: str, sub_matrix: dict, gap: int) -> tuple[list[list[int]], list[list[int]]]:
    S = [[0]]
    T = [[0]]

    for i in range(1, len(seq1)+1):
        S.append([gap*i])
        T.append([UP])
    for j in range(1, len(seq2)+1):
        S[0].append(gap*j)
        T[0].append(LEFT)

    # here i,j are the indices of the characters in seq1 and seq2, respectively
    # which are one less than the index of the row/column in the matrix S and T
    for i in range(0, len(seq1)):
        for j in range(0, len(seq2)):
            diag = S[i][j] + score_column(seq1[i], seq2[j], sub_matrix, gap)
            up = S[i][j+1] + gap
            left = S[i+1][j] + gap
            S[i+1].append(max(diag, up, left))
            T[i+1].append(best_move(diag, up, left))
    return S, T

def traceback_alignment(T: list[list[int]], seq1: str, seq2: str) -> tuple[str, str]:
    aligned1 = ""
    aligned2 = ""
    i = len(seq1)
    j = len(seq2)
    while i > 0 or j > 0:
        if T[i][j] == DIAG:
            aligned1 = seq1[i-1] + aligned1
            aligned2 = seq2[j-1] + aligned2
            i -= 1
            j -= 1
        elif T[i][j] == UP:
            aligned1 = seq1[i-1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1
    return aligned1, aligned2
