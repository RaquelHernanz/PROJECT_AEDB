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

