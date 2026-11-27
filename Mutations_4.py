import random as r
import utils
from utils import GAP, clean_alignment, print_msa, is_valid_msa

#-------------------------------------------- MUTATION OPERTATORS -------------------------------------------

def find_gap_blocks(sequence: str) -> list[tuple[int, int]]:
    """
    Find the blocks of consecutive gaps in an aligned sequence.

    Parameters
    - sequence: str. A single line of an MSA.

    Returns
    - list[tuple[int, int]]. List of (start, end) positions of each gap block.

    Example
    seq1:  AC--GT  (2,4)
    seq2:  A-CGGT  (1,2)
    seq3:  -ACGT-  (0,1) and (5,6)

    """
    gap_blocks: list[tuple[int, int]] = []
    start: int = None
    for i in range(len(sequence)):
        if sequence[i] == GAP and start == None:
            start = i
        elif sequence[i] != GAP and start != None:
            gap_blocks.append((start, i))
            start = None

    if start is not None:
        gap_blocks.append((start, len(sequence)))

    return gap_blocks


def insert_gap_block(msa: list[str], rng: r.Random) -> list[str]:
    """
    Insert a gap block of random size at a random position of a 
    single sequence.

    Parameters
    - msa : list[str]. Alignment to mutate (must not be modified).
    - rng : random.Random. Random number generator.

    Returns
    - list[str]. New, valid MSA.

    Example
      row 2:  -AC + -- + GT-   →  -AC--GT-
      rest:   +"--" al final   →  AC--GT--   A-CGGT--
    """

    #We are going to copy the msa sequences
    copy: list[str] = list(msa)
    #We chose a random sequence
    i: int = rng.randint(0,len(msa)-1)
    # And from that sequence we chose a random position
    position: int = rng.randint(0, len(msa[i]))

    #Number of GAPS to add randomly in each iteration of the function, 
    # between 1-10, not too large
    k: int = rng.randint(1,10) 

    #Add the GAPS in the sequences
    gaps = GAP*k
    copy[i] = copy[i][:position] + gaps + copy[i][position:]

    for j in range(len(copy)) :
        #We dont want to add again the GAPS from the sequence i
        if j!= i:
            copy[j] = copy[j] + gaps

    return clean_alignment(copy)

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
    return msa


def move_gap_block(msa: list[str], rng: r.Random) -> list[str]:
    """
    Move a gap block one position in either
    direction.

    Parameters
    - msa : list[str]. Alignment to mutate (must not be modified).
    - rng : random.Random. Random number generator.

    Returns
    list[str]. New, valid MSA.

    """

    copy: list[str] = list(msa)
    block_with_gaps: list[int] = []

    for j in range(len(msa)):
        if find_gap_blocks(msa[j]) != []:
            block_with_gaps.append(j)

    if (block_with_gaps == []):
        return copy

    #Chose a random position inside block_with_gaps
    pos: int = rng.randint(0, len(block_with_gaps)-1)
    #And take the real sequence index stored in that position
    i: int = block_with_gaps[pos]
    #Obtain all the GAPS positions inside a sequence
    block_positions: list[tuple[int, int]] = find_gap_blocks(copy[i])
    # K will be used for choosing a random GAP position from the
    # sequence
    k: int = rng.randint(0, len(block_positions)-1)
    (start, end) = block_positions[k]



    #We will move randomly to the left or right the GAP
    direction: str = ""
    x: int = rng.randint(0,1) #If 0 = LEFT, if 1 = RIGHT
    if x == 0:
        direction = "LEFT"
    else:
        direction = "RIGHT"

    """
    We should ckeck if the movement is posible:

    s = 0 → nothing to the left to move because it is the end
    of the sequence (without this check, seq[start-1] = seq[-1] 
    would silently grab the last character)

    with this ckeck: direction ← right
    We can move to the right instead
    """

    if (direction == "LEFT" and start == 0):
        direction = "RIGHT"
    elif (direction == "RIGHT" and end == len(copy[i])):
        direction = "LEFT" 

    seq_copy: str = copy[i]
    new_seq: str = ""


    """
    How the movement works: 
    left:     ""  + -  + A  + CGGT  →  -ACGGT
              (after de start-1) (block) (residue start-1) (desde e)

    right:    A   + C  + -  + GGT   →  AC-GGT
              (before the start) (residue e) (block) (desde end+1)
    """
    if (direction == "LEFT"):
        before = seq_copy[:start-1]
        residue = seq_copy[start-1]
        block = seq_copy[start:end]
        rest = seq_copy[end:]
        new_seq = before + block + residue + rest

    else:
        before = seq_copy[:start]
        block = seq_copy[start:end]
        residue = seq_copy[end]
        rest = seq_copy[end+1:]
        new_seq = before + residue + block + rest

    copy[i] = new_seq

    return clean_alignment(copy)
        



def mutate(msa: list[str], rng: r.Random) -> list[str]:
    """
    Apply one mutation operator chosen at random.

    Parameters
    - msa : list[str]. Alignment to mutate.
    - rng : random.Random. Random number generator.

    Returns
    list[str]. New, valid MSA.
    """
    #We choose randomly the mutation we apply
    # We can reuse the rng object
    i: int = rng.randint(0,2)
    if (i ==0):
        return insert_gap_block(msa, rng)
    elif (i ==1):
        return delete_gap_block(msa, rng)
    else:
        return move_gap_block(msa, rng)


# Rules for every mutation operator:
# - return a new MSA, never modify the original;
# - never change the order of the residues;
# - keep all rows with the same length;
# - finish with clean_alignment;
# - return an unchanged copy if the chosen sequence has no gaps.

if __name__ == "__main__":

    msa = ["AC--GT",
           "A-CGGT",
           "-ACGT-"]
    originales = ["ACGT", "ACGGT", "ACGT"]
    rng = r.Random(1) #For this we have a seed of 1 for this test

    # find_gap_blocks
    print(find_gap_blocks("--AC-G---"))   # expected: [(0, 2), (4, 5), (6, 9)]
    print(find_gap_blocks("ACGT"))        # expected: []

    print("\nOriginal:")
    print_msa(msa)

    print("\ninsert_gap_block:")
    resultado = insert_gap_block(msa, rng)
    print_msa(resultado)
    print("¿Valid?", is_valid_msa(resultado, originales))

    print("\ndelete_gap_block:")
    resultado = delete_gap_block(msa, rng)
    print_msa(resultado)
    print("¿Valid?", is_valid_msa(resultado, originales))

    print("\nmove_gap_block:")
    resultado = move_gap_block(msa, rng)
    print_msa(resultado)
    print("¿Valid?", is_valid_msa(resultado, originales))

    print("\nmutate:")
    resultado = mutate(msa, rng)
    print_msa(resultado)
    print("¿Valid?", is_valid_msa(resultado, originales))

    print("\n¿Original without changes?", msa == ["AC--GT", "A-CGGT", "-ACGT-"])
