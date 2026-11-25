
import random as r
import utils
from utils import GAP

#------------------------------ RANDOM MSA -----------------------------------------------------------

def random_msa(sequences: list[str], max_offset: int, rng: r.Random) -> list[str]:
  """
  Create the initial population of MSAs (Multiple Sequence Alignments) (SAGA approach).

    Parameters:
    - sequences: list[str]. Unaligned sequences (no gaps).
    - max_offset: int. Maximum number of leading gaps per sequence.
    - rng: random.Random. Random number generator (seeded, for reproducibility).

    Returns: a list[str] of random MSAs from the sequences present in the FASTA file. 
    - Invalid sequences do not occur (eg. columns containing only gaps).

    Example 
    with offsets of 1, 0, 2:
        ACGTAC   ->  -ACGTAC--
        AGTC     ->  AGTC-----
        ACCGTTA  ->  --ACCGTTA

  """
  #In case of max_offset is negative
  if max_offset < 0:
     raise ValueError("max_offset must be non-negative.")
  
  # 1. Add a random number of gaps (between 0 and 'max_offset', both included) at the beginning of each sequence.
  aligned_sequences = []

  for sequence in sequences:
    offset = rng.randint(0,max_offset)    #Picks a random integer between 0 and the defined maximum value (max_offset)
    new_sequence = GAP * offset + sequence  # Joins all the hyphens (- or GAP) with the original sequence.
    aligned_sequences.append(new_sequence)

    # 2. Finds the longest sequence and adds gaps to the end of the others until they all match its length: 
    #  - max() : finds the maximum length
    #  - len() : length of each sequence

  max_length = max(len(sequence) for sequence in aligned_sequences)
  padded_sequences = []

  #Subtracts the sequence length from the maximum length to find the number of gaps needed
  for sequence in aligned_sequences:
    gaps_needed = max_length - len(sequence)    
    new_sequence = sequence + "-" * gaps_needed
    padded_sequences.append(new_sequence)
    
  #3. Remove gap-only columns with 'clean_alignment'. (e.g. if every offset is >= 1, the first 
  # column only has gaps).
  aligned_sequences = utils.clean_alignment(padded_sequences)

  return aligned_sequences

#--------------------- POPULATION INITIALIZATION -------------------------------------------------

def initialize_population(sequences: list[str], population_size: int, max_offset: int, seed: int) -> list[list[str]]:
    """
    Generate an initial population of random MSAs. 

    Parameters:
    - sequences: list[str]. Unaligned sequences (no gaps).
    - population_size: int. Number of MSAs in the population.
    - max_offset: int. Maximum number of leading gaps per sequence.
    - seed: int. Seed for the random generator.

    Return: list[list[str]]. Population (list of MSAs)."""

  
    #1. Create a `random.Random(seed)` generator.
    rng = r.Random(seed)
    population = []

    if population_size < 1:
          raise ValueError("population_size must be at least 1.")

    #2. Call `random_msa`population_size` times with that generator.    
    for _ in range(population_size):
      msa = random_msa(sequences, max_offset, rng)
      population.append(msa)

    return population

if __name__ == "__main__":
    seqs = ["ACGTAC", "AGTC", "ACCGTTA"]

    # random_msa: valid MSA (rows, lenght, without column gaps, intact residues)
    msa = random_msa(seqs, 3, r.Random(0))
    utils.print_msa(msa)
    print(utils.is_valid_msa(msa, seqs)) # True

    # max_offset = 0 -> without initial gaps
    msa0 = random_msa(seqs, 0, r.Random(0))
    print(all(not row.startswith(GAP) for row in msa0[:1])) # True

    # initialize_population: sequences, size, offset and seed reproducibility
    pop = initialize_population(seqs, 50, 3, seed=1)
    print(len(pop))                                            # 50
    print(all(utils.is_valid_msa(m, seqs) for m in pop))       # True
    print(pop == initialize_population(seqs, 50, 3, seed=1))   # True  (same seed)
    print(pop == initialize_population(seqs, 50, 3, seed=2))   # False (other seed)

