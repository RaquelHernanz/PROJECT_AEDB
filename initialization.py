#import initialization as init
#import scoring as sc
#import mutations as mt
#import crossover as cso
#from crossover import GAP, clean_alignment, print_msa
#import utils


from Bio import SeqIO
import random as r


#-------------------------------------------- FASTA FILE -------------------------------------------

def read_fasta(filename): 
    """Reading the FASTA file. Returns a list with all sequences."""

    # SeqIO.parse lê qualquer extensão (.fasta, .txt, .fa) diretamente
    records = list(SeqIO.parse(filename, "fasta"))

    #Biopython stores sequences as Seq objects. When converting to str
    sequences = [str(record.seq) for record in records]

    return sequences 


#------------------------------ RANDOM MSA -----------------------------------------------------------

def random_msa(sequences: list[str], max_offset: int, rng: r.Random) -> list[str]:
  """SKELETON- Create the initial population of MSAs (Multiple Sequence Alignments) (SAGA approach).

    Parameters
    - sequences: list[str]. Unaligned sequences (no gaps).
    - max_offset: int. Maximum number of leading gaps per sequence.
    - rng: random.Random. Random number generator (seeded, for reproducibility).

    Returns: a list[str] of random MSAs from the sequences present in the FASTA file. 
    - Invalid sequences do not occur (eg. columns containing only gaps).
   
    """
    
  aligned_sequences = []

  # 1. Add a random number of gaps (between 0 and 'max_offset') at the beginning of each sequence.
  for sequence in sequences:

    offset = rng.choice(range(max_offset + 1))   #Picks a random integer between 0 and the defined maximum value (max_offset)
    new_sequence = "-" * offset + sequence       # Joins all the hyphens (-) with the original sequence.

    aligned_sequences.append(new_sequence)


    # 2.2. Finds the longest sequence and adds gaps to the end of the others until they all match its length: 
    #      - max() : finds the maximum length
    #      - len() : length of each sequence

  max_length = max(len(sequence) for sequence in aligned_sequences)
  padded_sequences = []

  #Subtracts the sequence length from the maximum length to find the number of gaps needed
  for sequence in aligned_sequences:
    gaps_needed = max_length - len(sequence)    
    new_sequence = sequence + "-" * gaps_needed
    padded_sequences.append(new_sequence)
    
  #2.3. Remove gap-only columns with 'clean_alignment'.
  aligned_sequences = remove_gap_columns(padded_sequences)

  return aligned_sequences

  #Create the population: Repeat the creation of MSAs until the desired population size is reached."""



#--------------------- 4. POPULATION INITIALIZATION -------------------------------------------------

def initialize_population(sequences: list[str], population_size: int, max_offset: int, seed: int) -> list[list[str]]:
    """
    SKELETON - Generate an initial population of random MSAs. 

    PARAMETERS:
    - sequences: list[str]. Unaligned sequences (no gaps).
    - population_size: int. Number of MSAs in the population.
    - max_offset: int. Maximum number of leading gaps per sequence.
    - seed: int. Seed for the random generator.

    RETURN: list[list[str]]. Population (list of MSAs).

  
    1. Create a `random.Random(seed)` generator. """
    rng = r.Random(seed)
    population = []

    for i in range(population_size):
      """ 2. Call `random_msa` `population_size` times with that generator."""

      msa = random_msa(sequences, max_offset, rng)

      population.append(msa)

    return population

