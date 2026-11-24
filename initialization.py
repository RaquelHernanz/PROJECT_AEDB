import random as r
import initialization as init
import scoring as sc
import mutations as mt
import crossover as cso
from crossover import GAP, clean_alignment, print_msa
import utils

### 4. POPULATION INITIALIZATION

def random_msa(sequences: list[str], max_offset: int, rng: r.Random) -> list[str]:
    """
    SKELETON - Create a random MSA from unaligned sequences (SAGA
    approach).

    To implement
    - Add a random number of gaps (between 0 and `max_offset`) at the
      beginning of each sequence.
    - Pad the end of each sequence with gaps until all of them have
      the length of the longest one.
    - Remove gap-only columns with `clean_alignment`.

    Parameters
    sequences: list[str]. Unaligned sequences (no gaps).
    max_offset: int. Maximum number of leading gaps per sequence.
    rng: random.Random. Random number generator (seeded, for reproducibility).

    Returns
    list[str]. A valid MSA of `sequences`.
    """
    raise NotImplementedError


def initialize_population(sequences: list[str], population_size: int,
                          max_offset: int, seed: int) -> list[list[str]]:
    """
    SKELETON - Generate an initial population of random MSAs.

    To implement
    - Create a `random.Random(seed)` generator.
    - Call `random_msa` `population_size` times with that generator.

    Parameters
    sequences: list[str]. Unaligned sequences (no gaps).
    population_size: int. Number of MSAs in the population.
    max_offset: int. Maximum number of leading gaps per sequence.
    seed: int. Seed for the random generator.

    Returns
    list[list[str]]. Population (list of MSAs).
    """
    raise NotImplementedError
