#!/usr/bin/env python3

import random
from sys import argv

# ----------------------------
# Parameters
# ----------------------------

SEQ_LENGTH = 100000
NUM_SEQUENCES = 20
EDIT_FRACTION = 0.30

DNA = "ACGT"

# Edit probabilities
P_SUB = 0.10
P_INS = 0.45
P_DEL = 0.45

random.seed()


def random_dna(length):
    """Generate a random DNA sequence."""
    return "".join(random.choices(DNA, k=length))


def mutate_sequence(sequence, edit_fraction=0.10):
    """
    Apply random substitutions, insertions and deletions.

    Each edit is chosen independently, so sequence lengths
    will vary naturally.
    """

    seq = list(sequence)

    total_edits = int(len(sequence) * edit_fraction)

    for _ in range(total_edits):

        r = random.random()

        # ------------------------
        # Substitution
        # ------------------------
        if r < P_SUB:

            if len(seq) == 0:
                continue

            pos = random.randrange(len(seq))
            old = seq[pos]
            seq[pos] = random.choice([b for b in DNA if b != old])

        # ------------------------
        # Insertion
        # ------------------------
        elif r < P_SUB + P_INS:

            pos = random.randrange(len(seq) + 1)
            seq.insert(pos, random.choice(DNA))

        # ------------------------
        # Deletion
        # ------------------------
        else:

            if len(seq) <= 1:
                continue

            pos = random.randrange(len(seq))
            del seq[pos]

    return "".join(seq)


def write_fasta(sequences, filename):
    with open(filename, "w") as f:
        for i, seq in enumerate(sequences):
            f.write(f">seq_{i}\n")
            f.write(seq + "\n")


def main():

    if len(argv) != 2:
        print(f"Usage: {argv[0]} output.fasta")
        return

    reference = random_dna(SEQ_LENGTH)

    sequences = [reference]

    for _ in range(NUM_SEQUENCES - 1):
        sequences.append(mutate_sequence(reference, EDIT_FRACTION))

    write_fasta(sequences, argv[1])

    print("Generated sequences:\n")

    for i, seq in enumerate(sequences):
        print(f"seq_{i:2d}: length = {len(seq)}")

    print(f"\nSaved to {argv[1]}")


if __name__ == "__main__":
    main()
