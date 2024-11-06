# dna toolkit 
import collections
Nucleotides = ["A","T","G","C"]
reverse_complement = {"A":"T", "T":"A", "G":"C", "C":"G"}

"""check the sequence to make sure it is a DNA stirng i.e weather it is a valid DNA string or not"""
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq 



def count_nuc_frequency(seq):
    # tmpfreqdict = {"A":0,"T": 0,"G":0,"C":0}
    # for nuc in seq:
    #     if nuc in tmpfreqdict:
    #         tmpfreqdict[nuc] += 1
    # return tmpfreqdict
    return dict(collections.Counter(seq))      # to avoid  all above line we used dict in collection as .counter  a one line code 


def transcription(seq):
    return seq.replace("T","U")


def reverse_complement(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement_dict[nuc] for nuc in seq])[::-1]


def gc_content(seq):
  """Calculates the GC content of a DNA sequence.

  Args:
    seq: A DNA sequence.

  Returns:
    The GC content as a percentage.
  """

  gc_count = 0
  total_count = len(seq)

  for nuc in seq:
    if nuc in ['G', 'C']:
      gc_count += 1

  gc_content = (gc_count / total_count) * 100
  return gc_content

def gc_content_subsec(seq, k=20):
    """Calculates the GC content of a sequence in subsections of size k.

    Args:
        seq: The input sequence.
        k: The size of each subsection.

    Returns:
        A list of GC content values for each subsection.
    """

    results = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        gc_content_value = gc_content(subseq)
        results.append(gc_content_value)

    return results

def translate_dna(dna_sequence):
  """Translates a DNA sequence into a protein sequence.

  Args:
    dna_sequence: The DNA sequence to translate.

  Returns:
    The corresponding protein sequence.
  """

  # Genetic code dictionary
  DNA_codons = { # here _ indicates stop codon 
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AGC': 'S', 'AGT': 'S',
    'AGA': 'R', 'AGG': 'R', 'AGC': 'S', 'AGT': 'S',
    'AAA': 'K', 'AAG': 'K', 'AAC': 'N', 'AAT': 'N',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'TAA': '_', 'TAG': '_', 'TGA': '_',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
    'TTT': 'F', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTG': 'L', 'TTT': 'F',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CCG': 'P', 'CCT': 'P', 'CCA': 'P', 'CCC': 'P',
    'CGT': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GCG': 'A', 'GCT': 'A', 'GCA': 'A', 'GCC': 'A',
    'GGG': 'G', 'GGT': 'G', 'GGC': 'G', 'GGG': 'G'
  }

  protein_sequence = ""
  for i in range(0, len(dna_sequence), 3):
    codon = dna_sequence[i:i+3]
    if codon in DNA_codons:
      amino_acid = DNA_codons[codon]
      if amino_acid == "Stop":
        break
      protein_sequence += amino_acid

  return protein_sequence


def codon_frequency(dna_sequence):
  """Calculates the frequency of each codon in a DNA sequence.

  Args:
    dna_sequence: The DNA sequence to analyze.

  Returns:
    A dictionary mapping codons to their frequencies.
  """

  codon_dict = {}
  for i in range(0, len(dna_sequence), 3):
    codon = dna_sequence[i:i+3]
    if codon in codon_dict:
      codon_dict[codon] += 1
    else:
      codon_dict[codon] = 1
  return codon_dict




import random

def generate_protein(length):
    """Generates a random amino acid sequence.

    Args:
        length: The desired length of the protein sequence.

    Returns:
        A string representing the generated protein sequence.
    """

    amino_acids = "ARNDCEQGHILKMFPSTWYV"  # All possible amino acids
    protein_sequence = ""

    for _ in range(length):
        random_index = random.randint(0, len(amino_acids) - 1)
        protein_sequence += amino_acids[random_index]

    return protein_sequence




