from DNAToolkit import *
import random
# creating a random DNA sequence for testing: 
rnd_dna_str = ''.join([random.choice(Nucleotides)
                       for nuc in range (6)])

DNA_str = validateSeq(rnd_dna_str)


print(f'\nSequence: {DNA_str}\n')
print(f'[1] + sequence length : {len(DNA_str)}\n')
print(f'[2] + Nucleotide frequency {count_nuc_frequency(DNA_str)}\n')
print(f'[3] + DNA_to_RNA transcription: {transcription(DNA_str)}\n')
print(f"[4]+ DNA_string + Reverse complement: \n5'{DNA_str} 3'")
print(f"   {''.join(['|' for c in range(len(DNA_str))])}")
print(f"3'{reverse_complement(DNA_str)[::-1]} 5'[complement]")
print(f"5'{reverse_complement(DNA_str)} 3' [rev.complement]\n")
print(f'[5] + GC content: {gc_content(DNA_str)}%\n')
print(f'[6] + GC content in subsection k=5:{gc_content_subsec(DNA_str, k=5)}\n')
print(f"[7]+ amino acid sequence form DNA: {translate_dna(DNA_str)}\n")
print(f"[8] + codon frequency (L): {codon_frequency(DNA_str)}\n")
protein_length = len(DNA_str) // 3  # Assuming 3 nucleotides per amino acid (codons)
print(f"[9] + generate proteins : {generate_protein(protein_length)}\n")