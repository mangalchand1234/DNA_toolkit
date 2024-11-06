# here is file whrere all testing is done
from DNAToolkit import *  # here indicates import every thing
import random
#your dna string
rnd_dna = "ATGCAGGCGC"
print(validateSeq(rnd_dna))


# creating a random DNA sequence for testing: 
rnd_dna_str = ''.join([random.choice(Nucleotides)
                       for nuc in range (50)])
print(validateSeq(rnd_dna_str))
# print the 2 def i.e count_nuc_freq
print(count_nuc_frequency(rnd_dna_str))
print(transcription(rnd_dna_str))
print(reverse_complement(rnd_dna_str))



