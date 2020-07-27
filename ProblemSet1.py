# %% [markdown]
# # Problem Set 1
# %%
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
import numpy as np

# %% Code Setup
my_seq = Seq("AGTACGCAACTGGTA")
my_seq2 = Seq("AGTACACGCGCGCTGGTA", IUPAC.unambiguous_dna)
my_prot = Seq("AGTACACTGGTA", IUPAC.protein)

# %% [markdown]
# Complete the codes below

# %%
# my_seq's variable type, my_prot's variable type
type_of_my_seq =
type_of_my_prot =


# %%
# Join sequences my_seq and my_seq2 into my_seq3.
my_seq3 =


# %%
# Find the number of base T in my_seq
t_myseq =

# %%
# Find the GC content of a DNA given the function GC
gc_myseq =


# %%
# Find the geometric mean of base A in my_seq and my_seq2
gmean =

# %%
# Find the abundance of A, T, C and G in my_seq. Report the answer as a list of abundances.
abundance_my_seq =

# %%
# Find the abundance of A, T, C and G in my_seq. Report the answer as a numpy array of abundances.
abundance_my_seq_array =

# %%
# Find the abundance of A, T, C and G in my_seq. Report the answer as a dictionary of abundances.
abundance_my_seq_dict =

# %%
# Transcribe my_seq to mRNA using the method .reverse_complement()
transcript_my_seq =


# %%
# Flip 5' and 3' end of my_seq
flipped_my_seq =


# %% [markdown]
# # Troubleshooting exercise
# The following code convert a DNA sequence to its protein sequence and displays the consensus sequence. Can you identify the errors and correct it?

# %%
# Import additional libraries
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import logomaker
import matplotlib.pyplot as plt
import pandas as pd

# %% [markdown]
# Perform a blast search of accession number [AY527219](https://www.ncbi.nlm.nih.gov/nuccore/AY527219). Saves the blast results in `blast_results.xml`

# %%
#result_handle = NCBIWWW.qblast("blastn", "nt", "AY527219")
#with open("blast_results.xml", "w") as out_handle:
#    out_handle.write(result_handle.read())

#result_handle.close()

# %% [markdown]
# Load blast results. 
#
# **Spot the error!**

# %%
result_handle = open(blast_results.xml)
blast_record = NCBIXML.read(result_handle)

# %% [markdown]
# The following code snippet extract the first 120 base pairs from the matching genes. 
# **Spot the error!**

# %%
# Extract first 120 bases of matching sequences.
base_pos1 = 0
base_pos2 = 120
matches = []
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        matches.append(hsp.sbjct[base_pos1:base_pos2])
        
pd.DataFrame(matches)

# %% [markdown]
# Translate the codon to their respective amino acids.

# %%
# Convert gene to protein sequences
matches_prot = []
expected_length_primary_sequence = 5 # This is incorrect, choose the right expected length?
for sequence in matches:
    pri_seq = Seq(sequence).translate() # Translate sequence
    try:
        assert len(pri_seq) == expected_length_primary_sequence
    except:
        raise AssertionError('The protein length obtained is different from the expected length. Change accordingly')
    matches_prot.append(list(pri_seq))

# Check if 
matches_df = pd.DataFrame(matches_prot)
matches_df

# %% [markdown]
# Count the number of each type of amino acid at each position.

# %%
prot_aa = []
for ind, pri_seq in matches_df.iteritems():
    aa_abundance = pri_seq.value_counts()
    prot_aa.append(aa_abundance)

prot_df = pd.DataFrame(prot_aa).fillna(0)
prot_df

# %% [markdown]
#

# %%
# Display Consensus Sequence
crp_logo = logomaker.Logo(prot_df, figsize=(10,2), color_scheme='chemistry')

# %% [markdown]
# Identify the error.

# %%
# style and show figure
crp_logo.ax.set_xlabel(Percentage)
crp_logo.ax.set_title(Primary consensus sequence)
crp_logo.fig

# %%
