# pepAligner
Matche peptide sequences from a list to protein sequences from a FASTA file.

This script matches peptide sequences from a list to protein sequences from a FASTA file. It processes unique peptide sequences from a CSV file and searches for fuzzy matches in a FASTA file using Levenshtein distance with up to 1 mismatch (adjustable via the `max_mismatches` parameter). 
⚠️Please, observe that the time will increase as you load larger lists of peptides and/or protein FASTA files.

## USAGE:
In your python terminal type:
> python PepAligner.py [path/to/fasta_file] [path/to/peptide_file] [path/to/output_file]

## The CSV output contains 6 columns:
- **Protein_ID:** Protein identifier.
- **Peptide_Count:** Total number of peptides matching this protein.
- **Unique_Peptide_Count** shows how many peptides are unique to each protein.
- **Peptide_Sequences:** Peptide sequences with Levenshtein distance in parentheses.
- **Amino_Acid_Substitutions:** shows mutations in format like "A1G;K5R" (Ala at position 1 to Gly, Lys at position 5 to Arg).
- **Grantham_Distances:** shows the Grantham distance for each substitution. Scores range from 5 (conservative, e.g., I/L) to 215 (radical, e.g., C/W).

## Key references:
> - **Levenshtein distance:** Levenshtein, V. I. Binary Codes Capable of Correcting Deletions, Insertions and Reversals. Soviet Physics Doklady, 1966 Vol. 10, p.707
> - **Grantham distance:** Grantham R. Amino acid difference formula to help explain protein evolution. Science. 1974 Sep 6;185(4154):862-4. doi: 10.1126/science.185.4154.862. PMID: 4843792.
> - **Grantham distance stratification:** Li WH, Wu CI, Luo CC. Nonrandomness of point mutation as reflected in nucleotide substitutions in pseudogenes and its evolutionary implications. J Mol Evol. 1984;21(1):58-71. doi: 10.1007/BF02100628. PMID: 6442359.

**Grantham's distance (d) between amino acids.**
- Conservative if 0 < d ≤ 50,
- Moderately conservative if 50 < d ≤ 100,
- Moderately radical if 100 < d ≤ 150,
- Radical if d > 150.
