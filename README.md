# PyCOTS: CRISPR Off-Target Search

## To run

1. Install poetry
2. In this folder, `poetry install`
2. `poetry run pycots <args>` or `pycots <args>`

## Assumptions

- Input FASTA sequence is converted to uppercase (ignoring soft-masking)
- "N" nucleotides in the input reference sequence prevent a match covering that position
- Population SNVs are provided and can modify the reference sequence in the target sequence (spacer) and the PAM.

## Matching Algorithm

1. Load one contig at a time from the FASTA file
2. Convert this n-length contig into a 4 x n binary array, where each row indicates whether the reference sequence may be A, C, G, or T.
  - 'N' is represented with all zeros to avoid matching in padded regions
3. Load population variants to modify the matrix, allowing a position to potentially be any or all of A,C,G, or T (multiple rows may contain a value of 1 in the same column)
4. Convert the spacer pattern into a similar matrix, optionally using ambiguous IUPAC nucleotide codes
  - 'N' is represented by all ones, indicating it may match any nucleotide
5. Count the number of matching spacer pattern positions at each position in the reference genome, keeping those with at most 1 mismatch
6. Adjacent to each of these, check for matches to a PAM pattern but convert the matrix back into possible sequences rather than counting the number of matching positions

### Matching Example:

"ACGAT" converts to a matrix: 

```
10010
01000
00100
00000
```

A pattern "AN" converts to a matrix:

```
11
01
01
01
``

Performing matrix multiplication at the first position gets

```
10  
01
00
00
```

Indicating 2 matched positions, and a possible sequence of AC

Likewise at the 2nd position, the result is

```
00
00
01
00
```

Indicating only one match (a G in the second position of the pattern)