---
title: HW 3
---

<link rel="stylesheet" href="../../assets/style.css">

## Homework problem 3: Neigbor-joining and the bootstrap

In this challenge you'll use the simple nucleotide sequence alignment for primates: [primates.dna](../../data/primates.dna). You're also welcome to test out your code with other datasets that we've used in the course!

Download the following data file from the class page: [primates.dna](../../data/primates.dna).

1. Compute a distance matrix for these DNA sequence data under the Jukes-Cantor model as we did in class. (You can also try other models if you like.)
2. Calculate a neighbor-joining tree from the original, full dataset.
3. Bootstrap resample your input `"DNAbin"` data object at least 100 times. For each resampled dataset, calculate a NJ tree.
4. Evaluate bootstrap proportions on your *original* tree from this set of trees.

Please note: if you elect to turn in your solution (this is not required) please do so in Rmarkdown and rendered (HTML or PDF) format.
