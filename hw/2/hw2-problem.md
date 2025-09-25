---
title: HW 2
---

<link rel="stylesheet" href="../../assets/style.css">

## Homework problem 2: Testing trees

In this challenge you'll use a dataset for gobioid fishes ([concat.goby_75pct_JE2.nex.fas](../../data/concat.goby_75pct_JE2.nex.fas)) from Johnson et al. ([2025](https://doi.org/10.1016/j.ympev.2025.108424)).

Download the following data file from the class page: [concat.goby_75pct_JE2.nex.fas](../../data/concat.goby_75pct_JE2.nex.fas).

1. Estimate a Maximum Parsimony tree for this dataset using *phangorn*. You can use the function `pratchet` to get a good tree, but I suggest trying both with & without a random starting tree to see if you converge on the same result.
2. Estimate a Maximum Likelihood tree for this dataset using *phangorn*. You should first select an appropriate model using `modelTest`. You can use the function `pml_bb` if you want, or `pml` followed by `optim.pml`. Note that this will probably take a while!
3. Both `pratchet` and `pml_bb` return phylogenies that contain estimated bootstrap proportions. Identify where this information is stored in each object and plot both trees with their bootstrap proportions shown.
4. Statistically compare your two trees using the Shimodaira-Hasegawa method. Are the two trees statistically different by this measure?

Please note: if you elect to turn in your solution (this is not required) please do so in Rmarkdown and rendered (HTML or PDF) format.

References:

1. Data for this exercise are here: [doi:10.5061/dryad.n02v6wx96](https://doi.org/10.5061/dryad.n02v6wx96) (via Data Dryad).
2. Full reference: Johnson K, Tornabene L, Li C, Ruber L, Schliewen U, Hogan D, Pezold F (2005) Exon-capture data resolve relationships resulting from a rapid radiation within family Gobiidae. Molecular Phylogenetics and Evolution 212: 108424. [PDF](https://doi.org/10.1016/j.ympev.2025.108424).
