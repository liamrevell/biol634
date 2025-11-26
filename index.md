---
title: Biol 634
---

<link rel="stylesheet" href="assets/style.css">

# Biol 634: Methods in Phylogenetic & Macroevolutionary Analysis (using R)

- **Professor:** Dr. Liam J. Revell ([liam.revell@umb.edu](mailto:liam.revell@umb.edu)).
- **Meeting time & place:** Wednesdays & Fridays from 1-2:15pm, Sep. 3 to Dec. 12, in Wheatley Hall W03-0022.
- **Zoom meeting room (for online participants):** [https://umassboston.zoom.us/j/93343863686](https://umassboston.zoom.us/j/93343863686) (classes will be recorded, but synchronous participation is strongly encouraged).
- **Office hours:** 1-2:30pm Mondays, 2:30-4pm Wednesdays, drop-in, & by appointment.

---

## Recommended Texts (no books are required)
1. Harmon, L. J. 2018. *Phylogenetic Comparative Methods: Learning from Trees*. ISBN: 978-1719584463. Free PDF (or buy a hard-copy on Amazon.com).  
2. Revell, L. J. and L. J. Harmon. 2022. *Phylogenetic Comparative Methods in R*. Princeton University Press, Princeton NJ. (Recommended.)  
3. Yang, Z. 2015. *Molecular Evolution: A Statistical Approach*. Oxford University Press, Oxford UK. (Recommended.)  

---

## Course Description
Evolutionary biology is unique in the biological sciences because it is, to a large extent, a historical discipline. The history of life on earth unfolds over thousands to millions of years – a timescale not particularly amenable to experimental manipulation or direct study. In the absence of paleontological data, which are sparse for many taxa, one of the best resources we have for studying evolutionary history over very long timescales is phylogenetic biology.  

In the first part of this course, students will survey the theory and application of modern phylogeny inferences. Modern phylogenetic biology, however, consists of much more than the mere estimation of evolutionary trees. Phylogenies are now used to investigate an enormous range of evolutionary questions. For instance, biologists have used phylogenies to ask questions about the tempo of species diversification or extinction; about the evolutionary processes that govern phenotypic changes in organisms through time; and about the characteristics of the long extinct ancestors of living species.  

In the second part of this course, students will learn the important theory and practical methods of phylogenetic comparative biology – the use of inferred trees in evolutionary inferences.

---

## Grading
- **Participation in discussions & exercises:** 40%  
- **In-class & take-home challenge problems:** 30%  
- **Final project presentation:** 30%  

---

## Presentations & Class Recordings

1. Sept. 3: Introductions & basic phylogeny terms & concepts. 
    - [video](https://umassboston.zoom.us/rec/share/eWXtqev8SGIrH9vhfkOGcpizbpXNhYRl_Me5r40LNo_6BugzhWFcV5jokEuylEJ7.ZPKcZxE4-ZfnX_KT?startTime=1756918436000) (Passcode required.)
    - [PDF](lec/1.IntroPhylogenies.pdf)

2. Sept. 5: Working with phylogenies in R.
    - [video](https://umassboston.zoom.us/rec/share/WRjayCcgdUd9OyGuCXP7RDtGlM1oYPtlwkKQVy8tMIpWxpv_wjQYK6WkVc7y4ntW.ssijcQ_trArQonAz) (Passcode required.)
    - [R code](ex/1/intro-R-phylogenetics.R)
    - [primates.phy](data/primates.phy)

3. Sept. 10: Introduction to phylogenetic tree inference.
    - [video](https://umassboston.zoom.us/rec/share/_E6mm6ln97vG_1p1_Uwns8bQW1fX47d1WYI8Qec7e0KD1HoFr48V0KTVUMVHm_99.xoQP4-zi2r_L6mDb) (Passcode required.)
    - [PDF](lec/2.IntroPhylogenyInference.pdf)
    - [R code](ex/2/inference-via-parsimony.R)
    - [Laurasiatherian.nex](data/Laurasiatherian.nex)
    
4. Sept. 12: Introduction to ML & model-based phylogeny inference.
    - [video](https://umassboston.zoom.us/rec/share/FGsHHVVp_Z4EL8icYeAII2ZAV6Wq2pxnRknyGxy578D24WBngK315bTcr_D2aPIA.LZgLlVATQ2TI-UIi) (Passcode required.)
    - [PDF](lec/3.MaximumLikelihood.pdf)
    - [Homework 1: Maximum Parsimony phylogeny estimation](hw/1/hw1-problem.md)
    
5. Sept. 17: Review of MP homework & more on Maximum Likelihood
    - [video 1](https://umassboston.zoom.us/rec/share/N6p5Z0IQK2uUIGCLBFhqJ-BGx4KdYw5n_BXClUAi-B6fAETDGNgMLGmfXb_iNLxG.fGRCL5Q-GPj1l3ts), [video 2](https://umassboston.zoom.us/rec/share/tScSghOmY6b7Y6zFn97Xg4tTJGJWJH_0eemRO9w5SVKg5u3aNfSQPDcK7uLDTdDw.D_azWz8uTF-cZBFD) (Passcode required.)
    - [R code](ex/3/parsimony-and-maximum-likelihood.R)
    - [Jackman-etal.nex](data/Jackman-etal.nex)
    - [Jackman-etal-tree.jpg](data/Jackman-etal-tree.jpg)

6. Sept. 19: Maximum Likelihood tree inference in R
    - [video](https://umassboston.zoom.us/rec/share/FBIRzoG7XS_roFcihQZcuZLuHfeuWs9HH9V1_bFsaj2kKoSgPzgpROjSLvTeIyam.jDNsz96P8vt9UOee) (Passcode required.)
    - [R code](ex/4/maximum-likelihood-continued.R)
    - [Google Sheet](https://docs.google.com/spreadsheets/d/1Hu0PmyClE75FrsSXHQ4XfRRMTh5G0A5HHaQW3UZ3dik/edit?usp=sharing) (Passcode required.)

7. Sept. 24: Morphological phylogeny inference (using MP or ML) and the nonparametric bootstrap
    - [video](https://umassboston.zoom.us/rec/share/V5POIJkdZZLzS35PBjtPGUXtCBOBT-5zIPU7LwT-Up8Fpu4G6BLcofF6qA1ob25k.RzDaOSP4vgxH7kWZ) (Passcode required.)
    - [craniata_tree.nwk](data/craniata_tree.nwk)
    - [R code](ex/5/morphological-phylogeny-inference.R)
    - [PDF](lec/4.Bootstrap-other-methods.pdf)

8. Sept. 26: Computing bootstrap proportions from split frequencies, consensus tree building, testing trees, and an explanation of likelihood
    - [video](https://umassboston.zoom.us/rec/share/IcEaAGyE_sG8HDRS5o-Lunh1D6KD6aVYxhVCQiPp_P63x8qohcXATZkvWkpZbnBh.8b4HhqPV4mNxKnim) (Passcode required.)
    - [PDF](lec/5.Consensus-trees-and-networks.pdf)
    - [Homework 2: Testing trees](hw/2/hw2-problem.md)
    - [concat.goby_75pct_JE2.nex.fas](data/concat.goby_75pct_JE2.nex.fas)

9. Oct. 1: Tree comparison, bootstraps, rooting, and tree distances
    - [video](https://umassboston.zoom.us/rec/share/nh0FM921dSkdu9yhip1DbdaS8KjcfFaT8eHKrQbXHG43y-fvFQbcb2myoutc23kT.1bIKjUf42l6jLkY6) (Passcode required.)
    - [Goby_MP.tre](data/Goby_MP.tre)
    - [Goby_ML.tre](data/Goby_ML.tre)
    - [R code](ex/6/tree-comparison-and-distance.R)
    
10. Oct. 3: Bayes theorem & Bayesian phylogenetic inference
    - [video](https://umassboston.zoom.us/rec/share/atIieAcNI9vWQryWIwu2eB0jnpfD9w5rFR2gz6wwcG5_ZRH3WupLtpmofL2OKWyE.nRKDn9WDkjPVXe7J) (Passcode required.)
    - [PDF](lec/6.Other-topics-inference.pdf)

11. Oct. 8: A bit about Bayesian MCMC & Distance matrix methods
    - [video](https://umassboston.zoom.us/rec/share/83yVuWYgmPW-Jz-dmK1lX7HzDfI0_MBXGSNBGXDeQBpbkHeTpVh2j4dSQksZcKKL.i32gErR3YonKhjF9) (Passcode required.)
    - [carnivore.distances.csv](data/carnivore.distances.csv)
    - [primates.dna](data/primates.dna)
    - [R code](ex/7/distance-matrix-methods.R)
    - [Homework 3: Distance matrix methods](hw/3/hw3-problem.md)
    
12. Oct. 10: Bootstrapping with neighbor-joining & The gene-tree, species-tree problem
    - [video](https://umassboston.zoom.us/rec/share/habxEf0iBSrrwRbMFrEQBGldUIslI8eAVmeYlhtxx6p0fxnz64S5-WotxVvD6U-i.VhyGRghIutJGivX-) (Passcode required.)
    - [R code](ex/8/nj-bootstrapping-analysis.R)

13. Oct. 15: The gene-tree, species-tree problem & Obtaining a chronogram (time-tree) using rate-smoothing
    - [video](https://umassboston.zoom.us/rec/share/cYrQbjEaNpOtjolp9IRSF672YtZSOuiVbhFHMHsU0c1fI6XTwRSkXKsOFo13I1FZ.588nmPc0rlk5r-Ps) (Passcode required.)
    - [LaurasiatherianML.nex](data/LaurasiatherianML.nex)
    - [R code](ex/9/time-calibration-with-rate-smoothing.R)

14. Oct. 17: Introduction to phylogenetic comparative methods (guest lecture)
    - [video](https://umassboston.zoom.us/rec/share/uXgdWdORjMVwpLMoua0ohrPUgKv4Axr-xDFicSBKQhsA3H3w-uxFR3cDJdw0jb-G.tdjLnfTehOzinXUF) (Passcode required.)
    - [PDF](lec/7.IntroComparativeMethods.pdf)
    
15. Oct. 22: Discrete phenotypic trait evolution (the M*k* model)
    - [video](https://umassboston.zoom.us/rec/share/yZ3Ap0tJcW5B2zcgDChnxLJRlppWXG1XpBarUPuJWi9J7omUUg8G9A90r2w6L5bR.reGyU9K-Qiv55PdR) (Passcode required.)
    - [PDF](lec/8.Disc-char-models.pdf)
    - [R code](ex/10/basic-mk-models.R)
    
16. Oct. 24: Custom discrete character M*k* models in R
    - [video](https://umassboston.zoom.us/rec/share/6-db14z92mq5-iyoAJqn78NTfc0qCOAPaHd0oiVABLjCFPagta5OZychUARx4Zvf.RNbTVcV3JDTXsZ0I) (Passcode required.)
    - [R code](ex/11/custom-mk-models.R)
    - [brandley_table.csv](data/brandley_table.csv)
    - [squamate.tre](data/squamate.tre)
    - [Homework 4: Custom M*k* models](hw/4/hw4-problem.md)
    
17. Oct. 29: More about M*k* model fitting in R
    - [video](https://umassboston.zoom.us/rec/share/CkbXl28S_hEXwCf-dxZw2igL8T0_23PPxADSTEyTDDrXE9mr5guetMQdvxw1Z1Jp.g244YnLuVTvK7ae-) (Passcode required.)
    - [R code](ex/12/more-about-mk-models.R)
    
18. Oct. 31: Brownian motion and phylogenetically independent contrasts
    - [video](https://umassboston.zoom.us/rec/share/KUWsFnTRGmaIV6osFJuUHGGo5McpIzMFmQ7sTa-8o_iKV5ZsnjCRVEEfWt3VVTzX.SQw5k-VtENhLy0Xd) (Passcode required.)
    - [PDF](lec/9.BrownianMotion-and-PICs.pdf)
    
19. Nov. 5: Phylogenetically independent contrasts & PGLS
    - [video](https://umassboston.zoom.us/rec/share/vlV2nHSo4v49UldPgLQNV1Oxd4jxzjfvDKZvS27A2Zk01OwkukQjMKeSvvBftb5r.bFTnY5LiCj8HGdL6) (Passcode required.)
    - [R code](ex/13/contrasts-and-pgls.R)
    
20. Nov. 7: Ancestral state estimation for discrete characters
    - [video](https://umassboston.zoom.us/rec/share/bxp24W4iwlbMDoqXCExE9NVkHRx05yRSeTQT87vH8oB3HkCC5c_0gHNdd_sEk0hB.yU55btuB-FxhXMQj) (Passcode required.)
    - [PDF](lec/10.Ancestral-States.pdf)
    
21. Nov. 12: Ancestral state reconstruction of discrete traits in R
    - [video](https://umassboston.zoom.us/rec/share/dKyWY6K7-Qr421SxNmKMZWiPZErv08NAWgcBtPeKcmFk47bmMeRZSOqwwfwVnJrK.4pnAOq1n2HM9VOs2) (Passcode required.)
    - [R code](ex/14/asr-discrete-characters.R)
    
22. Nov. 14: Ancestral state reconstruction of continuous characters
    - [video](https://umassboston.zoom.us/rec/share/gBghumgxTFdaPUu2Efs07v3OalPDwmw0emi3qFxvrqtyZxi2ZyUjVd9VXUHMZ1aR.UglYBKKsakCNc4Mw) (Passcode required.)
    - [R code](ex/15/asr-continuous-traits.R)
    
23. Nov. 19: Diversification
    - [video](https://umassboston.zoom.us/rec/share/u2shlVdTRrOZFon-YRDP21ZX_Om3XhL_nrlv_E8SltkIGf02MtDXVfqjZ0oyXuEc.M1h65MCnFOVcosK_) (Passcode required.)
    - [PDF](lec/11.Diversification.pdf)
    - [R code](ex/16/diversification-models.R)
    - [elapidae.nex](data/elapidae.nex)
    - [10kTrees_Primates.nex](data/10kTrees_Primates.nex)

24. Nov. 21: More diversification and SSE models
    - [video](https://umassboston.zoom.us/rec/share/y0Ia-Kaw9nxjkXq6A9JL4CfnKUP7-X7YgWPmJ3ZAcS6UsvQ6_-dSxZVER6AaSwlI.tyg44PZwdij12bhd) (Passcode required.)
    - [R code](ex/17/more-diversification-and-sse-models.R)
    - [grunts.csv](data/grunts.csv)
    - [grunts.phy](data/grunts.phy)
    - [pipefishes.phy](data/pipefishes.phy)
    
25. Nov. 26: Complex models of discrete & continuous trait evolution
    - [bonyfish.csv](data/bonyfish.csv)
    - [bonyfish.tre](data/bonyfish.tre)
    - [Liolaemidae.data.csv](data/Liolaemidae.data.csv)
    - [Liolaemidae.MCC.nex](data/Liolaemidae.MCC.nex)
    - [pagel94-example.csv](data/pagel94-example.csv)
    - [pagel94-example.tre](data/pagel94-example.tre)
    
Other data files:
  - [conraua_phylogeny.tre](data/conraua_phylogeny.tre)
  - [conraua_svl.csv](data/conraua_svl.csv)
  - [H3N2.csv](data/H3N2.csv)
  - [H3N2.fasta](data/H3N2.fasta)
  - [Jackman_trees.phy](data/Jackman_trees.phy)
  - [primates.csv](data/primates.csv)
  - [salamanders.nex](data/salamanders.nex)
  - [terapontidae.csv](data/terapontidae.csv)
  - [terapontidae.phy](data/terapontidae.phy)
  - [vertebrate-morphology.csv](data/vertebrate-morphology.csv)
  - [woodmouse.fas](data/woodmouse.fas)

---

## Topics & Schedule (subject to change)
We will cover only a subset of the topics & exercises given below, depending on class participant interests and prior experience.

1. Introduction to phylogenetic terms & concepts (lecture)  
2. Exercise 1: A brief introduction to phylogenetics in R (computer lab)  
3. Introduction to inferring phylogenies: the method of Maximum Parsimony (lecture)  
4. Exercise 2: Introduction to phylogeny inference using Maximum Parsimony (computer lab)  
5. Challenge problem 1: Estimating a tree from DNA sequences using MP (homework)  
6. Statistical properties of Parsimony (and comparison to likelihood) (lecture)  
7. Introduction to likelihood (and models of molecular evolution) (lecture)  
8. Exercise 3: Phylogeny inference using maximum likelihood (computer lab)  
9. Distance methods for phylogeny inference (lecture)  
10. Exercise 4: Distance methods (computer lab)  
11. The nonparametric bootstrap & other methods (lecture)  
12. Exercise 5: The non-parametric bootstrap for phylogenetic data (computer lab)  
13. Challenge problem 2: Testing trees (homework)  
14. Consensus trees & networks (lecture)  
15. Exercise 6: Consensus trees & networks (computer lab)  
16. Introduction to Bayesian phylogenetic inference (lecture)  
17. Exercise 7: Bayesian phylogeny inference (computer lab)  
18. Distances between trees (lecture)  
19. The gene-tree species-tree problem (lecture)  
20. Phylogenetically independent contrasts and PGLS (lecture)  
21. Exercise 8: Phylogenetic independent contrasts in R (computer lab)  
22. Exercise 9: Phylogenetic generalized least squares (PGLS) (computer lab)  
23. Modeling continuous character evolution on a phylogeny (lecture)  
24. Exercise 10: Fitting models of continuous character evolution (computer lab)  
25. Challenge problem 3: Fitting models of continuous character evolution in R (homework)  
26. Multi-rate, multi-regime, and multivariate models for continuous traits (lecture)  
27. Exercise 11: Multi-rate, multi-regime, and multivariate models (computer lab)  
28. Introduction to modeling discrete character evolution on a phylogeny (lecture)  
29. Exercise 12: Fitting discrete character evolution models to phylogenetic data in R (computer lab)  
30. Challenge problem 4: Discrete character evolution (homework)  
31. Other models of discrete character evolution (lecture)  
32. Exercise 13: Pagel’s model for correlated binary character evolution (computer lab)  
33. Exercise 14: Multi-regime, polymorphic trait, and hidden-rate models for discrete character evolution in R (computer lab)  
34. Overview of ancestral state reconstruction methods for continuous and discrete characters (lecture)  
35. Exercise 15: Ancestral character estimation in R (computer lab)  
36. Challenge problem 5: Stochastic character mapping (homework)  
37. Introduction to the analysis of diversification with phylogenies (lecture)  
38. Exercise 16: Studying diversification on phylogenies using R (computer lab)  
39. Challenge problem 6: Estimating diversification rates (homework)  
40. Complex models of diversification (lecture)  
41. Exercise 17: Time and density dependent diversification in R (computer lab)  
42. Exercise 18: Analyzing character-dependent diversification (computer lab)  
43. Biogeography and phylogenetic community ecology (lecture)  
44. Exercise 19: Ancestral area reconstruction on phylogenies in R (computer lab)  
45. Exercise 20: Phylogenetic community ecology (computer lab)  
46. Exercise 21: Plotting phylogenies and comparative data in R (computer lab)  
47. Challenge problem 7: Visualization (homework)  
48. Final presentations  

---

## Accommodations
UMass Boston is committed to creating learning environments that are inclusive and accessible. If you have a personal circumstance that will impact your learning and performance in this class, please let me know as soon as possible, so we can discuss the best ways to meet your needs and the requirements of the course. If you have a documented disability, or would like guidance about navigating support services, contact the Ross Center for Disability Services by email ([ross.center@umb.edu](mailto:ross.center@umb.edu)), phone (617-287-7430), or in person (Campus Center, UL Room 211). To receive accommodations, students must be registered with the Ross Center and must request accommodations each semester that they are in attendance at UMass Boston. For more information visit: [https://www.umb.edu/academics/seas/disability-services/](https://www.umb.edu/academics/seas/disability-services/). Please note that the Ross Center will provide a letter for your instructor with information about your accommodation only and not about your specific disability.

---

## Academic Integrity and Student Code of Conduct
Education at UMass Boston is sustained by academic integrity. Academic integrity requires that all members of the campus community are honest, trustworthy, responsible, respectful, and fair in academic work at the university. As part of being educated here, students learn, exercise, increase, and uphold academic integrity. Academic integrity is essential within all classrooms, in the many spaces where academic work is carried out by all members of the UMass Boston community, and in our local and global communities where the value of this education fulfills its role as a public good. Students are expected to adhere to the Student Code of Conduct, including policies about academic integrity, delineated in the University of Massachusetts Boston Graduate Studies Bulletin, Undergraduate Catalog, and relevant program student handbook(s), linked at www.umb.edu/academics/academic_integrity.

---

## Use of Artificial Intelligence (AI) Tools
AI is allowed with attribution: Use of AI tools, including ChatGPT, is permitted in this course on certain assignments, as outlined in class. To adhere to our scholarly values and to the Student Code of Conduct, students must cite any AI-generated material that informed their work; citations should include not 7 only in-text citations and listing in the references, but also the full text of cited ChatGPT (or other Large Language Model (LLM) generator) as an appendix to the assignment. Using an AI tool to generate content without proper attribution qualifies as academic dishonesty. Students are also responsible for making sure that any AI generated text does not contain false or erroneous information. If students are unsure about whether or not a source is appropriate to use in the assignment, they should contact the instructor.

---

## Health, Wellbeing, and Success
- UMass Boston is a vibrant, multi-cultural, and inclusive institution committed to ensuring that all members of our diverse campus community are able to thrive and succeed. The university provides a wide variety of resources to support students’ overall success. As we continue to deal with the ongoing impacts of the COVID-19 pandemic, these resources are more important than ever. 
- Are you in emotional distress? Call 617.287.5690 to speak with a licensed clinician 24/7 who can offer support, crisis recommendations, and assistance with finding resources.
- Have a campus question or issue? Use Here4U in the UMass Boston app or via umb.edu/here4U.
- Need help with food insecurity, legal consultation, financial counseling or emergency supplies? Contact U-ACCESS ([https://www.umb.edu/campus-life/current-students/u-access/](https://www.umb.edu/campus-life/current-students/u-access/))
- Want advice in navigating a university or life situation? Contact the Dean of Students Office at [https://www.umb.edu/deanofstudents](https://www.umb.edu/deanofstudents).
- Looking for identity-based community support? Find more resources at [https://www.umb.edu/all-of-us/](https://www.umb.edu/all-of-us/)
- Want to make the most of your academic experience? Visit [https://www.umb.edu/academics/seas/](https://www.umb.edu/academics/seas/).
- Unable to attend class on a specific date or participate in an exam or class requirement due to a religious observance? Fill out the excused absence form (requires 2-weeks’ notice) to request religious accommodation at [https://www.umb.edu/campus-life/current-students/policies/right-to-excused-absence-because-religious-belief/](https://www.umb.edu/campus-life/current-students/policies/right-to-excused-absence-because-religious-belief/).

---
