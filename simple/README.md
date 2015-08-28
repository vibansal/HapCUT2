
08/11/2015 

Modified HAPCUT algorithm that is designed to handle long reads (also works for mate-pair data) and does not use the read-haplotype graph for max-cut computations. Instead, it uses the same greedy heuristic to update the haplotypes using the exact haplotype likelihoods.

Features of this code:

1. Operates on likelihood function and also models chimeric fragments
2. storage requirements for graph are significantly less for long reads 
