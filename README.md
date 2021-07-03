### Non-Parametric Changepoint Detection

This repository stores the code for an algorithm developed for my master's thesis on the topic of offline univariate multiple changepoint detection.

In the early stages of the research phase, the non-parametric likelihood was found to be a powerful two-sample goodness of fit test. The results of the papers [1]-[2] were found to be 
reproducible in the sense that it performed favorably relative to other well-known non-parametric cost functions such as the Cramer-von-Mises and Anderson-Darling statistics. This was extended to a changepoint detection algorithm using Binary Segmentation [3], with properties of the likelihood and empirical distribution function being manipulated to improve its efficiency. 

In simulation studies it was found to improve upon the classification power and speed of an earlier algorithm [4] using a different variant of the statistic but with the PELT [5] optimisiation approach.

###### Acknowledgements:

With thanks to my supervisors Dr. Dean Bodenham and Professor Niall Adams for their patience, guidance, and for introducing me to an interesting problem.

###### References:
<sub> [1] 	J. Zhang (2002), "Powerful goodness-of-fit tests based on the likelihood ratio", Journal of the Royal Statistical Society, vol. 64, pp. 281 - 294 \
[2]   J. Einmahl \& I. Mckeague (2003), "Empirical likelihood based hypothesis testing", Bernoulli, vol. 9, pp. 267 - 290 \
[3] 	A. Scott \& M. Knott (1974), "A cluster analysis method for grouping means in the analysis of variance" , Biometrics, vol. 30, pp. 507 - 512 \
[4]   K. Haynes, P. Fearnhead, \& I. Eckley (2017), "A computationally efficient nonparametric approach for changepoint detection", Springer, Journal of Statistics and Computing, vol. 27, pp. 1293 - 1305 \
[5] 	R. Kilick, P. Fearnhead, \& I. A. Eckley (2012), "Optimal Detection of Changepoints With a Linear Computational Cost", Journal of the American Statistical Association, vol. 107, pp. 1590 - 1598 </sub>
