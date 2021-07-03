### Future Work

One conclusion of this project was that the algorithm presented is an effective classifier, which in many cases outperforms the current state of the art. The most notable downside is its speed. It has quadratic time complexity, relative to parametric approaches which are log-linear for example.

In part this is inherent to using a non-parametric cost function, which has less scope for memoization. Within the final chapter of the thesis, it was noted however that the existing steps taken to improve its efficiency could be generalised and in fact taken a step further. Specifically, alternative norms on the ECDF such as the Cramer-von-Mises which take a polynomial form have additional properties which could be exploited with the use of Fenwick [1] and Order Statistic Trees [2] to reduce the compute time to polylogarithmic-linear.

The core question therefore would be regarding the extent to which differences in classification performance justify differences in computational cost. As a matter of personal interest, this is a question to be answered in Q4 2021.

###### References:
<sub> [1] 	P. Fenwick (1994), "A new data structure for cumulative frequency tables", Software: Practice and Experience, vol 24, pp. 327 - 336 \
[2]   T. Cormen, C. Leiserson, R. Rivest, \& S. Clifford (2009), "Introduction to Algorithms", MIT Press and McGraw-Hill, ISBN 0-262-03384-4  </sub>
