### User Guide

An implementation of the main algorithm is contained within the NPL-BS.cpp file. It has been designed with an R interface and C++ backend through use of Rcpp [1]. 
An explanation of its inner mechanics can be found in Chapter 6 of the document contained within the Thesis subfolder. A minimal working example illustrating how to call the method using the available code follows. 


```

#Necessary dependency
library(Rcpp)

#Edit to personalised location of NPL_BS.cpp file
custom_directory = 'x'

#Compile & load function into R environment
Rcpp::sourceCpp(custom_directory)

#A simple example - mean shift of 1 after 100 points with t distributed noise and 1 degree of freedom
example_data = c(rep(0, 100), rep(1, 100)) + rt(200,1)
penalty = 0.5 * log(200) #BIC penalty

#Calling NPL BS
NPL_BS(experiment_data, penalty)

```
A simple comparison to other methods would illustrates the methods upside. Typically, it precisely identifies the change where existing approaches instead return only false positives.

```

#Wild Binary Segmentation - unsurprisingly fails due to normal assumptions
library(wbs)
changepoints(wbs(example_data))

#PELT - unsurprisingly fails due to normal assumptions
library(changepoint)
cpt.mean(example_data, method = "PELT", class = FALSE)

#ED-PELT - despite being non-parametric it is also unsatisfactory
library(changepoint.np)
print(cpt.np(example_data, nquantiles = 4 * log(200), class = FALSE))

```

While less interest was taken in its development throughout the course of the project, code used to create simulation studies contained within the broader thesis can also be found in the utilities.R file for the purposes of reproducibility. 

###### References:
<sub> [1] D. Eddelbuettel \& R. Francois (2011), “Rcpp: Seamless R and C++ Integration”, Journal of Statistical Software, vol. 40, pp. 1 - 18 </sub>
