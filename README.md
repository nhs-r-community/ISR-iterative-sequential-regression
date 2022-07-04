# ISR (iterative sequential regression)

This release contains the original code developed for the manuscript [manuscript by Schlackow et al on Surveillance of Infection Severity: A Registry Study of Laboratory Diagnosed Clostridium difficile](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1001279) as well as some extra helpful pieces. 

Download the code, and install the package cotained inside .tar.gz using the command in your R script install.packages("/.../isr_1.0.tar.gz", repos=NULL, type="source") .

In this fork we provide an R script librariesupdatingpackagesfunctions_2020.R that should be run, simplest using the R base function source() https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/source in order to update the code to work with current version of packages such as ggplot2.

This script updates 3 of the ISR functions so that they are compatible with 2022 versions of different required packages such as ggplot2. 

Note lines 62-65 have been adjusted to get annual rate ratios from Poisson/negative binomial models as the model uses posixCT objects it calculates a rate ratio per second. Adjust these according to your own needs. So perhaps take the exp away if not using a log-link, or multiply just by 24 if you want a daily rate ratio. 

Additionally, this also adds an extra parameter to the plot and master functions, namely ylim.max, which allows you to set the upper limit of the y-axis so that if you want to create multiple plots, you can have the same limits for all of them. 

Similarly, you can change these further as per your own project's needs.

