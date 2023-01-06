# ISR (iterative sequential regression)

This package has been forked from original code developed for the 
[manuscript by Schlackow et al on Surveillance of Infection Severity: A Registry Study of Laboratory Diagnosed Clostridium difficile](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1001279). 

It is being reformatted into a fully open package by NHS-R Community to continue 
its support and note that issues have been raised to assist in that reformatting. 
These are open to contributions from anyone who wishes to be involved.

## Installation

```
# install.packages("remotes")
remotes::install_github("https://github.com/nhs-r-community/ISR-iterative-sequential-regression-")
```

From the second fork of this repository the package relies on the code 
`source(librariesupdatingpackagesfunctions_2020.R)` in order to update the 
package. It specifically updates 3 of the ISR functions so that they are 
compatible with 2022 versions of packages such as {ggplot2}. 

From the [README](https://github.com/karinadorisvihta/ISR-iterative-sequential-regression-):

> Note lines 62-65 have been adjusted to get annual rate ratios from 
Poisson/negative binomial models as the model uses posixCT objects it calculates 
a rate ratio per second. Adjust these according to your own needs. So perhaps 
take the exp away if not using a log-link, or multiply just by 24 if you want 
a daily rate ratio. 

> Additionally, this also adds an extra parameter to the plot and master 
functions, namely ylim.max, which allows you to set the upper limit of the 
y-axis so that if you want to create multiple plots, you can have the same 
limits for all of them. 

> Similarly, you can change these further as per your own project's needs.

This package has further been used in the following publications:
https://www.thelancet.com/pdfs/journals/laninf/PIIS1473-3099(18)30353-0.pdf 

# Contribution

This is an NHS-R Community project that is open for anyone to contribute
to in any way that they are able. The project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

If you want to learn more about this project, please join the discussion
at [the NHS-R Community Slack group](https://nhsrcommunity.slack.com/).

The simplest way to contribute is to raise an issue detailing the
feature or functionality you would like to see added, or any unexpected
behaviour or bugs you have experienced.
