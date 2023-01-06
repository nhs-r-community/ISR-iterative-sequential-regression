\name{transf.resp.variable}
\alias{transf.resp.variable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ppp
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
transf.resp.variable(dataset, var, transformation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
}
  \item{var}{
%%     ~~Describe \code{var} here~~
}
  \item{transformation}{
%%     ~~Describe \code{transformation} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dataset, var, transformation) 
{
    output.dataset <- dataset
    transf.var <- paste("transf", var, sep = ".")
    if (transformation == "boxcox") {
        temp <- bc.transform(var = dataset[[var]])
        output.dataset[, transf.var] <- temp$transf.var
        lambda <- temp$lambda
    }
    else {
        output.dataset[, transf.var] <- dataset[, var]
        lambda <- "na"
    }
    return(list(output.dataset = output.dataset, lambda = lambda, 
        transf.var = transf.var))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
