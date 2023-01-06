\name{bc.transform}
\alias{bc.transform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ccc
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bc.transform(var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var}{
%%     ~~Describe \code{var} here~~
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
function (var) 
{
    require(MASS)
    bc.transf <- boxcox(var ~ 1, lambda = seq(-5, 5, 1/10))
    lambda <- bc.transf$x[which.max(bc.transf$y)]
    transf.var <- if (lambda == 0) 
        log(var)
    else (var^lambda - 1)/lambda
    return(list(lambda = lambda, transf.var = transf.var))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
