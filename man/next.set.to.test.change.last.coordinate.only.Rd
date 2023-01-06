\name{next.set.to.test.change.last.coordinate.only}
\alias{next.set.to.test.change.last.coordinate.only}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
kkk
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
next.set.to.test.change.last.coordinate.only(current.set, number.of.intervals, min.distance)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{current.set}{
%%     ~~Describe \code{current.set} here~~
}
  \item{number.of.intervals}{
%%     ~~Describe \code{number.of.intervals} here~~
}
  \item{min.distance}{
%%     ~~Describe \code{min.distance} here~~
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
function (current.set, number.of.intervals, min.distance) 
{
    n <- length(current.set)
    if (current.set[n] + 1 + min.distance > number.of.intervals) 
        next.set <- "max permitted"
    else {
        next.set <- current.set
        next.set[n] <- next.set[n] + 1
    }
    return(list(next.set = next.set))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
