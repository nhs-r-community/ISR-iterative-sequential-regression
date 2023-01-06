\name{add.extra.coordinate}
\alias{add.extra.coordinate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
aaa
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
add.extra.coordinate(current.set, number.of.intervals, min.distance)
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
    current.last.coordinate <- if (length(current.set) == 0) 
        0
    else current.set[n]
    next.last.coordinate <- current.last.coordinate + min.distance
    if (next.last.coordinate + min.distance > number.of.intervals) 
        next.set <- "max permitted"
    else next.set <- append(current.set, next.last.coordinate)
    return(list(next.set = next.set))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
