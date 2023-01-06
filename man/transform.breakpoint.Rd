\name{transform.breakpoint}
\alias{transform.breakpoint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
qqq
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
transform.breakpoint(start.point, breakpoint, interval.length, time.unit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start.point}{
%%     ~~Describe \code{start.point} here~~
}
  \item{breakpoint}{
%%     ~~Describe \code{breakpoint} here~~
}
  \item{interval.length}{
%%     ~~Describe \code{interval.length} here~~
}
  \item{time.unit}{
%%     ~~Describe \code{time.unit} here~~
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
function (start.point, breakpoint, interval.length, time.unit) 
{
    if (identical(breakpoint, numeric(0))) 
        transf.breakpoint <- character(0)
    else {
        distance <- breakpoint * interval.length
        transf.breakpoint <- lapply(distance, function(x) seq(from = as.Date(start.point), 
            length.out = 2, by = paste(x, time.unit))[2])
        transf.breakpoint <- as.POSIXct(as.Date(unlist(transf.breakpoint), 
            origin = "1970-01-01"))
    }
    return(list(transf.breakpoint = transf.breakpoint))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
