\name{find.smallest}
\alias{find.smallest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ggg
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
find.smallest(list.of.models, criterion)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list.of.models}{
%%     ~~Describe \code{list.of.models} here~~
}
  \item{criterion}{
%%     ~~Describe \code{criterion} here~~
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
function (list.of.models, criterion) 
{
    value <- list.of.models[[1]][[criterion]]
    index <- 1
    for (i in 1:length(list.of.models)) if (list.of.models[[i]][[criterion]] < 
        value) {
        index <- i
        value <- list.of.models[[i]][[criterion]]
    }
    return(list(value = value, index = index))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
