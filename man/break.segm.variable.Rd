\name{break.segm.variable}
\alias{break.segm.variable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
eee
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
break.segm.variable(angle.param = 0, dataset, segm.variable, breakpoints)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{angle.param}{
%%     ~~Describe \code{angle.param} here~~
}
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
}
  \item{segm.variable}{
%%     ~~Describe \code{segm.variable} here~~
}
  \item{breakpoints}{
%%     ~~Describe \code{breakpoints} here~~
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
function (angle.param = 0, dataset, segm.variable, breakpoints) 
{
    output.dataset <- dataset
    output.dataset[, paste("char", segm.variable, 0, sep = ".")] <- output.dataset[, 
        segm.variable]
    split.variables <- paste("char", segm.variable, 0, sep = ".")
    if (angle.param == 0) {
        if (length(breakpoints) > 0) 
            for (j in 1:length(breakpoints)) {
                output.dataset[, paste("char", segm.variable, 
                  j, sep = ".")] <- ifelse(output.dataset[, segm.variable] < 
                  breakpoints[j], difftime(dataset[, segm.variable], 
                  dataset[, segm.variable]), difftime(dataset[, 
                  segm.variable], breakpoints[j], units = "secs"))
                split.variables <- append(split.variables, paste("char", 
                  segm.variable, j, sep = "."))
            }
    }
    else {
        if (length(breakpoints) >= 1) {
            output.dataset[, paste("char", segm.variable, 0, 
                sep = ".")] <- pmin(output.dataset[, segm.variable], 
                breakpoints[1])
            for (j in 1:length(breakpoints)) {
                lower.bound <- breakpoints[j]
                upper.bound <- if (j < length(breakpoints)) 
                  breakpoints[j + 1]
                else max(dataset[, segm.variable])
                output.dataset[, paste("char", segm.variable, 
                  j, sep = ".")] <- difftime(pmin(pmax(output.dataset[, 
                  segm.variable], lower.bound), upper.bound), 
                  lower.bound)
                split.variables <- append(split.variables, paste("char", 
                  segm.variable, j, sep = "."))
            }
        }
    }
    return(list(output.dataset = output.dataset, split.variables = split.variables))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
