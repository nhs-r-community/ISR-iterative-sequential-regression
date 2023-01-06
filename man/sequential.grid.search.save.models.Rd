\name{sequential.grid.search.save.models}
\alias{sequential.grid.search.save.models}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
nnn
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sequential.grid.search.save.models(dataset, model.type, resp.variable, expl.variables.without.segm.variable, segm.variable, model.parameters, angle.param = 0, start.point, end.point, min.distance, interval.length, time.unit, criterion, criterion.difference)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
}
  \item{model.type}{
%%     ~~Describe \code{model.type} here~~
}
  \item{resp.variable}{
%%     ~~Describe \code{resp.variable} here~~
}
  \item{expl.variables.without.segm.variable}{
%%     ~~Describe \code{expl.variables.without.segm.variable} here~~
}
  \item{segm.variable}{
%%     ~~Describe \code{segm.variable} here~~
}
  \item{model.parameters}{
%%     ~~Describe \code{model.parameters} here~~
}
  \item{angle.param}{
%%     ~~Describe \code{angle.param} here~~
}
  \item{start.point}{
%%     ~~Describe \code{start.point} here~~
}
  \item{end.point}{
%%     ~~Describe \code{end.point} here~~
}
  \item{min.distance}{
%%     ~~Describe \code{min.distance} here~~
}
  \item{interval.length}{
%%     ~~Describe \code{interval.length} here~~
}
  \item{time.unit}{
%%     ~~Describe \code{time.unit} here~~
}
  \item{criterion}{
%%     ~~Describe \code{criterion} here~~
}
  \item{criterion.difference}{
%%     ~~Describe \code{criterion.difference} here~~
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
function (dataset, model.type, resp.variable, expl.variables.without.segm.variable, 
    segm.variable, model.parameters, angle.param = 0, start.point, 
    end.point, min.distance, interval.length, time.unit, criterion, 
    criterion.difference) 
{
    models <- list()
    fixed.breakpoints <- numeric()
    last.end.point <- 1
    found.new.breakpoint <- TRUE
    while (found.new.breakpoint) {
        alpha <- sequential.search.for.new.breakpoint.whole.interval(dataset = dataset, 
            model.type = model.type, resp.variable = resp.variable, 
            expl.variables.without.segm.variable = expl.variables.without.segm.variable, 
            segm.variable = segm.variable, model.parameters = model.parameters, 
            angle.param = angle.param, fixed.breakpoints = fixed.breakpoints, 
            last.end.point = last.end.point, start.point = start.point, 
            end.point = end.point, min.distance = min.distance, 
            interval.length = interval.length, time.unit = time.unit, 
            criterion = criterion, criterion.difference = criterion.difference)
        found.new.breakpoint <- alpha$found.new.breakpoint
        fixed.breakpoints <- alpha$new.fixed.breakpoints
        last.end.point <- alpha$new.end.point
        models <- append(models, list(list(best.model = alpha$new.best.model, 
            good.breakpoint.models = alpha$good.breakpoint.models, 
            start.point = start.point, end.point = transform.breakpoint(start.point = start.point, 
                breakpoint = alpha$new.end.point, interval.length = interval.length, 
                time.unit = time.unit)$transf.breakpoint)))
    }
    return(list(models = models))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
