\name{breakpoint.models.fixed.interval}
\alias{breakpoint.models.fixed.interval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fff
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
breakpoint.models.fixed.interval(dataset, model.type, resp.variable, expl.variables.without.segm.variable, segm.variable, model.parameters, angle.param = 0, fixed.breakpoints, start.point, end.point, min.distance, interval.length, time.unit)
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
  \item{fixed.breakpoints}{
%%     ~~Describe \code{fixed.breakpoints} here~~
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
    segm.variable, model.parameters, angle.param = 0, fixed.breakpoints, 
    start.point, end.point, min.distance, interval.length, time.unit) 
{
    number.of.intervals <- length(seq(from = start.point, to = end.point, 
        by = paste(interval.length, time.unit))) - 1
    list.of.breakpoints <- list(straight.line.breakpoints = fixed.breakpoints)
    next.breakpoints <- add.extra.coordinate(current.set = fixed.breakpoints, 
        number.of.intervals = number.of.intervals, min.distance = min.distance)$next.set
    while (!identical(next.breakpoints, "max permitted")) {
        list.of.breakpoints <- append(list.of.breakpoints, list(next.breakpoints))
        next.breakpoints <- next.set.to.test.change.last.coordinate.only(current.set = next.breakpoints, 
            number.of.intervals = number.of.intervals, min.distance = min.distance)$next.set
    }
    list.of.models <- list()
    for (b in list.of.breakpoints) {
        transf.breakpoints <- transform.breakpoint(start.point = start.point, 
            breakpoint = b, interval.length = interval.length, 
            time.unit = time.unit)$transf.breakpoint
        temp.fit <- fit.model.with.breakpoints(dataset = dataset, 
            model.type = model.type, resp.variable = resp.variable, 
            expl.variables.without.segm.variable = expl.variables.without.segm.variable, 
            segm.variable = segm.variable, breakpoints = transf.breakpoints, 
            model.parameters = model.parameters, angle.param = angle.param)
        list.of.models[[paste("breakpoints=", b, collapse = ",")]] <- list(number.of.breakpoints = length(b), 
            breakpoints = b, transf.breakpoints = transf.breakpoints, 
            segm.var.coeffs = temp.fit$segm.var.coeffs, segm.var.significance = temp.fit$segm.var.significance, 
            aic = AIC(temp.fit$model), bic = AIC(temp.fit$model, 
                k = log(nrow(dataset))), formula.to.model = temp.fit$formula.to.model)
    }
    straight.line.model <- list.of.models[[paste("breakpoints=", 
        fixed.breakpoints, collapse = ",")]]
    return(list(list.of.breakpoints = list.of.breakpoints, list.of.models = list.of.models, 
        straight.line.model = straight.line.model))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
