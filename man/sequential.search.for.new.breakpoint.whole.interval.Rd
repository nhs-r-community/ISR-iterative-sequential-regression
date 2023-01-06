\name{sequential.search.for.new.breakpoint.whole.interval}
\alias{sequential.search.for.new.breakpoint.whole.interval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ooo
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sequential.search.for.new.breakpoint.whole.interval(dataset, model.type, resp.variable, expl.variables.without.segm.variable, segm.variable, model.parameters, angle.param = 0, fixed.breakpoints, last.end.point, start.point, end.point, min.distance, interval.length, time.unit, criterion, criterion.difference)
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
  \item{last.end.point}{
%%     ~~Describe \code{last.end.point} here~~
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
    segm.variable, model.parameters, angle.param = 0, fixed.breakpoints, 
    last.end.point, start.point, end.point, min.distance, interval.length, 
    time.unit, criterion, criterion.difference) 
{
    last.fixed.breakpoint <- if (identical(fixed.breakpoints, 
        numeric())) 
        0
    else fixed.breakpoints[length(fixed.breakpoints)]
    found.new.breakpoint <- FALSE
    temp.end.point <- max(last.end.point, last.fixed.breakpoint + 
        2 * min.distance)
    transf.temp.end.point <- transform.breakpoint(start.point = as.POSIXct(start.point), 
        breakpoint = temp.end.point, interval.length = interval.length, 
        time.unit = time.unit)$transf.breakpoint
    temp.dataset <- eval(parse(text = paste("subset(dataset, ", 
        segm.variable, " <=transf.temp.end.point)")))
    good.breakpoint.models <- list()
    while (transf.temp.end.point <= end.point) {
        alpha <- search.for.new.breakpoint.fixed.interval(dataset = temp.dataset, 
            model.type = model.type, resp.variable = resp.variable, 
            expl.variables.without.segm.variable = expl.variables.without.segm.variable, 
            segm.variable = segm.variable, model.parameters = model.parameters, 
            angle.param = angle.param, fixed.breakpoints = fixed.breakpoints, 
            start.point = start.point, end.point = transf.temp.end.point, 
            min.distance = min.distance, interval.length = interval.length, 
            time.unit = time.unit, criterion = criterion, criterion.difference = criterion.difference)
        if (alpha$found.new.breakpoint == TRUE) {
            found.new.breakpoint <- TRUE
            new.best.model <- alpha$best.model
            new.best.model.dataset <- temp.dataset
            new.end.point <- temp.end.point
            new.fixed.breakpoints <- alpha$new.fixed.breakpoints
            good.breakpoint.models <- append(good.breakpoint.models, 
                alpha$good.breakpoint.models)
            break
        }
        temp.end.point <- temp.end.point + 1
        transf.temp.end.point <- transform.breakpoint(start.point = start.point, 
            breakpoint = temp.end.point, interval.length = interval.length, 
            time.unit = time.unit)$transf.breakpoint
        temp.dataset <- eval(parse(text = paste("subset(dataset, ", 
            segm.variable, " <transf.temp.end.point)")))
    }
    if (found.new.breakpoint == FALSE) {
        new.end.point <- temp.end.point
        transf.breakpoints <- transform.breakpoint(start.point = start.point, 
            breakpoint = fixed.breakpoints, interval.length = interval.length, 
            time.unit = time.unit)$transf.breakpoint
        new.best.model.dataset <- dataset
        temp.fit <- fit.model.with.breakpoints(dataset = new.best.model.dataset, 
            model.type = model.type, resp.variable = resp.variable, 
            expl.variables.without.segm.variable = expl.variables.without.segm.variable, 
            segm.variable = segm.variable, breakpoints = transf.breakpoints, 
            model.parameters = model.parameters, angle.param = angle.param)
        new.best.model <- list(number.of.breakpoints = length(fixed.breakpoints), 
            breakpoints = fixed.breakpoints, transf.breakpoints = transf.breakpoints, 
            segm.var.coeffs = temp.fit$segm.var.coeffs, segm.var.significance = temp.fit$segm.var.significance, 
            aic = AIC(temp.fit$model), bic = AIC(temp.fit$model, 
                k = log(nrow(dataset))), formula.to.model = temp.fit$formula.to.model)
        good.breakpoint.models <- list(new.best.model)
        new.fixed.breakpoints <- fixed.breakpoints
    }
    return(list(found.new.breakpoint = found.new.breakpoint, 
        new.end.point = new.end.point, new.fixed.breakpoints = new.fixed.breakpoints, 
        new.best.model.dataset = new.best.model.dataset, new.best.model = new.best.model, 
        good.breakpoint.models = good.breakpoint.models))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
