\name{master.isr}
\alias{master.isr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
jjj
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
master.isr(dataset, resp.variable, resp.variable.label, step = FALSE, t.boxcox, factor.expl.variables, cts.expl.variables, offset.variable = "", segm.variable, segm.variable.label, model.type, model.parameters, angle.param = 0, start.point, end.point, time.unit, interval.length, min.distance, criterion, criterion.difference)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
}
  \item{resp.variable}{
%%     ~~Describe \code{resp.variable} here~~
}
  \item{resp.variable.label}{
%%     ~~Describe \code{resp.variable.label} here~~
}
  \item{step}{
%%     ~~Describe \code{step} here~~
}
  \item{t.boxcox}{
%%     ~~Describe \code{t.boxcox} here~~
}
  \item{factor.expl.variables}{
%%     ~~Describe \code{factor.expl.variables} here~~
}
  \item{cts.expl.variables}{
%%     ~~Describe \code{cts.expl.variables} here~~
}
  \item{offset.variable}{
%%     ~~Describe \code{offset.variable} here~~
}
  \item{segm.variable}{
%%     ~~Describe \code{segm.variable} here~~
}
  \item{segm.variable.label}{
%%     ~~Describe \code{segm.variable.label} here~~
}
  \item{model.type}{
%%     ~~Describe \code{model.type} here~~
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
  \item{time.unit}{
%%     ~~Describe \code{time.unit} here~~
}
  \item{interval.length}{
%%     ~~Describe \code{interval.length} here~~
}
  \item{min.distance}{
%%     ~~Describe \code{min.distance} here~~
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
function (dataset, resp.variable, resp.variable.label, step = FALSE, 
    t.boxcox, factor.expl.variables, cts.expl.variables, offset.variable = "", 
    segm.variable, segm.variable.label, model.type, model.parameters, 
    angle.param = 0, start.point, end.point, time.unit, interval.length, 
    min.distance, criterion, criterion.difference) 
{
    transformation <- if (t.boxcox) 
        "boxcox"
    else "identity"
    temp <- transf.resp.variable(dataset = dataset, var = resp.variable, 
        transformation = transformation)
    new.dataset <- temp$output.dataset
    new.resp.variable <- temp$transf.var
    lambda <- temp$lambda
    expl.variables.without.segm.variable <- character()
    expl.variables.without.segm.variable <- append(expl.variables.without.segm.variable, 
        factor.expl.variables)
    if (!identical(step, FALSE)) {
        temp <- add.step(dataset = new.dataset, var = resp.variable, 
            time.var = segm.variable, step = step)
        new.dataset <- temp$output.dataset
        expl.variables.without.segm.variable <- append(expl.variables.without.segm.variable, 
            temp$step.var)
    }
    expl.variables.without.segm.variable <- append(expl.variables.without.segm.variable, 
        cts.expl.variables)
    if ((!identical(offset.variable, "")) & (!identical(offset.variable, 
        character()))) {
        log.offset.variable <- paste("offset(log(", offset.variable, 
            "))", sep = "")
        expl.variables.without.segm.variable <- append(expl.variables.without.segm.variable, 
            log.offset.variable)
    }
    expl.variables.without.segm.variable <- expl.variables.without.segm.variable[which(expl.variables.without.segm.variable != 
        "")]
    isr.models <- sequential.grid.search.save.models(dataset = new.dataset, 
        model.type = model.type, resp.variable = new.resp.variable, 
        expl.variables.without.segm.variable = expl.variables.without.segm.variable, 
        segm.variable = segm.variable, model.parameters = model.parameters, 
        angle.param = angle.param, start.point = start.point, 
        end.point = end.point, min.distance = min.distance, interval.length = interval.length, 
        time.unit = time.unit, criterion = criterion, criterion.difference = criterion.difference)
    output <- summary.of.output(dataset = new.dataset, list.of.models = isr.models$models, 
        model.type = model.type, resp.variable = new.resp.variable, 
        expl.variables.without.segm.variable = expl.variables.without.segm.variable, 
        segm.variable = segm.variable, model.parameters = model.parameters, 
        angle.param = angle.param)
    require(ggplot2)
    plot.fit <- plot.fit(dataset = output$dataset, model = output$best.model, 
        summary.dataset = output$summary.dataset, original.resp.variable = resp.variable, 
        resp.variable.label = resp.variable.label, lambda = lambda, 
        factor.expl.variables = factor.expl.variables, cts.expl.variables = cts.expl.variables, 
        offset.variable = offset.variable, segm.variable = segm.variable, 
        segm.variable.label = segm.variable.label)
    return(list(list.of.models = isr.models$models, output = output, 
        plot.fit = plot.fit))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
