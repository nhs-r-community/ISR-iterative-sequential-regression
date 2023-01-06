\name{fit.model.with.breakpoints}
\alias{fit.model.with.breakpoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
iii
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fit.model.with.breakpoints(dataset, model.type, resp.variable, expl.variables.without.segm.variable, segm.variable, breakpoints, model.parameters, angle.param = 0)
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
  \item{breakpoints}{
%%     ~~Describe \code{breakpoints} here~~
}
  \item{model.parameters}{
%%     ~~Describe \code{model.parameters} here~~
}
  \item{angle.param}{
%%     ~~Describe \code{angle.param} here~~
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
    segm.variable, breakpoints, model.parameters, angle.param = 0) 
{
    temp <- break.segm.variable(angle.param = angle.param, dataset = dataset, 
        segm.variable = segm.variable, breakpoints = breakpoints)
    dataset <- temp$output.dataset
    extra.expl.variables <- temp$split.variable
    expl.variables <- append(expl.variables.without.segm.variable, 
        extra.expl.variables)
    expl.variables <- expl.variables[which(expl.variables != 
        "")]
    fit <- fit.model(dataset = dataset, model.type = model.type, 
        resp.variable = resp.variable, expl.variables = expl.variables, 
        model.parameters = model.parameters)
    summary.of.output <- summary(fit$model, se = "iid")$coefficients
    segm.var.coeffs <- c()
    for (v in extra.expl.variables) segm.var.coeffs <- append(segm.var.coeffs, 
        fit$coefficients[v])
    segm.var.significance <- c()
    for (v in extra.expl.variables) segm.var.significance <- append(segm.var.significance, 
        fit$significance[v])
    return(list(dataset = dataset, coefficients = fit$coefficients, 
        significance = fit$significance, segm.var.coeffs = segm.var.coeffs, 
        segm.var.significance = segm.var.significance, formula.to.model = fit$formula.to.model, 
        model = fit$model, summary.of.output = summary.of.output))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
