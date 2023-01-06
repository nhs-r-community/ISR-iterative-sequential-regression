\name{summary.of.output}
\alias{summary.of.output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ooo
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.of.output(list.of.models, dataset, model.type, resp.variable, expl.variables.without.segm.variable, segm.variable, model.parameters, angle.param = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list.of.models}{
%%     ~~Describe \code{list.of.models} here~~
}
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
function (list.of.models, dataset, model.type, resp.variable, 
    expl.variables.without.segm.variable, segm.variable, model.parameters, 
    angle.param = 0) 
{
    final <- list.of.models[[length(list.of.models)]]$best.model
    temp <- fit.model.with.breakpoints(dataset = dataset, model.type = model.type, 
        resp.variable = resp.variable, expl.variables.without.segm.variable = expl.variables.without.segm.variable, 
        segm.variable = segm.variable, breakpoints = final$transf.breakpoints, 
        model.parameters = model.parameters, angle.param = angle.param)
    best.model <- temp$model
    best.dataset <- temp$dataset
    temp.colnames <- c("breakpoint", "left.CI", "right.CI", "when.noticed")
    summary.dataset <- as.data.frame(t(data.frame(rep(NA, length(temp.colnames)))))
    colnames(summary.dataset) <- temp.colnames
    temp.breakpoint <- min(dataset[, segm.variable])
    summary.dataset <- rbind(summary.dataset, c(temp.breakpoint, 
        rep(NA, length(temp.colnames) - 1)))
    if (length(list.of.models) > 1) {
        for (i in 1:(length(list.of.models) - 1)) {
            m <- list.of.models[[i]]
            number.of.breakpoints <- length(m$best.model$transf.breakpoints)
            temp.breakpoint <- m$best.model$transf.breakpoints[number.of.breakpoints]
            all.breakpoints <- temp.breakpoint
            for (m2 in m$good.breakpoint.models) all.breakpoints <- append(all.breakpoints, 
                (m2$transf.breakpoints[number.of.breakpoints]))
            temp.left.CI <- min(all.breakpoints)
            temp.right.CI <- max(all.breakpoints)
            temp.when.noticed <- m$end.point
            summary.dataset <- rbind(summary.dataset, c(temp.breakpoint, 
                temp.left.CI, temp.right.CI, temp.when.noticed))
        }
    }
    summary.dataset <- summary.dataset[2:nrow(summary.dataset), 
        ]
    for (c in c("breakpoint", "left.CI", "right.CI", "when.noticed")) summary.dataset[, 
        c] <- as.POSIXct(summary.dataset[, c], origin = "1970-01-01")
    best.coeffs <- as.data.frame(temp$summary.of.output[(nrow(temp$summary.of.output) - 
        length(final$transf.breakpoints)):nrow(temp$summary.of.output), 
        ])
    if (identical(final$transf.breakpoints, character(0))) 
        best.coeffs <- as.data.frame(t(best.coeffs))
    best.coeffs$left.CI <- best.coeffs[, 1] - 1.96 * best.coeffs[, 
        2]
    best.coeffs$right.CI <- best.coeffs[, 1] + 1.96 * best.coeffs[, 
        2]
    summary.dataset$angle <- best.coeffs[, 1]
    summary.dataset$angle.left.CI <- best.coeffs[, 5]
    summary.dataset$angle.right.CI <- best.coeffs[, 6]
    summary.dataset$pvalue <- best.coeffs[, 4]
    return(list(best.model = best.model, dataset = best.dataset, 
        summary.dataset = summary.dataset))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
