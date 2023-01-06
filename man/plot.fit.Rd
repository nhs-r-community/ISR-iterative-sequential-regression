\name{plot.fit}
\alias{plot.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
lll
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.fit(dataset, model, summary.dataset, original.resp.variable, resp.variable.label, factor.expl.variables, cts.expl.variables, offset.variable = "", segm.variable, segm.variable.label, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{
%%     ~~Describe \code{dataset} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{summary.dataset}{
%%     ~~Describe \code{summary.dataset} here~~
}
  \item{original.resp.variable}{
%%     ~~Describe \code{original.resp.variable} here~~
}
  \item{resp.variable.label}{
%%     ~~Describe \code{resp.variable.label} here~~
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
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
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
function (dataset, model, summary.dataset, original.resp.variable, 
    resp.variable.label, factor.expl.variables, cts.expl.variables, 
    offset.variable = "", segm.variable, segm.variable.label, 
    lambda) 
{
    if (!(identical(cts.expl.variables, ""))) 
        for (c in cts.expl.variables) dataset[, c] <- mean(dataset[, 
            c], na.rm = TRUE)
    if (!identical(factor.expl.variables, "") & !identical(factor.expl.variables, 
        character())) {
        factors <- list()
        for (f in factor.expl.variables) factors[[f]] <- levels(dataset[[f]])
        factors.df <- expand.grid(factors)
        factors.df.long <- factors.df[rep(seq(nrow(factors.df)), 
            nrow(dataset)), ]
        if (length(factor.expl.variables) == 1) {
            factors.df.long <- as.data.frame(factors.df.long)
            colnames(factors.df.long) <- factor.expl.variables
        }
        dataset <- dataset[rep(seq(nrow(dataset)), nrow(factors.df)), 
            ]
        for (f in factor.expl.variables) dataset[[f]] <- factors.df.long[[f]]
        dataset$factors <- eval(parse(text = paste("paste(", 
            paste("dataset[['", factor.expl.variables, "']]", 
                sep = "", collapse = ","), ")", sep = "")))
    }
    dataset$predict <- predict(model, dataset, type = "response")
    dataset <- untransf.resp.variable(dataset = dataset, var = "predict", 
        lambda = lambda)$output.dataset
    dataset$yvar <- as.numeric(as.character(dataset[, original.resp.variable]))
    dataset$xvar <- dataset[, segm.variable]
    if (!identical(offset.variable, "")) {
        dataset$untransf.predict <- dataset$untransf.predict/dataset[[offset.variable]]
        dataset$yvar <- dataset$yvar/dataset[[offset.variable]]
    }
    p <- ggplot(dataset, aes(x = xvar))
    p <- p + geom_point(aes(y = yvar))
    if (!identical(factor.expl.variables, "") & !identical(factor.expl.variables, 
        character())) 
        p <- p + geom_line(aes(y = untransf.predict, colour = factors, 
            group = factors), size = 0.5)
    else p <- p + geom_line(aes(y = untransf.predict), size = 0.5, 
        colour = "blue")
    p <- p + xlab(segm.variable.label)
    p <- p + opts(axis.text.x = theme_text(angle = 45))
    p <- p + ylab(resp.variable.label)
    p <- p + opts(axis.ticks = theme_blank())
    if (nrow(summary.dataset) > 1) {
        summary.dataset2 <- summary.dataset[2:nrow(summary.dataset), 
            ]
        summary.dataset2$ycoord <- min(dataset$yvar)
        p <- p + geom_point(data = summary.dataset2, aes(x = breakpoint, 
            y = ycoord), colour = "blue")
        p <- p + geom_errorbarh(data = summary.dataset2, aes(xmin = left.CI, 
            xmax = right.CI, x = breakpoint, y = ycoord), colour = "blue", 
            height = 0)
        p <- p + geom_point(data = summary.dataset2, aes(x = when.noticed, 
            y = ycoord), colour = "blue", shape = 3)
    }
    return(list(p = p))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
