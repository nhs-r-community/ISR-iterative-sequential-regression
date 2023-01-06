\name{fit.model}
\alias{fit.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
hhh
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fit.model(dataset, model.type, resp.variable, expl.variables, model.parameters)
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
  \item{expl.variables}{
%%     ~~Describe \code{expl.variables} here~~
}
  \item{model.parameters}{
%%     ~~Describe \code{model.parameters} here~~
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
function (dataset, model.type, resp.variable, expl.variables, 
    model.parameters) 
{
    formula.to.model <- paste(paste(resp.variable, collapse = "+"), 
        "~", paste(expl.variables, collapse = "+"))
    if (model.parameters != "" && !is.null(model.parameters)) 
        formula.to.model <- paste(formula.to.model, ",", model.parameters)
    formula.to.model <- paste(model.type, "(", formula.to.model, 
        ",", "data=dataset", ")")
    model <- eval(parse(text = formula.to.model))
    return(list(coefficients = summary(model, se = "iid")$coefficients[, 
        1], significance = summary(model, se = "iid")$coefficients[, 
        4], formula.to.model = formula.to.model, model = model))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
