% Generated by roxygen2 (4.0.2): do not edit by hand
\name{genWellNumbers}
\alias{genWellNumbers}
\title{Well Number Generator}
\usage{
genWellNumbers(padding = 2, nrows = 8, ncols = 12)
}
\arguments{
\item{padding=2}{(overall number of digits, leading digits are filled with zeros)}

\item{nrows=8}{(number of rows on plate)}

\item{ncols=12}{(number of rows on plate)}
}
\value{
as.factor(well_numbers) (factors of well numbers)
}
\description{
Generator for microtiter plate well coordinates
}
\details{
genWellNumbers
}
\examples{
genWellNumbers(padding=1)
             >  A1, A2, ...., H12
           genWellNumbers(padding=3)
             > A001, A002, ..., H012
}
\keyword{microtiter}
\keyword{plate}
\keyword{plates}
\keyword{readers,}

