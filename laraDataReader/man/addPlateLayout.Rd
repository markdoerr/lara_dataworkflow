% Generated by roxygen2 (4.0.2): do not edit by hand
\name{addPlateLayout}
\alias{addPlateLayout}
\title{Adding Plate Layout Information to reader data frame}
\usage{
addPlateLayout(reader_df = NULL, barcode = "0000", use_db_layout = FALSE,
  set_Value_NA = FALSE, set_Slope_NA = FALSE)
}
\arguments{
\item{reader_df=NULL}{(data.frame) - plate reader data frame from import reader file}

\item{barcode="0000"}{(string)}

\item{set_Value_NA=TRUE}{(boolean) - set values to NA for empty plates of type 0}

\item{set_Slope_NA=FALSE}{(boolean) - set slope and intercept values to NA for empty plates of type 0}
}
\description{
Convenient function to load plate layout from a layout file or database
             and merging it with the reader data frame.
}
\details{
addPlateLayout
}
\note{
todo : merging bug with multiple measurements per data frame (e.g. groth data)
}
\keyword{layout}
\keyword{plate}
\keyword{readers,}

