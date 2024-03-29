% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PACells.R
\name{getMotifs}
\alias{getMotifs}
\title{Get the TF Motifs}
\usage{
getMotifs(
  species = "Homo sapiens",
  database = c("JASPAR", "cisBP")[1],
  JASPAR = JASPAR2022
)
}
\arguments{
\item{species}{The species source for the sequences, in Latin (Homo sapiens) or NCBI tax IDs (9606).}

\item{database}{A database of motif.}
}
\value{
A PFMatrixList object.
}
\description{
This function gets the TF Motifs from "JASPAR" or "cisBP" database.
}
