\name{ffrequency.matrix.function}
\alias{ffrequency.matrix.function}
\title{Frequency Matrix}
\description{
Calculates the frequency of each dNTP in each position of a nucleotide count matrix.
}
\usage{
ffrequency.matrix.function(count.matrix, seqlength)
}
\arguments{
  \item{count.matrix}{
A matrix of the counts for each dNTP from a matrix of aligned sequences.
}
  \item{seqlength}{
Length of sequences.
}
}
\details{
The argument count.matrix can be calculated using the function count.function.}
\value{
A matrix of the frequencies for each dNTP in each position of the barcode sequence.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{#Nucleotide VLF analysis
data(birds)
Nuc.count <- count.function(birds, specimen.Number, 648)
frequency.matrix <- ffrequency.matrix.function(Nuc.count, 648)}
}