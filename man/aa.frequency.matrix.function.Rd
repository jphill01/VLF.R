\name{aa.frequency.matrix.function}
\alias{aa.frequency.matrix.function}
\title{
Amino Acid Frequency Matrix}
\description{
Calculates the frequency of each amino acid.
}
\usage{
aa.frequency.matrix.function(aa.count, seqlength)
}
\arguments{
  \item{aa.count}{
A matrix containing the number of each amino acid in each position.
}
  \item{seqlength}{
The length of the amino acid sequence
}
}
\details{
The aa.count argument can be calculated using the function aa.count.function
}
\value{
A matrix of the frequencies for each amino acid in each position of the barcode sequence.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{data(birds_aminoAcids)
birds_aminoAcid_count <- aa.count.function(birds_aminoAcids, 216)
aminoAcid_frequency.Matrix <- aa.frequency.matrix.function(birds_aminoAcid_count, 216)}
}