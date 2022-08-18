\name{aa.count.function}
\alias{aa.count.function}
\title{
Amino Acid Count
}
\description{
Counts the number of each amino acid in each positino of the barcode.
}
\usage{
aa.count.function(aminoAcids, seqlength)
}
\arguments{
  \item{aminoAcids}{
A matrix of barcode amino acid sequences.
}
  \item{seqlength}{
Length of the amino acid sequences.
}
}
\details{
The first and second column of the aminoAcid argument must contain the unique specimen identifier and the species name, respectively, followed by the amino acid sequence.}
\value{
A matrix containing the number of each amino acid in each position of the sequence.  Each row is a different amino acid count, while the columns represent the sequence position.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{data(birds_aminoAcids)
birds_aminoAcid_count <- aa.count.function(birds_aminoAcids, 216)}
}