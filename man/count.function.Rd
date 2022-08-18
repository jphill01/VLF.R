\name{count.function}
\alias{count.function}
\title{
Nucleotide Count
}
\description{
Counts the number of each dNTP in each position of an aligned barcode matrix.
}
\usage{
count.function(nucleotides, spec.no, seqlength)
}
\arguments{
  \item{nucleotides}{
A matrix of aligned DNA barcode sequences.  DNA sequences should start at the third column of the matrix, while the first column contains a unique specimen identifier and the second column contains the species name.
}
  \item{spec.no}{
The number of specimen/sequences in the nucldeotide matrix.
}
  \item{seqlength}{
The length of the nucleotide sequences.
}
}
\value{
A matrix containing the number of each nucleotide in each position of the sequence.  Each row is a different dNTP count, while the columns represent the sequence position.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{data(birds)
Nuc.count <- count.function(birds, specimen.Number, 648)}
}