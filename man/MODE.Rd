\name{MODE}
\alias{MODE}
\title{
Modal Sequence
}
\description{
Calculates the nucleotide sequence that occurs most often in a matrix of sequences.
}
\usage{
MODE(freq, seqlength)
}
\arguments{
  \item{freq}{
Frequenct matrix for nucleotides.
}
  \item{seqlength}{
Length of nucleotide sequence.
}
}
\details{
The argument freq can be calculated using the function ffrequency.matrix.function.
}
\value{
A vector containing the first modal sequence.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{data(birds)
species.names <- birds[,2]
specimen.Number <- nrow(birds)
rownames(birds) <- species.names
Nuc.count <- count.function(birds, specimen.Number, 648)
frequency.matrix <- ffrequency.matrix.function(Nuc.count, 648)
birdSpec.freq <- specimen.frequencies(frequency.matrix, birds, specimen.Number, species.names, 648)
nucleotide.modalSequence <- MODE(frequency.matrix, 648)}
}