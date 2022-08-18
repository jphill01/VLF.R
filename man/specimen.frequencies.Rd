\name{specimen.frequencies}
\alias{specimen.frequencies}
\title{
Specimen Nucleotide Frequencies
}
\description{
Converts a matrix of sequences into a matrix of nucleotide frequencies.
}
\usage{
specimen.frequencies(freq, seq.matrix, no.spec, spec.names, seqlength)
}
\arguments{
  \item{freq}{
Frequency matrix for nucleotides.
}
  \item{seq.matrix}{
Matrix of specimen sequences, where the sequence starts in the third position of the matrix and the first and second position contain the unique specimen identifier and the species name, respectively.
}
  \item{no.spec}{
The number of specimen in seq.matrix.
}
  \item{spec.names}{
A vector containing the names of the specimen in the seq.matrix, in the order they appear in the matrix.
}
  \item{seqlength}{
The length of the nucleotide sequence.
}
}
\details{
The argument freq can be calculated using the function ffrequency.matrix.function.  The number of specimen can be calculated by using the nrow() function on seq.matrix.
}
\value{
A matrix containing the unique specimen identifer in the first position, the species name in the second position, and the frequencies for each nucleotide in the sequences starting at the third position.
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
birdSpec.freq <- specimen.frequencies(frequency.matrix, birds, specimen.Number, species.names, 648)}
}