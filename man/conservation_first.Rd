\name{conservation_first}
\alias{conservation_first}
\title{
First Modal Conserved
}
\description{
Calculates the conservation of the nucleotides that occur most often in a matrix of sequences/
}
\usage{
conservation_first(modal, p, seqlength)
}
\arguments{
  \item{modal}{
A vector of the frequencies of the nucleotides in the first modal sequences.
}
  \item{p}{
A conservation value for the nucleotide frequencies to be compared to.
}
  \item{seqlength}{
The length of the nucleotide sequence.
}
}
\details{
The argument modal can be calculated using the MODE.freq function.
}
\value{
A vector that contains how many nucleotides from the first modal sequence are conserved at the specified conservation level for each codon position.
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
nucleotide.modalSequence <- MODE(frequency.matrix, 648)
Bird_first.modal.frequencies <- MODE.freq(frequency.matrix, 648)
First_conserved_100 <- conservation_first(Bird_first.modal.frequencies, 1, 648)
First_conserved_99.9 <- conservation_first(Bird_first.modal.frequencies, 0.999, 648)}
}