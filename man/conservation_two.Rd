\name{conservation_two}
\alias{conservation_two}
\title{
First and Second Modal Conserved}
\description{
Calculates the conservation of the nucleotides that occur first and second most often in a matrix of sequences.
}
\usage{
conservation_two(modal1, modal2, p, seqlength)
}
\arguments{
  \item{modal1}{
A vector of the frequencies for the nucleotides in the first modal sequence.
}
  \item{modal2}{
A vector of the frequencies for the nucleotides that occur second most often.
}
  \item{p}{
A conservation value for the nucleotide frequencies to be compared to.
}
  \item{seqlength}{
The nucleotide sequence length.
}
}
\details{
The argument modal1 can be calculated using the function MODE.freq.  The argument modal2 can be calculated using the function MODE.second.freq.
}
\value{
A vector that contains how many nucleotides from the first and second modal sequences are conserved at the specified conservation level for each codon position.
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
Bird_second.modal.frequencies <- MODE.second.freq(frequency.matrix, 648)
FirstAndSecond_conserved_99.9 <- conservation_two(Bird_first.modal.frequencies, 
    Bird_second.modal.frequencies, 0.999, 648)}
}
