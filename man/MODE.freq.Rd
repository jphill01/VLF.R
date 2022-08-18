\name{MODE.freq}
\alias{MODE.freq}
\title{
Modal Frequencies}
\description{
Returns the frequencies of the nucleotides in each position of nucleotide sequence that occurs most often.
}
\usage{
MODE.freq(freq, seqlength)
}
\arguments{
  \item{freq}{
Frequency matrix for nucleotides.
}
  \item{seqlength}{
Length of the nucleotide sequence.
}
}
\details{
The argument freq can be calculated using the function ffrequency.matrix.function.
}
\value{
A vector of frequencies for the first modal sequence.
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
Bird_first.modal.frequencies <- MODE.freq(frequency.matrix, 648)}
}