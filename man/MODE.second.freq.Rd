\name{MODE.second.freq}
\alias{MODE.second.freq}
\title{
Second Modal Frequency}
\description{
Calculates the frequencies of the nucleotides that occur second most often in a matrix of sequences.
}
\usage{
MODE.second.freq(freq, seqlength)
}
\arguments{
  \item{freq}{
Frequency matrix for nucleotides.
}
  \item{seqlength}{
Length of nucleotide sequences.
}
}
\details{
The argument freq can be calculated using the function ffrequency.matrix.function.
}
\value{
A vector containing the frequencies of the nucleotide sequence that occurs second most often.
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
Bird_second.modal.frequencies <- MODE.second.freq(frequency.matrix, 648)}
}