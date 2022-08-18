\name{VLF.count.pos}
\alias{VLF.count.pos}
\title{
VLF Count for Sequence Positions
}
\description{
Calculates the number of very low frequency variants (VLFs) in each position in a matrix of sequences.
}
\usage{
VLF.count.pos(freq, p, seqlength)
}
\arguments{
  \item{freq}{
A matrix of frequencies for each specimen.
}
  \item{p}{
A very low frequency variant designation cut off frequency.  Any frequency in the freq matrix below this value is considered to be a very low frequency variant.
}
  \item{seqlength}{
The length of the sequences.
}
}
\details{
The argument freq can be calculated using the specimen.frequencies function.
}
\value{
A vector containing the number of VLFs for each position in the sequence.
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
Bird_specimen_VLFcount <- VLF.count.spec(birdSpec.freq, 0.001, 648)
Bird_position_VLFcount <- VLF.count.pos(birdSpec.freq, 0.001, 648)}
}