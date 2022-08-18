\name{VLF.nucleotides}
\alias{VLF.nucleotides}
\title{
Nucleotide VLF Convert}
\description{
Converts a matrix of nucleotide frequencies for each specimen into a matrix of nucleotides for each specimen.
}
\usage{
VLF.nucleotides(convert.matrix, seq.matrix, seqlength)
}
\arguments{
  \item{convert.matrix}{
A matrix consisting of only very low frequency cariant frequencies for each specimen, and NAs in all other positions of the sequence.
}
  \item{seq.matrix}{
A matrix of DNA sequences.
}
  \item{seqlength}{
The length of the sequences.
}
}
\details{
The argument convert.matrix can be calculated using the function VLF.convert.matrix.
}
\value{
A matrix containing only ntVLFs in each position of the sequences, and NAs in all other positions.
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
bird_VLFconvert <- VLF.convert.matrix(birds, birdSpec.freq, 0.001, 648)
bird_VLFnuc <- VLF.nucleotides(bird_VLFconvert, birds, 648)}
}