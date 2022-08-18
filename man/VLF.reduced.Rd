\name{VLF.reduced}
\alias{VLF.reduced}
\title{
Reduced VLF Matrix
}
\description{
Reduces a matrix of very low frequency variants (VLFs) so that only those specimen that contain VLFs remain in the matrix.
}
\usage{
VLF.reduced(NA.matrix, sCount, seqlength)
}
\arguments{
  \item{NA.matrix}{
A matrix with values for very low frequency variants (VLFs) for each specimen and NAs in all other positions of the sequence.
}
  \item{sCount}{
A vector of the very low frequency variant (VLF) counts for each specimen in the NA.matrix.
}
  \item{seqlength}{
The length of the sequences.
}
}
\details{
The argument NA.matrix can be calculated using the function VLF.convert.matrix and VLF.nucleotides.  The argument sCount can be calculated using the function VLF.count.spec.
}
\value{
A matrix containing only those specimen that have VLFs, and only VLFs in their positions in the sequence, all other positions contain NAs.
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
bird_VLFconvert <- VLF.convert.matrix(birds, birdSpec.freq, 0.001, 648)
bird_VLFnuc <- VLF.nucleotides(bird_VLFconvert, birds, 648)
bird_VLFreduced <- VLF.reduced(bird_VLFnuc, Bird_specimen_VLFcount, 648)}
}