\name{Decile.Plot}
\alias{Decile.Plot}
\title{
VLF Decile Plot}
\description{
Creates a plot of VLF distributions summed for each decile segment.
}
\usage{
Decile.Plot(VLF, seqlength)
}
\arguments{
  \item{VLF}{
A list of VLFs in each barcode position.  May be a matrix containing vectors of singleton and shared VLFs, or can be a single vector of total VLFs.
}
  \item{seqlength}{
The length of the sequence.  Usually 648 for nucleotide sequences and 216 for amino acid sequences.
}
}
\value{
A barplot containing the sum of the VLFs for each decile barcode segment.
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
Bird_position_VLFcount <- VLF.count.pos(birdSpec.freq, 0.001, 648)
Decile.Plot(Bird_position_VLFcount, 648)}
}