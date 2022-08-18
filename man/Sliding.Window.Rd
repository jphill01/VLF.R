\name{Sliding.Window}
\alias{Sliding.Window}
\title{
Sliding Window
}
\description{
Creates a sliding window analysis plot for the VLFs in a matrix of sequences.
}
\usage{
Sliding.Window(VLF, seqlength, n = 30)
}
\arguments{
  \item{VLF}{
A vector of VLFs per position across the barcocde.  Can be a single vector of all VLFs per positions, or can be a matrix containing singleton and shared VLFs.
}
  \item{seqlength}{
Length of the barcode sequence.
}
  \item{n}{
The number of positions to average the window across (n = 30 by default).
}
}
\details{
The argument VLF can be calculated using the function VLF.count.pos for all VLFs, or find.singles for singleton and shared VLFs.
}
\value{
A sliding window plot for the VLFs in each position of the barcode averaged over n.
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
Sliding.Window(Bird_position_VLFcount, 648)}
}