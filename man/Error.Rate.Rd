\name{Error.Rate}
\alias{Error.Rate}
\title{
Error Rate}
\description{
Calculates shared, singleton, and total second codon position error rate in a matrix of sequences.
}
\usage{
Error.Rate(single, shared, spec, seqlength)
}
\arguments{
  \item{single}{
A vector of singleton very low frequency variant (VLF) counts for each position in the sequence.}
  \item{shared}{
A vector of shared very low frequenct variant (VLF) counts for each position in the sequence.}
  \item{spec}{
The number of specimen being considered in the dataset.}
  \item{seqlength}{
The length of the barcode sequence.}
}
\details{
The arguments single and shared can be calculated simultaneously using the find.singles function.  The spec argument can be calculated by using the nrow() function on the sequence matrix.}
\value{
A vector containing the single, shared, and total error rate based on the number of second position VLFs.
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
Bird_position_VLFcount <- VLF.count.pos(birdSpec.freq, 0.001, 648)
bird_VLFconvert <- VLF.convert.matrix(birds, birdSpec.freq, 0.001, 648)
bird_VLFnuc <- VLF.nucleotides(bird_VLFconvert, birds, 648)
bird_VLFreduced <- VLF.reduced(bird_VLFnuc, Bird_specimen_VLFcount, 648)
bird_species <- separate(bird_VLFreduced)
birds_singleAndShared <- find.singles(bird_species, 648)
Bird_error <- Error.Rate(birds_singleAndShared[1,], birds_singleAndShared[2,], specimen.Number, 648)}
}