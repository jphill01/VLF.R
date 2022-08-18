\name{find.singles}
\alias{find.singles}
\title{
Single and Shared VLF Find}
\description{
Calculates the number of singleton and shared VLFs for each position of the nucleotide, by first seeing if there is only one specimen for a species, and then calling on the compare() function to calculate the number of singleton and shared VLFs for those species with multiple specimen.}
\usage{
find.singles(species, seqlength)
}
\arguments{
  \item{species}{
A list of sequences separated byh species name.  Each entry in the list contains a matrix of sequences from the same species.
}
  \item{seqlength}{
Length of the nucleotide sequence.
}
}
\details{
The argument species can be calculated using the separate function.}
\value{
A matrix containing the number of singleton and shared ntVLFs in each position of the barcode.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{#Nucleotide VLF analysis
data(birds)
species.names <- birds[,2]
specimen.Number <- nrow(birds)
rownames(birds) <- species.names
Nuc.count <- count.function(birds, specimen.Number, 648)
frequency.matrix <- ffrequency.matrix.function(Nuc.count, 648)
birdSpec.freq <- specimen.frequencies(frequency.matrix, birds, specimen.Number, species.names, 648)
Bird_specimen_VLFcount <- VLF.count.spec(birdSpec.freq, 0.001, 648)
bird_VLFconvert <- VLF.convert.matrix(birds, birdSpec.freq, 0.001, 648)
bird_VLFnuc <- VLF.nucleotides(bird_VLFconvert, birds, 648)
bird_VLFreduced <- VLF.reduced(bird_VLFnuc, Bird_specimen_VLFcount, 648)
bird_species <- separate(bird_VLFreduced)
birds_singleAndShared <- find.singles(bird_species, 648)}
}