\name{separate}
\alias{separate}
\title{
Separate Specimen by Species Names}
\description{
Separates specimen into lists by species names.
}
\usage{
separate(x)
}
\arguments{
  \item{x}{
A matrix of sequences, usually reduced sequences containing only VLFs, where the second position of the matrix contains the species name for the specimen.
}
}
\details{
If the argument x needs to be a reduced matrix, it can be calculated using the function VLF.reduced.
}
\value{
A list containing a matrix of sequences for each species.
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
bird_species <- separate(bird_VLFreduced)}
}