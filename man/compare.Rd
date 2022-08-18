\name{compare}
\alias{compare}
\title{
Compare VLFs within Species
}
\description{
Compares VLFs between specimen of the same species.
}
\usage{
compare(x, seqlength)
}
\arguments{
  \item{x}{
A list of sequences separated by species name.  Each entry in the list contains a matrix of sequences from the same species.
}
  \item{seqlength}{
Length of the sequences.
}
}
\details{
List of sequences by species names, x, can be created using the separate function
}
\value{
A matrix containing two vectors, one with singleton VLF counts for each position of the sequence, and one with shared VLF counts for each position of the sequence.
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
bird_VLFreduced <- VLF.reduced(bird_VLFnuc, Bird_specimen_VLFcount, 648)
bird_species <- separate(bird_VLFreduced)

#The compare function is called on from within the find.singles function
birds_singleAndShared <- find.singles(bird_species, 648)}
}