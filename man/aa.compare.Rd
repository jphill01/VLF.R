\name{aa.compare}
\alias{aa.compare}
\title{
Amino Acid Comparison
}
\description{
Compares amino acid very low frequency variants between specimen of the same species
}
\usage{
aa.compare(x, seqlength)
}
\arguments{
  \item{x}{
A list of amino acid sequences separated by species name.
}
\item{seqlength}{
The length of the amino acid sequences within the list.
}
}
\details{
The argument x can be calculated using the separate function.
}
\value{
A matrix containing two vectors, one with singleton VLF counts for each position of the sequence, and one with shared VLF counts for each position of the sequence.
}
\author{
Taryn B.T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{data(birds_aminoAcids)
birds_aminoAcid_speciesNames <- birds_aminoAcids[,2]
aminoAcids_specimenNumber <- nrow(birds_aminoAcids)
birds_aminoAcid_count <- aa.count.function(birds_aminoAcids, 216)
aminoAcid_frequency.Matrix <- aa.frequency.matrix.function(birds_aminoAcid_count, 216)
bird_aminoAcid_frequencies <- aa.specimen.frequencies(aminoAcid_frequency.Matrix, birds_aminoAcids,
    birds_aminoAcid_speciesNames, 216)
birds_aminoAcid_specimenVLFcount <- aa.VLF.count.spec(bird_aminoAcid_frequencies, 0.001, 216)
birds_aaVLFconvert <- aa.VLF.convert.matrix(birds_aminoAcids, bird_aminoAcid_frequencies, 0.001, 
    216)
birds_aminoAcidVLFs <- VLF.aminoAcids(birds_aaVLFconvert, birds_aminoAcids, 216)
birds_aaVLFreduced <- aa.VLF.reduced(birds_aminoAcidVLFs, birds_aminoAcid_specimenVLFcount, 216)
birds_aaSpecies <- separate(birds_aaVLFreduced)
birds_aminoAcid_singleAndShared <- aa.find.singles(birds_aaSpecies, 216)}
}
