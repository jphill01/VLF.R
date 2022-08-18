\name{aminoAcid.matching.positions}
\alias{aminoAcid.matching.positions}
\title{
Matching Amino Acid Positions
}
\description{
Gives the position and residue of the amino acid VLFs in a matrix containing amino acid VLFs and NAs.}
\usage{
aminoAcid.matching.positions(matchAA, aalength)
}
\arguments{
  \item{matchAA}{
A matrix containing aaVLFs and NAs in all other positions of the sequences.}
  \item{aalength}{
Amino acid sequence length.}
}
\details{
The argument matchAA can be calculated usingthe find.matching function and taking the first argument from the returned value.}
\value{
A list for each aaVLF with a matching specimen identifier to a ntVLF.  The first position in each entry of the list contains the specimen identifier, the second position contains the species name, the third position contains the sequence position of the aaVLF, and the fourth position contains the aaVLF.
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
birds_singleAndShared <- find.singles(bird_species, 648)

#Amino Acid VLF Analysis
data(birds_aminoAcids)
birds_aminoAcid_speciesNames <- birds_aminoAcids[,2]
aminoAcids_specimenNumber <- nrow(birds_aminoAcids)
birds_aminoAcid_count <- aa.count.function(birds_aminoAcids, 216)
aminoAcid_frequency.Matrix <- aa.frequency.matrix.function(birds_aminoAcid_count, 216)
bird_aminoAcid_frequencies <- aa.specimen.frequencies(aminoAcid_frequency.Matrix, birds_aminoAcids,
    birds_aminoAcid_speciesNames, 216)
aminoAcid_Modal <- aa.MODE(aminoAcid_frequency.Matrix, 216)
birds_aminoAcid_specimenVLFcount <- aa.VLF.count.spec(bird_aminoAcid_frequencies, 0.001, 216)
birds_aaVLFconvert <- aa.VLF.convert.matrix(birds_aminoAcids, bird_aminoAcid_frequencies, 0.001, 
    216)
birds_aminoAcidVLFs <- VLF.aminoAcids(birds_aaVLFconvert, birds_aminoAcids, 216)
birds_aaVLFreduced <- aa.VLF.reduced(birds_aminoAcidVLFs, birds_aminoAcid_specimenVLFcount, 216)
birds_aaSpecies <- separate(birds_aaVLFreduced)
birds_aminoAcid_singleAndShared <- aa.find.singles(birds_aaSpecies, 216)

#Concordance Analysis
VLF_match <- find.matching(bird_VLFreduced, birds_aaVLFreduced, 648, 216)
position_matchingNuc <- nucleotide.matching.positions(VLF_match[[2]], 648)
position_matchingAA <- aminoAcid.matching.positions(VLF_match[[1]], 216)}
}
