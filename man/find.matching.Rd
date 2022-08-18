\name{find.matching}
\alias{find.matching}
\title{
Find Matching ntVLF and aaVLF Specimen.
}
\description{
Compares a list of aaVLF and ntVLF matrices for common specimen identifiers.
}
\usage{
find.matching(NucleotideList, AminoAcidList, nuclength, aalength)
}
\arguments{
  \item{NucleotideList}{
Matrix of VLF nucleotide sequences containing the only the nucleotidies that are VLFs and NAs in the other positions of the sequences.
}
  \item{AminoAcidList}{
Matrix of VLF amino acid sequences containing only the aaVLFs and NAs in the other positions of the sequences.
}
  \item{nuclength}{
Length of the nucleotide sequence (should be 3X the length of the amino acid sequence).
}
  \item{aalength}{
Length of the amino acid sequence (should be 1/3 the length of the nucleotide sequence).}
}
\details{
The argument NucleotideList can be calculated using the VLF.convert.matrix, VLF.nucleotides, and VLF.reduced functions.  The argument AminoAcidList can be calculated using the aa.VLF.convert.matrix, VLF.aminoAcids, and aa.VLF.reduced functions.}
\value{
A list containing matrices of aaVLFs in the first position and ntVLFs in the second position who have matching specimen identifiers.
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
VLF_match <- find.matching(bird_VLFreduced, birds_aaVLFreduced, 648, 216)}
}
