\name{matched.codon.position}
\alias{matched.codon.position}
\title{
Codon Position of Matching aa and ntVLFs
}
\description{
Counts which codon positions of ntVLFs lead to the concordant aaVLF.
}
\usage{
matched.codon.position(matched)
}
\arguments{
  \item{matched}{
A list of the nucleotide position of concordant ntVLFs and their associated aaVLFs.
}
}
\details{
The argument matched can be calculated using the function overall.matched.
}
\value{
A vector containing the number of concordant VLFs caused by ntVLFs in each codon position.
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
position_matchingAA <- aminoAcid.matching.positions(VLF_match[[1]], 216)
matching_comparison <- overall.matched(position_matchingNuc, position_matchingAA, 648, 216)
matching_codons <- matched.codon.position(matching_comparison)}
}
