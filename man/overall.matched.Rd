\name{overall.matched}
\alias{overall.matched}
\title{
Final Matching
}
\description{
Compares the ntVLFs and aaVLFs with the same specimen identifier, and determines which ntVLFs are concordant with aaVLFs.
}
\usage{
overall.matched(positionNuc, positionAA, nuclength, aalength)
}
\arguments{
  \item{positionNuc}{
A list containing the names of the specimen and the ntVLF positions for specimens that have both aaVLFs and ntVLFs.
}
  \item{positionAA}{
A list containing the names of the specimen, the aaVLF, and the position of the aaVLF for specimens that have both aaVLFs and ntVLFs.
}
  \item{nuclength}{
The length of the nucleotide sequence (should by 3X the length of the amino acid sequence)
}
  \item{aalength}{
The length of the amino acid sequence (should be 1/3 the length of the nucleotide sequence)
}
}
\details{
The argument positionNuc can be calculated using the function nucleotide.matching.positions.  The argument positionAA can be calculated using the function aminoAcid.matching.positions.
}
\value{
A list of each ntVLF containing the specimen identifier in the first position of each list entry, the species name in the second position of each list entry, the aaVLF in the third position of each entry, the amino acid position of the aaVLF in the fourth entry, and the codon position of the concordant ntVLF in each position of the entry.  If multiple ntVLFs have concordance with one aaVLF, then that aaVLF may contain multiple entries in the list, one for each ntVLF.
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
birds_aa_speciesNames <- birds_aminoAcids[,2]
aminoAcids_specimenNumber <- nrow(birds_aminoAcids)
birds_aminoAcid_count <- aa.count.function(birds_aminoAcids, 216)
aa_freq.Mat <- aa.frequency.matrix.function(birds_aminoAcid_count, 216)
bird_aa_freq <- aa.specimen.frequencies(aa_freq.Mat, birds_aminoAcids, birds_aa_speciesNames, 216)
aminoAcid_Modal <- aa.MODE(aa_freq.Mat, 216)
birds_aminoAcid_specimenVLFcount <- aa.VLF.count.spec(bird_aa_freq, 0.001, 216)
birds_aaVLFconvert <- aa.VLF.convert.matrix(birds_aminoAcids, bird_aa_freq, 0.001, 216)
birds_aminoAcidVLFs <- VLF.aminoAcids(birds_aaVLFconvert, birds_aminoAcids, 216)
birds_aaVLFreduced <- aa.VLF.reduced(birds_aminoAcidVLFs, birds_aminoAcid_specimenVLFcount, 216)
birds_aaSpecies <- separate(birds_aaVLFreduced)
birds_aminoAcid_singleAndShared <- aa.find.singles(birds_aaSpecies, 216)

#Concordance Analysis
VLF_match <- find.matching(bird_VLFreduced, birds_aaVLFreduced, 648, 216)
position_matchingNuc <- nucleotide.matching.positions(VLF_match[[2]], 648)
position_matchingAA <- aminoAcid.matching.positions(VLF_match[[1]], 216)
matching_comparison <- overall.matched(position_matchingNuc, position_matchingAA, 648, 216)}
}