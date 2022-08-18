\name{aa.find.singles}
\alias{aa.find.singles}
\title{
Find amino acid singles}
\description{
Determines the number of shared and singleton amino acid VLFs.}
\usage{
aa.find.singles(aaSpecies, seqlength)
}
\arguments{
  \item{aaSpecies}{
List of amino acid sequences separated by species name.}
  \item{seqlength}{
Length of amino acid sequences.
}
}
\details{
The argument aaSpecies contains only amino acid VLFs, and NAs in any other position in the sequence.  The list can be created using the separate function.}
\value{
A matrix containing the number of singleton and shared aaVLFs in each position of the barcode.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
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
