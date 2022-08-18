\name{VLF.aminoAcids}
\alias{VLF.aminoAcids}
\title{
Amino Acid VLFs}
\description{
Converts a matrix of amino acid frequencies into a matrix of amino acids.
}
\usage{
VLF.aminoAcids(convert.matrix, seq.matrix, seqlength)
}
\arguments{
  \item{convert.matrix}{
A matrix consisting of only aaVLF frequencies for each specimen, and NAs in every other position of the sequence.
}
  \item{seq.matrix}{
A matrix of amino acid sequences.
}
  \item{seqlength}{
The length of the amino acid sequence.
}
}
\details{
The argument convert.matrix can be calculated using the function aa.VLF.convert.matrix}
\value{
A matrix containing only aaVLFs and NAs in every other position of the sequence.
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
aminoAcid_Modal <- aa.MODE(aminoAcid_frequency.Matrix, 216)
birds_aminoAcid_specimenVLFcount <- aa.VLF.count.spec(bird_aminoAcid_frequencies, 0.001, 216)
birds_aaVLFconvert <- aa.VLF.convert.matrix(birds_aminoAcids, bird_aminoAcid_frequencies, 0.001, 
    216)
birds_aminoAcidVLFs <- VLF.aminoAcids(birds_aaVLFconvert, birds_aminoAcids, 216)}
}
