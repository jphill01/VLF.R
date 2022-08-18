\name{VLF.convert.matrix}
\alias{VLF.convert.matrix}
\title{
VLF Matrix Convert
}
\description{
Converts a matrix of nucleotide frequencies for each specimen into a matrix consisting entirely of very low frequency variant (VLF) frequencies and NAs in each other position.
}
\usage{
VLF.convert.matrix(seq.matrix, freq, p, seqlength)
}
\arguments{
  \item{seq.matrix}{
A matrix of aligned DNA barcode sequences.
}
  \item{freq}{
A matrix of nucleotide frequencies for each specimen.
}
  \item{p}{
A very low frequency variant designation cut off frequency.  Any frequency in the freq matrix below this value is considered to be a very low frequency variant.
}
  \item{seqlength}{
Length of nucleotide sequence.}
}
\details{
The argument freq can be calculated using the function specimen.frequencies.
}
\value{
A matrix of VLF nucleotide frequencies, containing only those nucleotide frequencies that occur less than the designation p value, and NAs in each other position of the matrix.
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
birds_aaVLFconvert <- aa.VLF.convert.matrix(birds_aminoAcids, bird_aminoAcid_frequencies,0.001,216)}
}
