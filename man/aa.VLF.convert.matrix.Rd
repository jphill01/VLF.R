\name{aa.VLF.convert.matrix}
\alias{aa.VLF.convert.matrix}
\title{
Convert Amino Acid Matrix
}
\description{
Converts a matrix of amino acid frequencies for each specimen into a matrix consisting of only VLF values and NAs in every non-VLF position.
}
\usage{
aa.VLF.convert.matrix(seq.matrix, freq, p, seqlength)
}
\arguments{
  \item{seq.matrix}{
A matrix of aligned DNA barcode amino acid sequences.
}
  \item{freq}{
A matrix of amino acid frequencies for each specimen.
}
  \item{p}{
A very low frequency variant cut-off frequency.  Any frequency in the freq matrix below this value is considered to be a very low frequency variant.
}
  \item{seqlength}{
The length of the amino acid sequences.}
}
\value{
A matrix of VLF amino acid frequencies, containing only those nucleotide frequencies that occur less than the designation p value, and NAs in each other position of the matrix.
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
birds_aaVLFconvert <- aa.VLF.convert.matrix(birds_aminoAcids, bird_aminoAcid_frequencies, 0.001, 
    216)}
}
