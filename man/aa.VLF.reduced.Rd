\name{aa.VLF.reduced}
\alias{aa.VLF.reduced}
\title{
Amino Acid Reduced
}
\description{
Reduces a matrix of amino acid very low frequency variants (aaVLFs) so that only those specimen that contain aaVLFs remain
}
\usage{
aa.VLF.reduced(NA.matrix, sCount, seqlength)
}
\arguments{
  \item{NA.matrix}{
A matrix with values for amino acid very low frequency variants for each specimen and NAs in all other positions.
}
  \item{sCount}{
A vector for amino acid very low frequency variant (VLF) counts for each specimen in the NA.matrix.}
  \item{seqlength}{
Length of the amino acid sequences.
}
}
\details{
The argument NA.matrix can be calculated using aa.VLF.convert.matrix and VLF.aminoAcids, the sCount argument can be calculaed using the aa.VLF.count.spec function.
}
\value{
A matrix containing only specimen with aaVLFs, and only the aaVLF values in the sequences.  All other positions of the sequence contain NAs.
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
birds_aaVLFreduced <- aa.VLF.reduced(birds_aminoAcidVLFs, birds_aminoAcid_specimenVLFcount, 216)}
}
