\name{aa.specimen.frequencies}
\alias{aa.specimen.frequencies}
\title{
Specimen Nucleotide Frequencies
}
\description{
Converts a matrix of amino acid sequences into a matrix of amino acid frequencies.
}
\usage{
aa.specimen.frequencies(freq, seq.matrix, spec.names, seqlength)
}
\arguments{
  \item{freq}{
Frequency matrix for amino acids.
}
  \item{seq.matrix}{
Matrix of specimen amino acid sequences.
}
  \item{spec.names}{
A vector of the species names for each specimen in aminoAcids in the ordfer they appear in the matrix.
}
  \item{seqlength}{
Length of amino acid sequences.
}
}
\details{
The argument freq can be calculated using the function aa.frequency.matrix.function.
}
\value{
A matrix containing the frequencies of each amino acid in the sequence.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{data(birds_aminoAcids)
birds_aminoAcid_speciesNames <- birds_aminoAcids[,2]
aminoAcids_specimenNumber <- nrow(birds_aminoAcids)
birds_aminoAcid_count <- aa.count.function(birds_aminoAcids,216)
aminoAcid_frequency.Matrix <- aa.frequency.matrix.function(birds_aminoAcid_count,216)
bird_aminoAcid_frequencies <- aa.specimen.frequencies(aminoAcid_frequency.Matrix, birds_aminoAcids,
    birds_aminoAcid_speciesNames, 216)}
}
