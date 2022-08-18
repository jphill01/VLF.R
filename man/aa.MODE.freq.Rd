\name{aa.MODE.freq}
\alias{aa.MODE.freq}
\title{
Amino Acid Modal Frequencies
}
\description{
Returns the frequencies of the amino acids that occur most often in each position of the sequence.
}
\usage{
aa.MODE.freq(freq.matrix, seqlength)
}
\arguments{
  \item{freq.matrix}{
Frequency matrix for amino acids.
}
  \item{seqlength}{
Length of the amino acid sequences
}
}
\details{
The argument freq.matrix can be calculated using the function aa.frequency.matrix.function
}
\value{
A vector of frequencies for the first modal sequence.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{data(birds_aminoAcids)
birds_aminoAcid_count <- aa.count.function(birds_aminoAcids, 216)
aminoAcid_frequency.Matrix <- aa.frequency.matrix.function(birds_aminoAcid_count, 216)
aminoAcid_firstModalFreq <- aa.MODE.freq(aminoAcid_frequency.Matrix, 216)}
}