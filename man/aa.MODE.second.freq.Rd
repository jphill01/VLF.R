\name{aa.MODE.second.freq}
\alias{aa.MODE.second.freq}
\title{
Amino Acid Second Modal Frequency
}
\description{
Returns the frequencies of the amino acids that occur second most often in each position of a matrix of amino acid sequences.
}
\usage{
aa.MODE.second.freq(freq.matrix, seqlength)
}
\arguments{
  \item{freq.matrix}{
Frequency matrix for amino acids.
}
  \item{seqlength}{
Length of amino acid sequences.
}
}
\details{
The argument freq.matrix can be calculated using the function aa.frequency.matrix.function
}
\value{
A vector containing the frequencies of the second modal amino acid sequence.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{data(birds_aminoAcids)
birds_aminoAcid_count <- aa.count.function(birds_aminoAcids, 216)
aminoAcid_frequency.Matrix <- aa.frequency.matrix.function(birds_aminoAcid_count, 216)
aminoAcid_secondModalFreq <- aa.MODE.second.freq(aminoAcid_frequency.Matrix, 216)}
}