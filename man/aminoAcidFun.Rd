\name{aminoAcidFun}
\alias{aminoAcidFun}
\title{Amino Acid VLF Analysis Function
}
\description{
Runs the full amino acid VLF analysis for the user and outputs results.
}
\usage{
aminoAcidFun(x, p = 0.001, seqlength = 216, own = NULL)
}
\arguments{
  \item{x}{
A matrix of amino acid sequences with the first column containing the unique specimen identifier, the second column containing the species name and the remaining columns containing the amino acid sequence.
}
  \item{p}{
A VLF designation frequency cut-off to be used within the analysis.  By default p = 0.001.
}
  \item{seqlength}{
The length of the amino acid sequence.  By default seqlength = 216.
}
\item{own}{
If the user wants to compare their own sequences separate from reference sequences, then argument can be used.  Similar to x, this argument is a matrix of amino acid sequences with the first column containing the unique specimen identifier, the second column contains the species name and the remaining columns containing the nucleotide sequence.  By default own = NULL.)
}
}
\value{
\item{modal}{A vector containing the amino acid sequence that occurs most often in the dataset.}
\item{con100}{The number of amino acid positions that are 100\% conserved in the sequence}
\item{conp}{The number of amino acid positions that are (1-p)\% conserved in the sequence}
\item{combine}{The number of amino acid positions that are (1-p)\% conserved when combining the first and second modal sequences.}
\item{specimen}{A vector containing the number of VLFs for each specimen in the dataset.}
\item{position}{A vector containing the number of VLFs for each position in the sequences.}
\item{sas}{A matrix containing vectors of single and shared amino acid VLF counts for each position of the sequence.}
\item{VLFmatrix}{A matrix containing only those specimen that have VLFs as well as the amino acid at the positions that contain VLFs and NAs in all other positions.}
\item{ownSpecCount}{A vector containing the number of VLFs for each specimen in the users own specified dataset.  Only appears if own is not NULL.}
\item{ownPosCount}{A vector containing the number of VLFs for each position in the sequences of the users own specified dataset.  Only appears if own is not NULL.}
\item{ownVLFMatrix}{A matrix containing only those amino acids at the positions that contain VLFs and NAs in all other positions of the sequence.  Only appears if own is not NULL.}
\item{ownVLFreduced}{ownVLFMatrix}{A matrix containing only those specimen that have VLFs as well as the amino acids at the positions that contain VLFs and NAs in all other positions of the sequence.  Only appears if own is not NULL.}
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{data(birds_aminoAcids)
bird_aaAnalysis <- aminoAcidFun(birds_aminoAcids)}
}