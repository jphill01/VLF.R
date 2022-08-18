\name{vlfFun}
\alias{vlfFun}
\title{
Nucleotide VLF Assessment Function}
\description{
Runs the full nucleotide VLF analysis for the user and outputs the results. 
}
\usage{
vlfFun(x, p = 0.001, seqlength = 648, own = NULL)
}
\arguments{
  \item{x}{
A matrix of nucleotide sequences with the first column containing the unique specimen identifier, the second column containing the species name and the remaining columns containing the nucleotide sequence.}
  \item{p}{
A VLF designation frequency cut-off to be used within the analysis.  By default p = 0.001.
}
  \item{seqlength}{
The length of the nucleotide sequence.  By default seqlength = 648.
}
\item{own}{
If the user wants to compare their own sequences separate from reference sequences, then argument can be used.  Similar to x, this argument is a matrix of nucleotide sequences with the first column containing the unique specimen identifier, the second column contains the species name and the remaining columns containing the nucleotide sequence.  By default own = NULL.)
}
}
\value{
\item{modal}{A vector containing the nucleotide sequence that occurs most often in the dataset.}
\item{con100}{The number of nucleotide positions that are 100\% conserved in the sequence, separated by codon position.}
\item{conp}{The number of nucleotide positions that are (1-p)\% conserved in the sequence, separated by codon position.}
\item{combine}{The number of nucleotide positions that are (1-p)\% conserved when combining the first and second modal sequences.}
\item{specimen}{A vector containing the number of VLFs for each specimen in the dataset.}
\item{position}{A vector containing the number of VLFs for each position in the sequences.}
\item{sas}{A matrix containing vectors of single and shared ntVLF counts for each position in the sequences.}
\item{VLFmatrix}{A matrix containing only those specimen that have VLFs as well as the nucleotides at the positions that contain VLFs and NAs in all other positions of the sequence.}
\item{ownSpecCount}{A vector containing the number of VLFs for each specimen in the users own specified dataset.  Only appears if own is not NULL.}
\item{ownPosCount}{A vector containing the number of VLFs for each position in the sequences of the users own specified dataset.  Only appears if own is not NULL.}
\item{ownVLFMatrix}{A matrix containing only those nucleotides at the positions that contain VLFs and NAs in all other positions of the sequence.  Only appears if own is not NULL.}
\item{ownVLFreduced}{ownVLFMatrix}{A matrix containing only those specimen that have VLFs as well as the nucleotides at the positions that contain VLFs and NAs in all other positions of the sequence.  Only appears if own is not NULL.}
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{data(birds)
bird_vlfAnalysis <- vlfFun(birds)}
}