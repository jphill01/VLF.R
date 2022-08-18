\name{concordanceFun}
\alias{concordanceFun}
\title{
VLF Concordance Check Function
}
\description{
Compares ntVLFs to aaVLFs to see if they are concordant (i.e., if the ntVLF causes the aaVLF).
}
\usage{
concordanceFun(nuc, aa, nuclength = 648, aalength = 216, aminoAcid_Modal)
}
\arguments{
  \item{nuc}{
A matrix of ntVLFs, that contains only those specimen with VLFs, and a sequence with only VLF nucleotides and NAs in all other positions of the nucleotide sequences.
}
  \item{aa}{
A matrix of aaVLFs, that contains only those specimen with VLFs, and a sequence with only VLF amino acids and NAs in all other positions of the amino acid sequence.
}
  \item{nuclength}{
The length of the nucleotide sequence.  By default is 648.
}
  \item{aalength}{
The length of the amino acid sequence.  By default is 216.
}
  \item{aminoAcid_Modal}{
The modal amino acid sequence (i.e., the amino acid sequence that occurs most often in the given sequences)
}
}
\details{
The argument nuc can be taken from the VLFmatrix output from the vlfFun function.  The argument aa can be taken from the VLFmatrix output from the aminoAcidFun function.  The argument aminoAcid_Modal can be taken from the modal output from the aminoAcidFun function.}
\value{
\item{matched}{A list of the concordant ntVLFs and aaVLFs.  Contains the specimen identifier, the species name, the concordant amino acid, the amino acid position, and the concordant amino acid position.  There may be mulitple entries for the same aaVLF if that VLF is concordant to more than one ntVLF.}
\item{codons}{A vector containing calculations for how many of the concordant amino acids were caused by changes in each of the nucleotide codon positions.}
\item{concordantType}{Contains information on how many of the concordant aaVLFs had a change in amino acid residue type and how many remained in the same amino acid residue category.}
\item{aminoAcidType}{Contains information on how many of the aaVLFs had a change in amino acid residue type and how many remained in the same amino acid residue category.}
\item{concordNuc}{Gives the number of ntVLFs that showed concordance to aaVLFs.}
\item{concordAA}{Gives the number of aaVLFs that showed concordance to ntVLFs.}
\item{sequences}{Gives the number of sequences that had both ntVLFs and aaVLFs.}
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}
\examples{
\dontrun{#VLF analysis
data(birds)
bird_vlfAnalysis <- vlfFun(birds)

#Amino Acid analysis
data(birds_aminoAcids)
bird_aaAnalysis <- aminoAcidFun(birds_aminoAcids)

#Concordance analysis
bird_Concordance <- concordanceFun(bird_vlfAnalysis$VLFmatrix, bird_aaAnalysis$VLFmatrix, 648, 216,
    bird_aaAnalysis$modal)}
}
