\name{VLF-package}
\alias{VLF-package}
\alias{VLF}
\docType{package}
\title{
Frequency Matrix Approach for Assessing Very Low Frequency Variants in Sequence Records
}
\description{
Using frequency matrices, very low frequency variants (VLFs) are assessed for amino acid and nucleotide sequences.  The VLFs are then compared to see if they occur in only one member of a species, singleton VLFs, or if they occur in multiple members of a species, shared VLFs.  The amino acid and nucleotide VLFs are then compared to see if they are concordant with one another.  Amino acid VLFs are also assessed to determine if they lead to a change in amino acid residue type, and potential changes to protein structures.
}
\details{
vlfFun()
aminoAcidFun()
concordanceFun()
}
\author{
Taryn B. T. Athey and Paul D. McNicholas

Maintainer: Taryn B. T. Athey <tathey@uoguelph.ca>
}
\keyword{VLFs}
\examples{
\dontrun{#VLF analysis
data(birds)
bird_vlfAnalysis <- vlfFun(birds)

#Amino Acid analysis
data(birds_aminoAcids)
bird_aaAnalysis <- aminoAcidFun(birds_aminoAcids)

#Concordance analysis
nuc_matrix <- bird_vlfAnalysis$VLFmatrix
aa_matrix <- bird_aaAnalysis$VLFmatrix
aa_modal <- bird_aaAnalysis$modal
bird_Concordance <- concordanceFun(nuc_matrix, aa_matrix, 648, 216, aa_modal)}}
