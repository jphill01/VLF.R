\name{aaVLFs.to.modalchanges}
\alias{aaVLFs.to.modalchanges}
\title{
Amino Acid Changes
}
\description{
Determines how many aaVLFs have changed ``type'' of amino acid from the modal amino acid sequence.  Amino acid types are polar charged, polar uncharged, non-polar, and those with a unique side group.}
\usage{
aaVLFs.to.modalchanges(modal, AminoAcidList, aalength)
}
\arguments{
  \item{modal}{
The modal amino acid sequence (i.e., the amino acid sequence that occurs most often based on the amino acid frequency matrix)}
  \item{AminoAcidList}{
Matrix of VLF amino acid sequences containing only aaVLFs and NAs anywhere else}
  \item{aalength}{
Amino Acid sequence length.}
}
\details{
The argument modal can be created using the MODE function.  The argument AminoAcidList can be created using the aa.VLF.convert.matrix, VLF.aminoAcids, and aa.VLF.reduced functions.}
\value{
A sameAll value representative of the number of amino acids that were the same type as the modal, a changedAll value representative of the number of amino acids that changed amino acid type from the modal.
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
birds_aaVLFconvert <- aa.VLF.convert.matrix(birds_aminoAcids, bird_aminoAcid_frequencies, 0.001, 
    216)
birds_aminoAcidVLFs <- VLF.aminoAcids(birds_aaVLFconvert, birds_aminoAcids, 216)
birds_aaVLFreduced <- aa.VLF.reduced(birds_aminoAcidVLFs, birds_aminoAcid_specimenVLFcount, 216)
All_aaType_change <- aaVLFs.to.modalchanges(aminoAcid_Modal, birds_aaVLFreduced, 216)}
}
