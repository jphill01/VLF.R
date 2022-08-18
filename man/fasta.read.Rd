\name{fasta.read}
\alias{fasta.read}
\title{
Read Fasta Files}
\description{
Reads in fasta files and converts into sequence matrix
}
\usage{
fasta.read(file, seqlength = 648, pos1 = 1, pos2 = 3)
}
\arguments{
  \item{file}{
A fasta file to be read in.
}
  \item{seqlength}{
Length of sequence.
}
  \item{pos1}{
The position within the fasta title of the unique specimen identifier.  By default pos1 = 1.
}
  \item{pos2}{
The position within the fasta title of the species name.  By default pos2 = 3.
}
}
\value{
A matrix of sequences, with the unique specimen identifers in the first column, the species names in the second column, and the sequence starting in the third column.
}
\author{
Taryn B. T. Athey and Paul D. McNicholas
}