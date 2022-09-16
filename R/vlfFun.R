vlfFun <-
function(x, p=0.001, seqlength=648, own = NULL){
	species.names <- x[,2]

	specimen.Number <- nrow(x)

	rownames(x) <- species.names

	Nuc.count <- count.function(x, specimen.Number,seqlength)

	frequency.matrix <- ffrequency.matrix.function(Nuc.count,seqlength)

	spec.freq <- specimen.frequencies(frequency.matrix, x, specimen.Number, species.names,seqlength)

	nucleotide.modalSequence <- MODE(frequency.matrix,seqlength)

	first.modal.frequencies <- MODE.freq(frequency.matrix,seqlength)

	second.modal.frequencies <- MODE.second.freq(frequency.matrix,seqlength)

	First_conserved_100 <- conservation_first(first.modal.frequencies, 1,seqlength)

	First_conserved_99.9 <- conservation_first(first.modal.frequencies, (1-p),seqlength)

	FirstAndSecond_conserved_99.9 <- conservation_two(first.modal.frequencies, second.modal.frequencies, (1-p),seqlength)

	specimen_VLFcount <- VLF.count.spec(spec.freq, p,seqlength)

	position_VLFcount <- VLF.count.pos(spec.freq, p,seqlength)

	VLFconvert <- VLF.convert.matrix(x, spec.freq, p,seqlength)

	VLFnuc <- VLF.nucleotides(VLFconvert, x,seqlength)

	VLFreduced <- VLF.reduced(VLFnuc, specimen_VLFcount, seqlength)

	species <- separate(VLFreduced)

	singleAndShared <- find.singles(species,seqlength)
	
	if(is.null(own)){
		foo<-list(modal=nucleotide.modalSequence, con100=First_conserved_100, conp=First_conserved_99.9,combine=FirstAndSecond_conserved_99.9,specimen=specimen_VLFcount, position=position_VLFcount, sas=singleAndShared, VLFmatrix = VLFreduced)
		class(foo)<-"vlf"
		foo
	}
	else{
		ownspec.freq <- specimen.frequencies(frequency.matrix, own, nrow(own), own[,2], seqlength)
		ownspec.VLFcount <- VLF.count.spec(ownspec.freq, p, seqlength)
		ownpos.VLFcount <- VLF.count.pos(ownspec.freq, p, seqlength)
		own.VLFconvert <- VLF.convert.matrix(own, ownspec.freq, p, seqlength)
		own.VLFnuc <- VLF.nucleotides(own.VLFconvert, own, seqlength)
		own.VLFreduced <- VLF.reduced(own.VLFnuc, ownspec.VLFcount, seqlength)
		
		foo<-list(modal=nucleotide.modalSequence, con100=First_conserved_100, conp=First_conserved_99.9,combine=FirstAndSecond_conserved_99.9,specimen=specimen_VLFcount, position=position_VLFcount, sas=singleAndShared, VLFmatrix = VLFreduced, ownSpecCount = ownspec.VLFcount, ownPosCount = ownpos.VLFcount, ownVLFMatrix = own.VLFnuc, ownVLFreduced = own.VLFreduced)
		class(foo)<-"vlf"
		foo
	}	
}

### Code altered by Jarrett when own is not NULL ###

vlfFun <- function (x, p = 0.001, seqlength = 648, own = NULL) 
{
  species.names <- x[, 2]
  specimen.Number <- nrow(x)
  rownames(x) <- species.names
  Nuc.count <- count.function(x, specimen.Number, seqlength)
  frequency.matrix <- ffrequency.matrix.function(Nuc.count, 
                                                 seqlength)
  spec.freq <- specimen.frequencies(frequency.matrix, x, specimen.Number, 
                                    species.names, seqlength)
  nucleotide.modalSequence <- MODE(frequency.matrix, seqlength)
  first.modal.frequencies <- MODE.freq(frequency.matrix, seqlength)
  second.modal.frequencies <- MODE.second.freq(frequency.matrix, 
                                               seqlength)
  First_conserved_100 <- conservation_first(first.modal.frequencies, 
                                            1, seqlength)
  First_conserved_99.9 <- conservation_first(first.modal.frequencies, 
                                             (1 - p), seqlength)
  FirstAndSecond_conserved_99.9 <- conservation_two(first.modal.frequencies, 
                                                    second.modal.frequencies, (1 - p), seqlength)
  specimen_VLFcount <- VLF.count.spec(spec.freq, p, seqlength)
  position_VLFcount <- VLF.count.pos(spec.freq, p, seqlength)
  VLFconvert <- VLF.convert.matrix(x, spec.freq, p, seqlength)
  VLFnuc <- VLF.nucleotides(VLFconvert, x, seqlength)
  VLFreduced <- VLF.reduced(VLFnuc, specimen_VLFcount, seqlength)
  species <- separate(VLFreduced)
  singleAndShared <- find.singles(species, seqlength)
  if (is.null(own)) {
    foo <- list(modal = nucleotide.modalSequence, con100 = First_conserved_100, 
                conp = First_conserved_99.9, combine = FirstAndSecond_conserved_99.9, 
                specimen = specimen_VLFcount, position = position_VLFcount, 
                sas = singleAndShared, VLFmatrix = VLFreduced)
    class(foo) <- "vlf"
    foo
  }
  else {
    own.species.names <- own[, 2]
    own.specimen.Number <- nrow(own)
    rownames(own) <- own.species.names
    own.Nuc.count <- count.function(own, own.specimen.Number, seqlength)
    own.frequency.matrix <- ffrequency.matrix.function(own.Nuc.count, 
                                                       seqlength)
    ownspec.freq <- specimen.frequencies(frequency.matrix, 
                                         own, nrow(own), own[, 2], seqlength)
    own.nucleotide.modalSequence <- MODE(own.frequency.matrix, seqlength)
    own.first.modal.frequencies <- MODE.freq(own.frequency.matrix, seqlength)
    own.second.modal.frequencies <- MODE.second.freq(own.frequency.matrix, 
                                                     seqlength)
    own.First_conserved_100 <- conservation_first(own.first.modal.frequencies, 
                                                  1, seqlength)
    own.First_conserved_99.9 <- conservation_first(own.first.modal.frequencies, 
                                                   (1 - p), seqlength)
    own.FirstAndSecond_conserved_99.9 <- conservation_two(own.first.modal.frequencies, 
                                                          own.second.modal.frequencies, (1 - p), seqlength)
    ownspec.VLFcount <- VLF.count.spec(ownspec.freq, p, seqlength)
    ownpos.VLFcount <- VLF.count.pos(ownspec.freq, p, seqlength)
    own.VLFconvert <- VLF.convert.matrix(own, ownspec.freq, 
                                         p, seqlength)
    own.VLFnuc <- VLF.nucleotides(own.VLFconvert, own, seqlength)
    own.VLFreduced <- VLF.reduced(own.VLFnuc, ownspec.VLFcount, 
                                  seqlength)
    own.species <- separate(own.VLFreduced)
    own.singleAndShared <- find.singles(own.species, seqlength)
    foo <- list(own.modal = own.nucleotide.modalSequence, own.con100 = own.First_conserved_100, 
                own.conp = own.First_conserved_99.9, own.combine = own.FirstAndSecond_conserved_99.9, 
                own.sas = own.singleAndShared, own.VLFmatrix = own.VLFreduced, ownSpecCount = ownspec.VLFcount, 
                ownPosCount = ownpos.VLFcount, ownVLFMatrix = own.VLFnuc, 
                ownVLFreduced = own.VLFreduced)
    class(foo) <- "vlf"
    foo
  }
}