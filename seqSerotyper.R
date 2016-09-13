# File seqSerotyper.R
# Author: Andries J. van Tonder 
#
# Copyright (C) 2016 University of Oxford
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
#  This is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this software package; if not, write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA  02110-1301  USA

###############################
## Output of seqSerotyper -help
###############################

help <- paste(
  "",
  "seqSerotyper.R assigns serotypes to pneumococcal WGS",
  "Author: Andries J. van Tonder",
  "",
  "Usage: Rscript seqSerotyper.R -dataFile -srcDir",
  "	e.g. Rscript seqSerotyper.R -dataFile pneumo_fasta_ids.txt -srcDir /home/user/Software/avt_scripts/seqSerotyper/bin",
  "",
  "Options:",
  "",
  "-dataFile",
  " The path of a file containing sample ids and file paths to fasta files",
  "-srcDir",
  " Path of the directory that contains the reference files and related scripts",
  "",
  sep="\n"
)

##########################################################################################
## Create tmp directory to work in and save main seqSerotyper directory to seqSerotyperDir
##########################################################################################

system2("mkdir", "tmp")
seqSerotyperDir <- getwd()

######################  
## Load seqinr library
######################

library(seqinr)

###################
## Define functions
###################

####################################################
## Format the dir appropriately.
## @dir: The directory path to be formated properly.
####################################################
formatDir = function(dir = NULL){
  dirLength = nchar(dir)
  if(substr(dir, dirLength, dirLength) !='/'){
    dir = paste(dir,"/",sep="")
  }
  return(dir)
}

################################################
## checkExist: checks file paths exist
## @filePath: the path of the file to be checked
################################################

checkExist <- function(filePath = NULL){
  message(filePath)
  message(class(filePath))
	doesNotExist = which(!file.exists(filePath))
  message("Done!")
	if(length(doesNotExist) > 0){
		stop(paste("The file", filePath[doesNotExist]," does not exist!", collapse="\n", sep=""))
	}
}

####################################################
## Turn command line inputs into a n x 2 matrix
## @args: A string vector of command line arguments.
####################################################
getCommandLineInputMatrix = function(args = NULL){
  argsCount = length(args)/2
  inputs = matrix(nrow=argsCount, ncol=2)
  inputs[,1] = args[c(1:argsCount)*2 -1]
  
  if(length(which(regexpr("-",inputs[,1]) != 1)) > 0){
    stop("Argument names must start with '-'! E.g. -dataFile.")
  }
  
  inputs[,1] = gsub("-", "", inputs[,1])
  inputs[,2] = args[c(1:argsCount)*2]
  return(inputs)
  
}

#####################################################################################
## Extract input arguments.
## @argName: Name of the argument to be extracted
## @commandLineInputs: A n x 2 matrix, where n is the number of command line inputs.
## @default: default value of the input argument if not provided
## The first column of the matrix contains the names of the arguments, and the second
## column contains the arguments values.
#####################################################################################
extractInputArgument = function(argName = NULL, commandLineInputs = NULL,
                                default = NULL, canBeNULL = FALSE){
  
  argIndex = which(commandLineInputs[,1] == argName)
  argIndexCount = length(argIndex)
  if(argIndexCount == 0){
    if(canBeNULL | !is.null(default)){
      return(default)
    }else{
      stop(paste(c("The argument", argName, "must be specified!"), collapse=" "))
    }
  }else if(argIndexCount > 1){
    stop(paste(c("The argument", argName, "has been specified multiple times!"), collapse=" "))
  }else{
    return(commandLineInputs[argIndex,2])
  }
  
}

###########################################################################
## RunBlat: runs blat using appropriate query and reference sequences
## @querySeqs: sequence/s to be compared against reference sequence
## @refNuc: reference sequence against which query sequences compared
## @blatOutput: name/s of blat output files in standard blast output format
###########################################################################

RunBlat <- function(querySeqs, refNuc, blatOutput) {
  system2(paste(blatPath),args=c("-t=dna",querySeqs,"-q=dna", refNuc, blatOutput,"out=blast"))
}

###########################################################################
## RunBlatBlast8: runs blat using appropriate query and reference sequences
## @querySeqs: sequence/s to be compared against reference sequence
## @refNuc: reference sequence against which query sequences compared
## @blatOutput: name/s of blat output files in blast8 output format
###########################################################################

RunBlatBlast8 <- function(querySeqs, refNuc, blatOutput) {
  system2(paste(blatPath),args=c("-t=dna",querySeqs,"-q=dna", refNuc, blatOutput,"out=blast8"))
}

#############################################################################
## RunBlastall: runs blastall using appropriate query and reference sequences
## @BARef: reference sequence against which query sequences compared
## @BAquerySeqs: sequence/s to be compared against reference sequence
## @blastAllOut: names of blastall output files
#############################################################################
  
RunBlastall <- function(BARef, BAquerySeqs, blastAllOut) {
  system2(paste(blastallPath),args=c("-p blastn", "-d", BARef, "-i", BAquerySeqs,"-m 8","-o",blastAllOut))
}

##########################################################################
## runMview: runs mview which extracts sequence from blast output as fasta
## @mviewIn: input file in blast format
## @mviewOut: output file in pearson (fasta) format
##########################################################################

runMview <- function(mviewIn, mviewOut) {
  system2(paste(mviewPath),args=c("-in","blast",mviewIn,"-out","pearson", "-top", "1", ">",mviewOut))
}

#############################################################################
## runSplitMultifasta: splits multifasta file into individual fasta sequences
## @multifasta: multifasta file to be split
## @outDir: output directory for individual fasta files
#############################################################################

runSplitMultifasta <- function(multifasta, outDir) {
  system2("perl",args=c(splitMultifastapath,"-i", multifasta,"-o",outDir))
}

##############################################################################
## runTranseq: runs transeq which translates nucleotide sequence to amino acid
## @TranseqIn: fasta nucleotide sequences
## @TranseqOut: fasta amino acid sequences
##############################################################################

runTranseq <- function(TranseqIn, TranseqOut) {
  system2(paste(transeqPath),args=c("-sequence", TranseqIn, "-outseq", TranseqOut))
}

###############################################
##  Read in the arguments from the command line
###############################################

args <- commandArgs(trailingOnly = TRUE)
if(length(args!=0)){
  if(args[1]=="-help" | args[1]=="-h"){
    cat(help,sep="\n")
    q("no")
  }
}

if((length(args)%%2)!=0 | length(args)==0) {
  cat(help,sep="\n")
  stop("\nIncorrect usage\n")
}

####################################
##  Get inputs from the command line
####################################
 
message("Arg length:")
message(length(args))
inputs <- getCommandLineInputMatrix(args = args)
message(nrow(inputs))
message(ncol(inputs))

############################
##  Input arguments required
############################

dataFilePath <- extractInputArgument(argName = "dataFile", commandLineInputs = inputs)
srcDir <- formatDir(extractInputArgument(argName = "srcDir", commandLineInputs = inputs, default = getwd()))
message("Input arguments have been read.")

###########################
## Reference files location
###########################

cat(paste0("Reference directory: ",srcDir))

#######################################
## Get the paths of the reference files
#######################################

assign11Eref <- paste(srcDir, "11A_wcjE_nuc.fas", sep="")
assign11ADref <- paste(srcDir, "11A_wcrL_nuc.fas", sep="")
assign18BCref <- paste(srcDir, "18B_18C_wciX_nuc.fas", sep="")
assign20ABref <- paste(srcDir, "20B_whaF_nuc.fas", sep="")
assign22AFref <- paste(srcDir, "22F_wcwA_nuc.fas", sep="")
assign25AFref <- paste(srcDir, "25F_wzg_nuc.fas", sep="")
assign33AFref <- paste(srcDir, "33A_33F_wcjE_nuc.fas", sep="")
assign37ref <- paste(srcDir, "37_tts_nuc.fas", sep="")
assign6ABref <- paste(srcDir, "6A_6B_wciP_nuc.fas", sep="")
assign6CEGref <- paste(srcDir, "6A_6B_wzy_nuc.fas", sep="")
assign6Fref <- paste(srcDir, "6A_wciN_alpha_nuc.fas", sep="")
assign6CDref <- paste(srcDir, "6D_wciN_beta_nuc.fas", sep="")
assign7AFref <- paste(srcDir, "7A_7F_wcwD_nuc.fas", sep="")
assign9AVref <- paste(srcDir, "9A_9V_wcjE_nuc.fas", sep="")
mainRefPath <- paste(srcDir, "Pneumo_serotypes.fasta", sep="")
splitMultifastapath <- paste(srcDir, "split_multifasta.pl", sep="")
blatPath <- paste(srcDir, "blat", sep="")
blastallPath <- paste(srcDir, "blastall", sep="")
transeqPath <- paste(srcDir, "transeq", sep="")
mviewPath <- paste(srcDir, "mview", sep="")

##############################
## Check reference files exist
###############################

checkExist(c(assign11Eref, assign11ADref, assign18BCref, assign20ABref, assign22AFref, assign25AFref, assign33AFref,
          assign37ref, assign6ABref, assign6CEGref, assign6CEGref, assign6Fref, assign6CDref, assign7AFref, assign9AVref,
          mainRefPath, splitMultifastapath, blatPath, blastallPath, transeqPath, mviewPath))

############################################################
## Read in data file and sort alphabetically based on the id
############################################################

data.df <- read.table(file=dataFilePath, header=T, as.is=T)
orderdata.df <- data.df[order(data.df$id),]
fastaFile <- orderdata.df$filePath
fastaID <- orderdata.df$id
numberSeqs <- length(fastaFile)
message(paste(numberSeqs, "sequence files found.", "\n"))

###############################
## Set working directory to tmp
###############################

setwd("tmp")

######################################################################
## Use blat to compare each genome sequence to serotype reference file
######################################################################


fastaFileTab <- paste(fastaID,".tab",sep="")

for (i in 1:length(fastaID)) {
  RunBlatBlast8(querySeqs = fastaFile[[i]], refNuc = mainRefPath, blatOutput = fastaFileTab[[i]])  
}

#############################
##  Process blat output files
#############################

blatTabOut <- NULL
for (i in 1:length(fastaFileTab)) blatTabOut[[i]] <- read.delim(fastaFileTab[i],header=F,sep="\t")

sortBlatTabOut <- NULL
for (i in 1:length(blatTabOut)) sortBlatTabOut[[i]] <- blatTabOut[[i]][order(-blatTabOut[[i]]$V12),]

topBlatHit <- NULL
for (i in 1:length(sortBlatTabOut)) topBlatHit[[i]] <- sortBlatTabOut[[i]][1,]

######################################
## Write out best hits for each genome
######################################

topBlatHitOut <- file("topBlatHit.out",open="a")
for (i in 1:length(topBlatHit)) {
  write.table(topBlatHit[[i]],file=topBlatHitOut,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
}
close(topBlatHitOut)

###########################################################
## Read best hits file and process by serotype and best hit
###########################################################

blatResults <- read.delim("topBlatHit.out",header=F,sep="\t")
addtopBlatHitData <- cbind(orderdata.df,blatResults)
fastablatID <- addtopBlatHitData$id

####################################################################
## Assign NT if length of alignment is <3000 i.e. no hit to database
####################################################################

NTresults <- subset(addtopBlatHitData,V4<3000)
nonNTresults <- subset(addtopBlatHitData,V4>3000)

if (nrow(NTresults) > 0) {
  NTresults$V13 <- "NT"
  NTout <- NTresults[,c("id","V13")]
} else if (nrow(NTresults)==0) {
  cat("No NT serotypes found","\n")
}

##########################################
## Extract serotypes for separate analyses
##########################################

sero6Results<-subset(nonNTresults, V1=="6A"|V1=="6B"|V1=="6Bii(6E)"|V1=="6F"|V1=="6G"|V1=="6C"|V1=="6D")
sero9Results<-subset(nonNTresults, V1=="9V"|V1=="9A")
sero7Results<-subset(nonNTresults, V1=="7A"|V1=="7F")
sero18Results<-subset(nonNTresults, V1=="18B"|V1=="18C")
sero33Results<-subset(nonNTresults, V1=="33A"|V1=="33F")
sero22Results<-subset(nonNTresults, V1=="22A"|V1=="22F")
sero11Results<-subset(nonNTresults, V1=="11A"|V1=="11D")
sero20Results<-subset(nonNTresults, V1=="20")
sero25Results<-subset(nonNTresults, V1=="25A"|V1=="25F")

###########################################################################
## Save remaining serotype predictions to own variable for saving out later
###########################################################################

resttoBlat <- nonNTresults[nonNTresults$V1!="6A"&nonNTresults$V1!="6B"&nonNTresults$V1!="6Bii(6E)"&nonNTresults$V1!="6F"&nonNTresults$V1!="6G"&nonNTresults$V1!="6C"&nonNTresults$V1!="6D"&nonNTresults$V1!="9V"&nonNTresults$V1!="9A"&nonNTresults$V1!="7A"&nonNTresults$V1!="7F"&nonNTresults$V1!="18B"&nonNTresults$V1!="18C"&nonNTresults$V1!="33A"&nonNTresults$V1!="33F"&nonNTresults$V1!="22A"&nonNTresults$V1!="22F"&nonNTresults$V1!="11A"&nonNTresults$V1!="11D"&nonNTresults$V1!="20"&nonNTresults$V1!="25A"&nonNTresults$V1!="25F",]
resttoBlatOut <- resttoBlat[,c("id","V1")]

#############################
## Serogroup 6 identification
#############################

if (nrow(sero6Results)>0) {
  
  fasta6 <- sero6Results$filePath
  fasta6ID <- sero6Results$id
  id6 <- sero6Results[,c("id","filePath")]
  tab6 <- paste(fasta6ID,".6all.tab",sep="")
  n6 <- length(fasta6)

###############################  
## Identify serotypes 6C and 6D
###############################

  for (i in 1:n6) {
    RunBlastall(BARef = assign6CDref, BAquerySeqs = fasta6[[i]], blastAllOut = tab6[[i]])    
  }
  
  out6CDin <- dir(pattern=".6all.tab")

  out6CDRead <- NULL
  for (i in 1:length(out6CDin)) out6CDRead[[i]] <- read.delim(out6CDin[i],header=F, sep="\t")

  sort6CDReadOut <- NULL
  for (i in 1:length(out6CDRead)) sort6CDReadOut[[i]] <- out6CDRead[[i]][order(-out6CDRead[[i]]$V12),]

  head6CDwciN <- NULL
  for (i in 1:length(sort6CDReadOut)) head6CDwciN[[i]] <- sort6CDReadOut[[i]][1,]

  head6CDwciNout <- file("6CD_head_out", open="a")
  for (i in 1:length(head6CDwciN)) {
    write.table(head6CDwciN[[i]],file=head6CDwciNout,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  
  close(head6CDwciNout)
 
  sero6CDresults <- read.delim("6CD_head_out",header=F,sep="\t")
  addSero6CDdata <- cbind(id6,sero6CDresults)
  assignOtherSero6 <- subset(addSero6CDdata, V4<500)
  assignSero6CD <- subset(addSero6CDdata, V4>500)
  
  if (nrow(assignSero6CD) > 0) {
  
    checkSero6CD <- assignSero6CD$filePath
    checkSero6CDid <- assignSero6CD$id
    checkSero6CDTab <- paste(checkSero6CDid,".6CD.tab",sep="")
    
    for (i in 1:length(checkSero6CD)) {
      RunBlat(querySeqs = checkSero6CD[[i]], refNuc = assign6CEGref, blatOutput = checkSero6CDTab[[i]])  
    }
    
    listSero6CDTab<-dir(pattern=".6CD.tab")
    listSero6CDTabFasta <- paste(listSero6CDTab,".fasta",sep="")
    
    for (i in 1:length(listSero6CDTab)) {
      runMview(mviewIn = listSero6CDTab[[i]], mviewOut = listSero6CDTabFasta[[i]])
    }
    
    listSero6CDFasta <- dir(pattern=".6CD.tab.fasta")
    system2("mkdir","6CD")

    for (i in 1:length(listSero6CDFasta)) {
      runSplitMultifasta(multifasta = listSero6CDFasta[[i]], outDir = "6CD")       
    }

    system2("rm","6CD/wzy.fsa")
    fsaSero6CDListT <- list.files("6CD",full.names=TRUE)
    fsaSero6CDListF <- list.files("6CD",full.names=FALSE)
    pepSero6CD <- paste(fsaSero6CDListF,".pep",sep="")

    for (i in 1:length(fsaSero6CDListT)){
      runTranseq(TranseqIn = fsaSero6CDListT[[i]], TranseqOut = pepSero6CD[[i]])
    }

    listSero6CDPep <- dir(pattern=".pep")
    readlistSero6CDPep <- NULL

    for (i in 1:length(listSero6CDPep)) {
      readlistSero6CDPep[[i]] <-read.fasta(listSero6CDPep[[i]],seqtype="AA",set.attributes=FALSE)
    }

    wzy117 <- NULL
    sero6CDOut<-NULL

    for (i in 1:length(readlistSero6CDPep)){
      wzy117[[i]] <- getFrag(readlistSero6CDPep[[i]],begin = 117, end=117)
    }

    for (i in 1:length(readlistSero6CDPep)){
      sero6CDOut[[i]]<-cbind(listSero6CDPep[[i]],wzy117[[i]])
    }

    sero6CDOutWrite <- file("6CD_out",open="a")
    for (i in 1:length(sero6CDOut)){
      write.table(sero6CDOut[[i]],file=sero6CDOutWrite,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
    close(sero6CDOutWrite)

    assignSero6CDin <- read.delim("6CD_out",header=F,sep="\t",colClasses="character")
    for (i in 1:length(assignSero6CDin$V2)){
      if (assignSero6CDin$V2[i]=="T")
        {assignSero6CDin$V3[i]="6C"}
      else if (assignSero6CDin$V2[i]=="I")
        {assignSero6CDin$V3[i]="6D"}
    }

  addIDassignSero6CDin <- cbind(checkSero6CDid,assignSero6CDin)
  assignSero6CDfinal <- addIDassignSero6CDin[,c("checkSero6CDid","V3")]

  } else if (nrow(assignSero6CD)==0) {
    cat("No serotypes 6C or 6D found","\n")
  }

###############################
## Identify serotypes 6E and 6G
###############################
  
  if (nrow(assignOtherSero6) > 0) {
        
    checknotSero6CD <- assignOtherSero6$filePath
    checknotSero6CDid <- assignOtherSero6$id
    checknotSero6CDTab <- paste(checknotSero6CDid,".6AG.tab",sep="")

    for (i in 1:length(checknotSero6CD)) {
      RunBlat(querySeqs = checknotSero6CD[[i]], refNuc = assign6CEGref, blatOutput = checknotSero6CDTab[[i]])
    }

    listSero6AGTab<-dir(pattern=".6AG.tab")
    listSero6AGTabFasta <- paste(listSero6AGTab,".fasta",sep="")

    for (i in 1:length(listSero6AGTab)) {
      runMview(mviewIn = listSero6AGTab[[i]], mviewOut = listSero6AGTabFasta[[i]])      
    }

    listSero6AGFasta <- dir(pattern=".6AG.tab.fasta")
    system2("mkdir","6AG")

    for (i in 1:length(listSero6AGFasta)) {
      runSplitMultifasta(multifasta = listSero6AGFasta[[i]], outDir = "6AG") 
    }

    system2("rm","6AG/wzy.fsa")
    fsaSero6AGListT <- list.files("6AG",full.names=TRUE)
    fsaSero6AGListF <- list.files("6AG",full.names=FALSE)
    pepSero6AG <- paste(fsaSero6AGListF,".6AG.pep",sep="")

    for (i in 1:length(fsaSero6AGListT)){
      runTranseq(TranseqIn = fsaSero6AGListT[[i]], TranseqOut = pepSero6AG[[i]]) 
    }

    listSero6AGpep <- dir(pattern=".6AG.pep")
    readlistSero6AGpep <- NULL

    for (i in 1:length(listSero6AGpep)) {
      readlistSero6AGpep[[i]] <-read.fasta(listSero6AGpep[[i]],seqtype="AA",set.attributes=FALSE)
    }

    wzy220 <- NULL
    wzy225 <- NULL
    sero6AGOut<-NULL

    for (i in 1:length(readlistSero6AGpep)){
      wzy220[[i]] <- getFrag(readlistSero6AGpep[[i]],begin = 220, end=220)
    }

    for (i in 1:length(readlistSero6AGpep)){
      wzy225[[i]] <- getFrag(readlistSero6AGpep[[i]],begin = 225, end=225)
    }

    sero6AGOut <- NULL
    for (i in 1:length(listSero6AGpep)) {
      sero6AGOut[[i]]<-cbind(listSero6AGpep[[i]],wzy220[[i]],wzy225[[i]])
    }

    sero6AGOutWrite <- file("6AG_out",open="a")
    for (i in 1:length(sero6AGOut)){
      write.table(sero6AGOut[[i]],file=sero6AGOutWrite,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
    close(sero6AGOutWrite)

    assignSero6E <- read.delim("6AG_out",header=F,sep="\t",colClasses="character")
    addIDassignSero6Ein <- cbind(checknotSero6CDid, checknotSero6CD, assignSero6E)

    for (i in 1:length(addIDassignSero6Ein$V2)){
      if (addIDassignSero6Ein$V2[i]=="T"&addIDassignSero6Ein$V3[i]=="Y")
        {addIDassignSero6Ein$V4[i]="6Bii(6E)"}
      else if (addIDassignSero6Ein$V2[i]=="S"&addIDassignSero6Ein$V3[i]=="Y")
        {addIDassignSero6Ein$V4[i]="6A/F"}
      else if (addIDassignSero6Ein$V2[i]=="A"&addIDassignSero6Ein$V3[i]=="Y")
        {addIDassignSero6Ein$V4[i]="6B"}
      else if (addIDassignSero6Ein$V2[i]=="S"&addIDassignSero6Ein$V3[i]=="D")
        {addIDassignSero6Ein$V4[i]="6G"}
      else if (addIDassignSero6Ein$V2[i]!="S"|addIDassignSero6Ein$V2[i]!="A"|addIDassignSero6Ein$V2[i]!="T")
        {addIDassignSero6Ein$V4[i]="6unk"}
    }

    idSero6E <- subset(addIDassignSero6Ein,V4=="6Bii(6E)")
    idSero6G <- subset(addIDassignSero6Ein,V4=="6G")
    idSero6AFB <- subset(addIDassignSero6Ein,V4=="6A/F"|V4=="6B"|V4=="6unk")
    assignSero6Efinal <- idSero6E[,c("checknotSero6CDid","V4")]
    assignSero6Gfinal <- idSero6G[,c("checknotSero6CDid","V4")]
    file6AFlist <- idSero6AFB[,c("checknotSero6CDid", "checknotSero6CD")]
          
##############    
## Identify 6F
##############
    if (nrow(file6AFlist) > 0) {

      listSero6AFfasta <- paste(file6AFlist$checknotSero6CD)
      listSero6AFfastaID <-paste(file6AFlist$checknotSero6CDid)
      listSero6AFTab <- paste(file6AFlist$checknotSero6CDid,".6AF.tab",sep="")

      for (i in 1:length(listSero6AFfasta)) {
        RunBlat(querySeqs = listSero6AFfasta[[i]], refNuc = assign6Fref, blatOutput = listSero6AFTab[[i]])
      }

      listSero6AFtab<-dir(pattern=".6AF.tab")
      listSero6AFtabFasta <- paste(listSero6AFtab,".fasta",sep="")

      for (i in 1:length(listSero6AFtab)) {
        runMview(mviewIn = listSero6AFtab[[i]], mviewOut = listSero6AFtabFasta[[i]])   
      }

      listSero6AFfastaIn <- dir(pattern=".6AF.tab.fasta")
      system2("mkdir","6AF")
    
      for (i in 1:length(listSero6AFfastaIn)) {
        runSplitMultifasta(multifasta = listSero6AFfastaIn[[i]], outDir = "6AF") 
      }

      system2("rm","6AF/6A_wciNalpha.fsa")
      fsaSero6AFlistT <- list.files("6AF",full.names=TRUE)
      fsaSero6AFlistF <- list.files("6AF",full.names=FALSE)
      pepSero6AF <- paste(fsaSero6AFlistF,".6AF.pep",sep="")

      for (i in 1:length(fsaSero6AFlistT)){
        runTranseq(TranseqIn = fsaSero6AFlistT[[i]], TranseqOut = pepSero6AF[[i]]) 
      }

      listpepSero6AF <- dir(pattern=".6AF.pep")
      readListPepSero6AF <- NULL

      for (i in 1:length(listpepSero6AF)) {
        readListPepSero6AF[[i]] <- read.fasta(listpepSero6AF[[i]],seqtype="AA",set.attributes=FALSE)
      }

      wciN150 <- NULL
      sero6AFout <- NULL

      for (i in 1:length(readListPepSero6AF)){
        wciN150[[i]] <- getFrag(readListPepSero6AF[[i]],begin=150,end=150)
      }

      for (i in 1:length(listpepSero6AF)){
        sero6AFout[[i]]<-cbind(listpepSero6AF[[i]],wciN150[[i]])
      }

      sero6AFoutWrite <- file("6AF_out",open="a")
      for (i in 1:length(sero6AFout)){
        write.table(sero6AFout[[i]],file=sero6AFoutWrite,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
      }
      close(sero6AFoutWrite)
 
      assignSero6F <- read.delim("6AF_out",header=F,sep="\t",colClasses="character")
      addIDassignSero6Fin <- cbind(listSero6AFfastaID, listSero6AFfasta, assignSero6F)
      for (i in 1:length(addIDassignSero6Fin$V2)) {
        if (addIDassignSero6Fin$V2[i]=="T")
          {addIDassignSero6Fin$V3[i]="6F"}
        else if (addIDassignSero6Fin$V2[i]!="T")
          {addIDassignSero6Fin$V3[i]="6A/6B"}
      }

      assignSero6Fout <- subset(addIDassignSero6Fin,V3=="6F")    
      assignSero6Ffinal <- assignSero6Fout[,c("listSero6AFfastaID","V3")]
      assignSero6ABout <- subset(addIDassignSero6Fin,V3=="6A/6B")
    } else if (nrow(file6AFlist)==0) {
      cat("No serotype 6F found","\n")
    }

#####################
## Identify 6A and 6B    
#####################

    if (exists("assignSero6ABout") == TRUE) {
    
      checkSero6AB <- paste(assignSero6ABout$listSero6AFfasta)
      checkSero6ABid <- assignSero6ABout$listSero6AFfastaID
      checkSero6ABtab <- paste(checkSero6ABid,".6AB.tab",sep="")

      for (i in 1:length(checkSero6AB)) {
        RunBlat(querySeqs = checkSero6AB[[i]], refNuc = assign6ABref, blatOutput = checkSero6ABtab[[i]])
      }

      listSero6ABtab<-dir(pattern=".6AB.tab")
      listSero6ABtabFasta <- paste(listSero6ABtab,".fasta",sep="")

      for (i in 1:length(listSero6ABtab)) {
        runMview(mviewIn = listSero6ABtab[[i]], mviewOut = listSero6ABtabFasta[[i]])
      }

      listSero6ABfasta <- dir(pattern=".6AB.tab.fasta")
      system2("mkdir","6AB")

      for (i in 1:length(listSero6ABfasta)) {
        runSplitMultifasta(multifasta = listSero6ABfasta[[i]], outDir = "6AB") 
      }

      system2("rm","6AB/wciP.fsa")
      fsaSero6ABlistT <- list.files("6AB",full.names=TRUE)
      fsaSero6ABlistF <- list.files("6AB",full.names=FALSE)
      pepSero6AB <- paste(fsaSero6ABlistF,".6AB.pep",sep="")
      
      for (i in 1:length(fsaSero6ABlistT)){
        runTranseq(TranseqIn = fsaSero6ABlistT[[i]], TranseqOut = pepSero6AB[[i]])
      }

      listSero6ABpep <- dir(pattern="6AB.pep")
      readlistSero6ABpep <- NULL
      for (i in 1:length(listSero6ABpep)) {
        readlistSero6ABpep[[i]] <-read.fasta(listSero6ABpep[[i]],seqtype="AA",set.attributes=FALSE)
      }

      wciP195 <- NULL
      outSero6AB<-NULL

      for (i in 1:length(readlistSero6ABpep)){
        wciP195[[i]] <- getFrag(readlistSero6ABpep[[i]],begin = 195, end=195)
      }

      for (i in 1:length(readlistSero6ABpep)){
        outSero6AB[[i]]<-cbind(listSero6ABpep[[i]],wciP195[[i]])
      }

      outSero6ABwrite <- file("6AB_out",open="a")
      for (i in 1:length(outSero6AB)){
        write.table(outSero6AB[[i]], file = outSero6ABwrite, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
      }
      close(outSero6ABwrite)

      assignSero6AB <- read.delim("6AB_out",header=F,sep="\t",colClasses="character")
      addIDassignSero6Ein <- cbind(checkSero6ABid, checkSero6AB, assignSero6AB)

      for (i in 1:length(addIDassignSero6Ein$V2)) {
        if (addIDassignSero6Ein$V2[i]=="S")
          {addIDassignSero6Ein$V3[i]="6A"}
        else if (addIDassignSero6Ein$V2[i]=="N")
          {addIDassignSero6Ein$V3[i]="6B"}
        else if (addIDassignSero6Ein$V2[i]!="N"|addIDassignSero6Ein$V2[i]!="S")
          {addIDassignSero6Ein$V3[i]="6unk"}
      }
      assignSero6ABfinal <- addIDassignSero6Ein[,c("checkSero6ABid","V3")]
  
    } else if (exists("assignSero6ABout") != TRUE) {
        cat("No serotypes 6A or 6B found","\n")
    }
  } else if (nrow(assignOtherSero6)==0) {
      cat("No serotypes 6A, 6B, 6Bii(6E), 6F or 6G found","\n")
  }
} else if (nrow(sero6Results)==0) {
  cat("No serotype 6 sequences found","\n")
}

###############################
## Identify serotypes 7A and 7F
###############################

if (nrow(sero7Results)>0) {

  fasta7 <- sero7Results$filePath
  fasta7ID <- sero7Results$id
  tab7 <- paste(fasta7ID,".7all.tab",sep="")
  n7 <- length(fasta7)

  for (i in 1:n7) {
    RunBlatBlast8(querySeqs = fasta7[[i]], refNuc = assign7AFref, blatOutput = tab7[[i]])
  }

  tabinSero7 <- dir(pattern=".7all.tab")

  outSero7<-NULL
  for (i in 1:length(tabinSero7)) outSero7[[i]] <- read.delim(tabinSero7[i],header=F,sep="\t")

  sortoutSero7<-NULL
  for (i in 1:length(tabinSero7)) sortoutSero7[[i]] <- outSero7[[i]][order(-outSero7[[i]]$V12),]

  headSero7<-NULL
  for (i in 1:length(tabinSero7)) headSero7[[i]] <- sortoutSero7[[i]][1,]

  headSero7out <- file("7_head_out",open="a")
  for (i in 1:length(tabinSero7)) {
    write.table(headSero7[[i]],file=headSero7out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  close(headSero7out)

  blatresSero7 <- read.delim("7_head_out",header=F,sep="\t")
  
  for (i in 1:length(blatresSero7$V7)) {
    if (blatresSero7$V4[i]>900)
      {blatresSero7$V13[i]="7F"}
    else if (blatresSero7$V4[i]<700)
      {blatresSero7$V13[i]="7A"}
  }

  addSero7data <- cbind(fasta7ID,blatresSero7)
  assignSero7final <- addSero7data[,c("fasta7ID","V13")]
  
} else if (nrow(sero7Results)==0) {
  cat("No serotype 7 sequences found","\n")
}

###############################
## Identify serotypes 9A and 9V
###############################

if (nrow(sero9Results)>0) {

  fasta9 <- sero9Results$filePath
  fasta9ID <- sero9Results$id
  tab9 <- paste(fasta9ID,".9all.tab", sep="")
  n9 <- length(fasta9)
  
  for (i in 1:n9) {
    RunBlatBlast8(querySeqs = fasta9[[i]], refNuc = assign9AVref, blatOutput = tab9[[i]])
  }

  tabinSero9 <- dir(pattern=".9all.tab")

  sero9out <- NULL
  for (i in 1:length(tabinSero9)) sero9out[[i]] <- read.delim(tabinSero9[i],header=F,sep="\t")

  sortSero9out <- NULL
  for (i in 1:length(tabinSero9)) sortSero9out[[i]] <- sero9out[[i]][order(-sero9out[[i]]$V12),]

  headSero9 <- NULL
  for (i in 1:length(tabinSero9)) headSero9[[i]] <- sortSero9out[[i]][1,]

  headSero9out <- file("9_head_out",open="a")
  for (i in 1:length(tabinSero9)) {
    write.table(headSero9[[i]],file=headSero9out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  close(headSero9out)

  blatresSero9 <- read.delim("9_head_out",header=F,sep="\t")

  for (i in 1:length(blatresSero9$V7)) {
    if (blatresSero9$V4[i]>1000)
      {blatresSero9$V13[i]="9V"}
    else if (blatresSero9$V4[i]<900)
      {blatresSero9$V13[i]="9A"}
  }

  addSero9data <- cbind(fasta9ID,blatresSero9)
  assignSero9final <- addSero9data[,c("fasta9ID","V13")]

} else if (nrow(sero9Results)==0) {
  cat("No serotype 9 sequences found","\n")
}

######################################
## Identify serotypes 11A, 11D and 11E
######################################

if (nrow(sero11Results)>0) {

  fasta11 <- sero11Results$filePath
  fasta11ID <- sero11Results$id
  id11 <- sero11Results[,c("id","filePath")]
  tab11 <- paste(fasta11ID,".11AD.tab",sep="")
  n11 <- length(fasta11)

#################################
## Identify serotypes 11A and 11D
#################################

  for (i in 1:n11) {
    RunBlat(querySeqs = fasta11[[i]], refNuc = assign11ADref, blatOutput = tab11[[i]])
  }

  listSero11ADtab <- dir(pattern=".11AD.tab")
  listSero11ADtabFasta <- paste(listSero11ADtab,".fasta",sep="")

  for (i in 1:length(listSero11ADtab)) {
    runMview(mviewIn = listSero11ADtab[[i]], mviewOut = listSero11ADtabFasta[[i]]) 
  }

  listSero11ADfasta <- dir(pattern=".11AD.tab.fasta")
  system2("mkdir","11AD")

  for (i in 1:length(listSero11ADfasta)) {
    runSplitMultifasta(multifasta = listSero11ADfasta[[i]], outDir = "11AD")
  }

  system2("rm","11AD/11A_wcrL.fsa")
  fsaSero11ADlistT <- list.files("11AD",full.names=TRUE)
  fsaSero11ADlistF <- list.files("11AD",full.names=FALSE)
  pepSero11AD <- paste(fsaSero11ADlistF,".11AD.pep",sep="")

  for (i in 1:length(fsaSero11ADlistT)){
    runTranseq(TranseqIn = fsaSero11ADlistT[[i]], TranseqOut = pepSero11AD[[i]]) 
  }

  listSero11ADpep <- dir(pattern=".11AD.pep")
  readListSero11ADpep <- NULL

  for (i in 1:length(listSero11ADpep)) {
    readListSero11ADpep[[i]] <- read.fasta(listSero11ADpep[[i]],seqtype="AA",set.attributes=FALSE)
  }

  wcrL112 <- NULL
  outSero11AD <- NULL

  for (i in 1:length(readListSero11ADpep)){
    wcrL112[[i]] <- getFrag(readListSero11ADpep[[i]],begin=112,end=112)
  }

  for (i in 1:length(listSero11ADpep)){
    outSero11AD[[i]]<-cbind(listSero11ADpep[[i]],wcrL112[[i]])
  }

  outSero11ADwrite <- file("11AD_out",open="a")
  for (i in 1:length(outSero11AD)){
    write.table(outSero11AD[[i]],file=outSero11ADwrite,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  close(outSero11ADwrite)

  assignSero11D <- read.delim("11AD_out",header=F,sep="\t",colClasses="character")
  addIDassignSero11Din <- cbind(id11, assignSero11D)

  for (i in 1:length(addIDassignSero11Din$V2)) {
    if (addIDassignSero11Din$V2[i]=="N")
      {addIDassignSero11Din$V3[i]="11A"}
    else if (addIDassignSero11Din$V2[i]=="S")
      {addIDassignSero11Din$V3[i]="11D"}
  }

  idSero11D <- subset(addIDassignSero11Din,V3=="11D")
  idSero11AE <- subset(addIDassignSero11Din,V3=="11A")
  assignSero11Dfinal <- idSero11D[,c("id","V3")]

###############
## Identify 11E
###############

  if (nrow(idSero11AE)>0) {  

    sero11AEfasta <- idSero11AE$filePath
    sero11AEid <- idSero11AE$id
    sero11AEtab <- paste(sero11AEid,".11AE.tab",sep="")

    for (i in 1:length(sero11AEfasta)) {
      RunBlatBlast8(querySeqs = sero11AEfasta[[i]], refNuc = assign11Eref, blatOutput = sero11AEtab[[i]])
    }

    tabinSero11AE <- dir(pattern=".11AE.tab")
    Sero11AEout<-NULL
    for (i in 1:length(tabinSero11AE)) Sero11AEout[[i]] <- read.delim(tabinSero11AE[i],header=F,sep="\t")

    sortSero11AEout<-NULL
    for (i in 1:length(tabinSero11AE)) sortSero11AEout[[i]] <- Sero11AEout[[i]][order(-Sero11AEout[[i]]$V12),]

    headSero11AE<-NULL
    for (i in 1:length(tabinSero11AE)) headSero11AE[[i]] <- sortSero11AEout[[i]][1,]

    headSero11AEout <- file("11AE_head_out",open="a")

    for (i in 1:length(tabinSero11AE)) {
      write.table(headSero11AE[[i]],file=headSero11AEout,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
    close(headSero11AEout)

    blatresSero11AE <- read.delim("11AE_head_out",header=F,sep="\t")
    addSero11AEdata <- cbind(sero11AEid, blatresSero11AE)

    for (i in 1:length(addSero11AEdata$V7)) {
      if (addSero11AEdata$V4[i]>1000)
        {addSero11AEdata$V13[i]="11A"}
      else if (addSero11AEdata$V4[i]<900)
        {addSero11AEdata$V13[i]="11E"}
    }

    assignSero11AEfinal <- addSero11AEdata[,c("sero11AEid","V13")]

  } else if (nrow(idSero11AE)==0) {
    cat("No serotype 11A/E found","\n")
  }
} else if (nrow(sero11Results)==0) {
  cat("No serotype 11 sequences found","\n")
}

#################################
## Identify serotypes 18B and 18C
#################################

if (nrow(sero18Results)>0) {

  fasta18 <- sero18Results$filePath
  fasta18ID <- sero18Results$id
  tab18 <- paste(fasta18ID,".18all.tab",sep="")
  n18 <- length(fasta18)

  for (i in 1:n18) {
    RunBlatBlast8(querySeqs = fasta18[[i]], refNuc = assign18BCref, blatOutput = tab18[[i]])
  }
  
  tabinSero18 <- dir(pattern=".18all.tab")
  sero18out<-NULL
  for (i in 1:length(tabinSero18)) sero18out[[i]] <- read.delim(tabinSero18[i],header=F,sep="\t")

  sortSero18out<-NULL
  for (i in 1:length(tabinSero18)) sortSero18out[[i]] <- sero18out[[i]][order(-sero18out[[i]]$V12),]

  headSero18<-NULL
  for (i in 1:length(tabinSero18)) headSero18[[i]] <- sortSero18out[[i]][1,]

  headSero18out <- file("18_head_out",open="a")
  for (i in 1:length(tabinSero18)) {
    write.table(headSero18[[i]],file=headSero18out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  close(headSero18out)

  blatresSero18 <- read.delim("18_head_out",header=F,sep="\t")

  for (i in 1:length(blatresSero18$V7)) {
    if (blatresSero18$V4[i]>900)
      {blatresSero18$V13[i]="18C"}
    else if (blatresSero18$V4[i]<900)
      {blatresSero18$V13[i]="18B"}
  }

  addSero18data <- cbind(fasta18ID,blatresSero18)
  assignSero18final <- addSero18data[,c("fasta18ID","V13")]

} else if (nrow(sero18Results)==0) {
  cat("No serotype 18 sequences found","\n")
}

#################################
## Identify serotypes 20A and 20B
#################################

if (nrow(sero20Results)>0) {

  fasta20 <- sero20Results$filePath
  fasta20ID <- sero20Results$id
  tab20 <- paste(fasta20ID,".20all.tab",sep="")
  n20 <- length(fasta20)

  for (i in 1:n20) {
    RunBlatBlast8(querySeqs = fasta20[[i]], refNuc = assign20ABref, blatOutput = tab20[[i]])
  }

  tabinSero20 <- dir(pattern=".20all.tab")

  sero20out<-NULL
  for (i in 1:length(tabinSero20)) sero20out[[i]] <- read.delim(tabinSero20[i],header=F,sep="\t")

  sortSero20out<-NULL
  for (i in 1:length(tabinSero20)) sortSero20out[[i]] <- sero20out[[i]][order(-sero20out[[i]]$V12),]

  headSero20<-NULL
  for (i in 1:length(tabinSero20)) headSero20[[i]] <- sortSero20out[[i]][1,]

  headSero20out <- file("20_head_out",open="a")
  for (i in 1:length(tabinSero20)) {
    write.table(headSero20[[i]],file=headSero20out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  close(headSero20out)

  blatresSero20 <- read.delim("20_head_out",header=F,sep="\t")
  for (i in 1:length(blatresSero20$V7)) {
    if (blatresSero20$V4[i]<900)
      {blatresSero20$V13[i]="20A"}
    else if (blatresSero20$V4[i]>900)
      {blatresSero20$V13[i]="20B"}
  }

  addSero20data <- cbind(fasta20ID, blatresSero20)
  assignSero20final <- addSero20data[,c("fasta20ID", "V13")]

} else if (nrow(sero20Results)==0) {
  cat("No serotype 20 sequences found","\n")
}

#################################
## Identify serotypes 22A and 22F
#################################

if (nrow(sero22Results)>0) {

  fasta22 <- sero22Results$filePath
  fasta22ID <- sero22Results$id
  tab22 <- paste(fasta22ID,".22all.tab",sep="")
  n22 <- length(fasta22)

  for (i in 1:n22) {
    RunBlastall(BARef = assign22AFref, BAquerySeqs = fasta22[[i]], blastAllOut = tab22[[i]])
  }

  file22list<-dir(pattern=".22all.tab")
  
  out22Read <- NULL
  for (i in 1:length(file22list)) out22Read[[i]] <- read.delim(file22list[i],header=F, sep="\t")

  sort22ReadOut <- NULL
  for (i in 1:length(out22Read)) sort22ReadOut[[i]] <- out22Read[[i]][order(-out22Read[[i]]$V12),]

  head22wcwA <- NULL
  for (i in 1:length(sort22ReadOut)) head22wcwA[[i]] <- sort22ReadOut[[i]][1,]

  head22wcwAout <- file("22_head_out", open="a")
  for (i in 1:length(head22wcwA)) {
    write.table(head22wcwA[[i]], file=head22wcwAout, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
  close(head22wcwAout)
  
  assignSero22AFin <- read.delim("22_head_out",header=F,sep="\t")
  addIDassignSero22AFin <- cbind(fasta22ID, assignSero22AFin)

  for (i in 1:length(addIDassignSero22AFin$fasta22ID)) {
    if (addIDassignSero22AFin$V4[i]<500)
      {addIDassignSero22AFin$V13[i]="22A"}
    else if (addIDassignSero22AFin$V4[i]>500)
      {addIDassignSero22AFin$V13[i]="22F"}
  }

  assignSero22AFfinal <- addIDassignSero22AFin[,c("fasta22ID","V13")]
  
} else if (nrow(sero22Results)==0) {
  cat("No serotype 22 sequences found","\n")
}

#####################################
## Identify serotypes 33A, 33F and 37
#####################################

if (nrow(sero33Results)>0) {

  fasta33 <- sero33Results$filePath
  fasta33ID <- sero33Results$id
  id33 <- sero33Results[,c("id","filePath")]

#######################
## Identify serotype 37  
#######################

  out37 <- paste(fasta33ID,".37all.out",sep="")
  n33 <- length(fasta33)

  for (i in 1:n33) {
    RunBlastall(BARef = assign37ref, BAquerySeqs = fasta33[[i]], blastAllOut = out37[[i]])
  }

  sero37outIn <- dir(pattern=".37all.out")

  sero37outRead <- NULL
  for (i in 1:length(sero37outIn)) sero37outRead[[i]] <- read.delim(sero37outIn[i],header=F, sep="\t")

  sortsero37outRead <- NULL
  for (i in 1:length(sero37outRead)) sortsero37outRead[[i]] <- sero37outRead[[i]][order(-sero37outRead[[i]]$V12),]

  headSero37tts <- NULL
  for (i in 1:length(sortsero37outRead)) headSero37tts[[i]] <- sortsero37outRead[[i]][1,]

  headSero37ttsOut <- file("37_head_out", open="a")
  for (i in 1:length(headSero37tts)) {
    write.table(headSero37tts[[i]],file=headSero37ttsOut,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  close(headSero37ttsOut)

  resSero37 <- read.delim("37_head_out",header=F,sep="\t")
  addIDresSero37 <- cbind(id33, resSero37)
  
  assignSero33AF <- subset(addIDresSero37, V4<500)
  assignSero37 <- subset(addIDresSero37, V4>500)
  
  if (nrow(assignSero37) > 0) {
    assignSero37$V13 <- "37"
    assignSero37out <- assignSero37[,c("id","V13")]
  } else if (nrow(assignSero37)==0) {
    cat("No serotype 37 found","\n")
  }

#################################
## Identify serotypes 33A and 33F  
#################################

  fasta33AF <- assignSero33AF$filePath
  fasta33AFid <- assignSero33AF$id
  tab33 <- paste(fasta33AFid, ".33all.tab", sep="")
  n33 <- length(fasta33AF)

  for (i in 1:n33) {
    RunBlatBlast8(querySeqs = fasta33AF[[i]], refNuc = assign33AFref, blatOutput = tab33[[i]])
  }

  tabinSero33 <- dir(pattern=".33all.tab")

  sero33out<-NULL
  for (i in 1:length(tabinSero33)) sero33out[[i]] <- read.delim(tabinSero33[i],header=F,sep="\t")

  sortSero33out<-NULL
  for (i in 1:length(tabinSero33)) sortSero33out[[i]] <- sero33out[[i]][order(-sero33out[[i]]$V12),]

  headSero33<-NULL
  for (i in 1:length(tabinSero33)) headSero33[[i]] <- sortSero33out[[i]][1,]

  headSero33out <- file("33_head_out",open="a")
  for (i in 1:length(tabinSero33)) {
    write.table(headSero33[[i]],file=headSero33out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  close(headSero33out)

  blatresSero33 <- read.delim("33_head_out",header=F,sep="\t")
  addIDresSero33 <- cbind(fasta33AFid, blatresSero33)

  for (i in 1:length(addIDresSero33$V7)) {
    if (addIDresSero33$V4[i]<400)
      {addIDresSero33$V13[i]="33F"}
    else if (addIDresSero33$V4[i]>1000)
      {addIDresSero33$V13[i]="33A"}
  }

  assignSero33AFfinal <- addIDresSero33[,c("fasta33AFid","V13")]

} else if (nrow(sero33Results)==0) {
  cat("No serotype 33 sequences found","\n")
}

#################################
## Identify serotypes 25A and 25F
#################################

if (nrow(sero25Results)>0) {

  fasta25 <- sero25Results$filePath
  fasta25ID <- sero25Results$id
  tab25 <- paste(fasta25ID, ".25all.tab", sep="")
  n25 <- length(fasta25)

  for (i in 1:n25) {
    RunBlatBlast8(querySeqs = fasta25[[i]], refNuc = assign25AFref, blatOutput = tab25[[i]])
  }
  
  tabinSero25 <- dir(pattern=".25all.tab")

  sero25out <- NULL
  for (i in 1:length(tabinSero25)) sero25out[[i]] <- read.delim(tabinSero25[i],header=F,sep="\t")

  sortSero25out <- NULL
  for (i in 1:length(tabinSero25)) sortSero25out[[i]] <- sero25out[[i]][order(-sero25out[[i]]$V12),]

  headSero25 <- NULL
  for (i in 1:length(tabinSero25)) headSero25[[i]] <- sortSero25out[[i]][1,]

  headSero25out <- file("25_head_out",open="a")
  for (i in 1:length(tabinSero25)) {
    write.table(headSero25[[i]],file=headSero25out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }
  close(headSero25out)

  blatresSero25 <- read.delim("25_head_out",header=F,sep="\t")

  for (i in 1:length(blatresSero25$V5)) {
    if (blatresSero25$V5[i]==0)
      {blatresSero25$V13[i]="25F"}
    else if (blatresSero25$V4[i]!=0)
      {blatresSero25$V13[i]="25A"}
  }

  addSero25data <- cbind(fasta25ID, blatresSero25) 
  assignSero25AFfinal <- addSero25data[,c("fasta25ID","V13")]

} else if (nrow(sero25Results)==0) {
  cat("No serotype 25 sequences found","\n")
}

####################################################################################
## Final output: combine all serotype prediction files into Serotype_predictions.txt
####################################################################################

final_out <- file("Serotype_predictions_unsorted.txt",open="a")
if (exists("NTout")==TRUE) {
  write.table(NTout,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("NTout")!=TRUE) {
  cat("","\n")
}
if (exists("resttoBlatOut")==TRUE) {
  write.table(resttoBlatOut,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("resttoBlatOut")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero6CDfinal")==TRUE) {
  write.table(assignSero6CDfinal,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero6CDfinal")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero6Efinal")==TRUE) {
  write.table(assignSero6Efinal,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero6Efinal")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero6Ffinal")==TRUE) {
  write.table(assignSero6Ffinal,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero6Ffinal")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero6Gfinal")==TRUE) {
  write.table(assignSero6Gfinal,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero6Gfinal")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero6ABfinal")==TRUE) {
  write.table(assignSero6ABfinal,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero6ABfinal")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero22AFfinal")==TRUE) {
  write.table(assignSero22AFfinal,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero22AFfinal")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero37out")==TRUE) {
  write.table(assignSero37out,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero37out")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero33AFfinal")==TRUE) {
  write.table(assignSero33AFfinal,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero33AFfinal")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero20final")==TRUE) {
  write.table(assignSero20final,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero20final")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero11Dfinal")==TRUE) {
  write.table(assignSero11Dfinal,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero11Dfinal")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero11AEfinal")==TRUE) {
  write.table(assignSero11AEfinal,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero11AEfinal")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero9final")==TRUE) {
  write.table(assignSero9final,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero9final")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero7final")==TRUE) {
  write.table(assignSero7final,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero7final")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero18final")==TRUE) {
  write.table(assignSero18final,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero18final")!=TRUE) {
  cat("","\n")
}
if (exists("assignSero25AFfinal")==TRUE) {
  write.table(assignSero25AFfinal,file= final_out,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
} else if (exists("assignSero25AFfinal")!=TRUE) {
  cat("","\n")
}
close(final_out)

#########################################
## Sort final serotype file by isolate id
#########################################

finalSerotypes <- read.delim("Serotype_predictions_unsorted.txt", header=F, sep="\t")
sortFinalSerotypes <- finalSerotypes[order(finalSerotypes$V1),]
names(sortFinalSerotypes) <- c("id", "Serotype")
write.table(sortFinalSerotypes, file = "Serotype_predictions.txt", sep="\t", row.names=FALSE,col.names=TRUE,quote=FALSE)

########################################################################################
## Copy Serotype_predictions.txt to main seqSerotyper directory and delete tmp directory
########################################################################################

seroPredCopy <- file.copy(from = "Serotype_predictions.txt", to = seqSerotyperDir)

setwd(seqSerotyperDir)

tmpDelete <- unlink("tmp", recursive = TRUE)

cat("Serotype prediction complete: see Serotype_predictions.txt","\n")
