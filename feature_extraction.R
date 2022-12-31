#Note: commented peptider library functions beacuse they are too slow
#Pre-requisites: Column name should be same as the provided input file

library(tools)

#exception handling
if (length(commandArgs(TRUE)) != 1) {
  stop("There must be 1 argument, name of the csv file name after R script name.", call.=FALSE)
} else {
  file <- commandArgs(TRUE)
}

if (file_ext(file) != "csv") {
  stop("File extension must be csv", call.=FALSE)
}

if (file.exists(file)) {
  datas = read.csv(file)
} else {
  stop("The file does not exist", call.=FALSE)
}

library(Peptides)
library(peptider)

# Loading a property to evaluate its autocorrelation
data(AAdata)
  
#input file
data = read.csv(file)

#renaming columns of input file
names(data) <- c("Peptide.Sequence", "Target")

#checking for irregularities in input file

#log file
log_txt <- file("log.txt")
logs <- c("Logs")

ind <- 0
for (i in data$Peptide.Sequence) {
  ind <- ind + 1
  if (i == "") {
    logs <- append(logs, paste("Empty cell in column 1 row", ind))
  }
}

for (i in data$Target) {
  if (is.na(i)) {
    logs <- append(logs, paste("Empty cell in column 2 row", ind))
  }
}
writeLines(logs, log_txt)
close(log_txt)

#removing irregularities from data (if any)
data <- data[(data$Peptide.Sequence != "" & !is.na(data$Target)),]

#list of empty vectors
number_of_features <- 56
features = list()
for (i in 1:number_of_features) {
  features[i] <- c()
}

#calculations
k <- 0
for (seq in data$Peptide.Sequence){
  #list index
  i <- 0
  
  #sequence name
  features[[i <- i + 1]] <- append(features[[i + 1]], seq)
  
  #aliphatic index
  features[[i <- i + 1]] <- append(features[[i + 1]], aIndex(seq))
  
  #boman index
  features[[i <- i + 1]] <- append(features[[i + 1]], boman(seq))
  
  #instaIndex
  features[[i <- i + 1]] <- append(features[[i + 1]], instaIndex(seq))
  
  #Probability of detection of peptides (takes too much time to compile)
  #features[[i <- i + 1]] <- append(features[[i + 1]], ppeptide(seq, libscheme = "NNK", N=10^8))
  
  #auto-correlation index for lag=1
  features[[i <- i + 1]] <- append(features[[i + 1]], autoCorrelation(sequence = seq, lag = 1, property = AAdata$Hydrophobicity$KyteDoolittle, center = TRUE))
  
  #auto-covariance index for lag=1
  features[[i <- i + 1]] <- append(features[[i + 1]], autoCovariance(sequence = seq, lag = 1, property = AAdata$Hydrophobicity$KyteDoolittle, center = TRUE))
  
  #computed average of BLOSUM indices of all the amino acids
  for (j in 1:10) {
    features[[i <- i + 1]] <- append(features[[i + 1]], blosumIndices(seq)[[1]][[j]])
  }
  
  #net charge of a protein sequence 
  features[[i <- i + 1]] <- append(features[[i + 1]], charge(seq, pH = 7, pKscale = "Lehninger"))
  
  #cross-covariance index for a lag=1
  features[[i <- i + 1]] <- append(features[[i + 1]], crossCovariance(seq, lag = 1, property1 = AAdata$Hydrophobicity$KyteDoolittle, property2 = AAdata$Hydrophobicity$Eisenberg, center = TRUE))
  
  #amino acid composition: Number
  for (j in 1:9) {
    features[[i <- i + 1]] <- append(features[[i + 1]], aaComp(seq)[[1]][j, 1])
  }
  #amino acid composition: % of moles
  for (j in 1:9) {
    features[[i <- i + 1]] <- append(features[[i + 1]], aaComp(seq)[[1]][j, 2])
  }
  
  #Cruciani properties of an amino-acids sequence (PP1, PP2, PP3)
  for (j in 1:3) {
    features[[i <- i + 1]] <- append(features[[i + 1]], crucianiProperties(seq)[[1]][[j]])
  }
  
  #FASGAI vectors (Hydrophobicity index, Alpha and turn propensities, Bulky properties, Compositional characteristic index, Local flexibility, Electronic properties)
  for (j in 1:6) {
    features[[i <- i + 1]] <- append(features[[i + 1]], fasgaiVectors(seq)[[1]][[j]])
  }
  
  #hydrophobic moment
  features[[i <- i + 1]] <- append(features[[i + 1]], hmoment(seq, angle = 100, window = 11))
  
  #amino acid length
  features[[i <- i + 1]] <- append(features[[i + 1]], lengthpep(seq))
  
  
  #molecular weight of a protein sequence
  features[[i <- i + 1]] <- append(features[[i + 1]], mw(seq, monoisotopic = FALSE))
  
  #number of codon representations (takes too much time to compile)
  #features[[i <- i + 1]] <- append(features[[i + 1]], codons(seq, libscheme="NNK"))
  
  #Number of possible neighbours (takes too much time to compile)
  #features[[i <- i + 1]] <- append(features[[i + 1]], getNofNeighbors(seq))
  
  #isoelectic point 
  features[[i <- i + 1]] <- append(features[[i + 1]], pI(seq, pKscale= "Murray"))
  
  #T-scales
  for (j in 1:5) {
    features[[i <- i + 1]] <- append(features[[i + 1]], tScales(seq)[[1]][[j]])
  }
  
  #target
  features[[i <- i + 1]] <- append(features[[i + 1]], data$Target[k <- k + 1])
}

#column naming of list
names(features) <- c('Sequence', 'aIndex', 'boman', 'instaIndex', 'autoCorrelation', 'autoCovariance', 'blosum1', 'blosum2', 'blosum3', 
                     'blosum4', 'blosum5', 'blosum6', 'blosum7', 'blosum8', 'blosum9', 'blosum10', 'charge', 'crossCovariance', 'Tiny', 
                     'Small', 'Aliphatic', 'Aromatic', 'NonPolar', 'Polar', 'Charged', 'Basic', 'Acidic', 'Tiny mole%', 'Small mole%', 
                     'Aliphatic mole%', 'Aromatic mole%', 'NonPolar mole%', 'Polar mole%', 'Charged mole%', 'Basic mole%', 'Acidic mole%', 
                     'PP1', 'PP2', 'PP3', 'Hydrophobicity index', 'Alpha and turn propensities', 'Bulky properties', 'Compositional characteristic index', 
                     'Local flexibility', 'Electronic properties', 'hydrophobic moment', 'AA length', 'molecular weight', 'isoelectic point', 'T1', 'T2', 'T3', 
                     'T4', 'T5', 'target')

#list of vectors to dataframe
df <- as.data.frame(do.call(cbind, features))

#writing dataframe in csv file
write.csv(df, "output.csv", row.names = FALSE)

cat("Done")





#checking if given sequence has letters that belongs to Amino Acids
#would take a lot of runtime
#aaComp("AC")
#if ("a" in aaList()[[1]]) {
#  print("yes")
#}
#"a" == aaList()[[1]]
#toupper('ac') %in% aaList()