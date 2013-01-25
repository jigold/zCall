#!/usr/bin/env Rscript

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# April 5th, 2012


## Input is means and standard deviations of the homozygote clusters from common sites created from findMeanSD.py


args <- commandArgs(TRUE)

X <- read.table(args[1],header=TRUE,sep="\t") #read in file as a dataframe
outFile <- args[2] # Specify output file
weighted <- args[3] # If 1, use weighted linear regression

w <- ((1/ X$nMinorHom) + (1 / X$nCommonHom))^-0.5
print(summary(w))
#-------------------------------------
# LM of meanY as a function of meanX
if ( weighted == 1 ){
m <- lm(X$meanY ~ X$meanX, weights = w)}
if ( weighted == 0 ){
m <- lm(X$meanY ~ X$meanX)}

a <- as.numeric(summary(m)$coefficients[1]) # beta0
b <- as.numeric(summary(m)$coefficients[2]) # beta1
pa <- as.numeric(summary(m)$coefficients[7]) # pvalue for beta0
pb <- as.numeric(summary(m)$coefficients[8]) # pvalue for beta1

### Correlation -- uncomment
#print(summary(m, correlation = TRUE))

#-------------------------------------
# LM of meanX as a function of meanY
if ( weighted == 1 ){
x <- lm(X$meanX ~ X$meanY, weights = w)}
if ( weighted == 0 ){
x <- lm(X$meanX ~ X$meanY)}

c <- as.numeric(summary(x)$coefficients[1]) # beta0
d <- as.numeric(summary(x)$coefficients[2]) # beta1
pc <- as.numeric(summary(x)$coefficients[7]) # pvalue for beta0
pd <- as.numeric(summary(x)$coefficients[8]) # pvalue for beta1

### Correlation -- uncomment
#print(summary(x, correlation = TRUE))

#-------------------------------------
# LM of sdY as a function of sdX
if ( weighted == 1 ){
l <- lm(X$sdY ~ X$sdX, weights = X$nMinorHom)}
if ( weighted == 0 ){
l <- lm(X$sdY ~ X$sdX)}

e <- as.numeric(summary(l)$coefficients[1]) # beta0
f <- as.numeric(summary(l)$coefficients[2]) # beta1
pe <- as.numeric(summary(l)$coefficients[7]) # pvalue for beta0
pf <- as.numeric(summary(l)$coefficients[8]) # pvalue for beta1

### Correlation -- uncomment
#print(summary(l, correlation = TRUE))

#-------------------------------------
# LM of sdX as a function of sdY
if ( weighted == 1 ){
y <- lm(X$sdX ~ X$sdY, weights = w)}
if ( weighted == 0 ){
y <- lm(X$sdX ~ X$sdY)}

g <- as.numeric(summary(y)$coefficients[1]) # beta0
h <- as.numeric(summary(y)$coefficients[2]) # beta1
pg <- as.numeric(summary(y)$coefficients[7]) # pvalue for beta0
ph <- as.numeric(summary(y)$coefficients[8]) # pvalue for beta1

### Correlation -- uncomment
#print(summary(y, correlation = TRUE))

#-------------------------------------
## Output to betas.txt

labels <- c("meanY~meanX","meanX~meanY","sdY~sdX","sdX~sdY")

beta0 <- c(a,c,e,g)
beta1 <- c(b,d,f,h)
p0 <- c(pa,pc,pe,pg)
p1 <- c(pb,pd,pf,ph)

head <- c("Model","Beta0","Beta1","pBeta0","pBeta1")

t <- cbind(labels, beta0, beta1, p0, p1)
r <- rbind(head, t)

write.table(as.matrix(r), file = paste(outFile,sep=""),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

