############################################################
#                                                          #
# Calculation of GWAS sig. threshold based on heritability #
#                                                          #
############################################################

# This code is adapted from https://doi.org/10.1186/s12864-019-5992-7 and altered by the author of this project

#rICE
D<- read.big.matrix("GNrice.txt", type="char", sep="\t", head = TRUE)
dim(D)
D=D[,2:34849]
D1=as.data.frame(as.matrix(D))
t(D1)
D2=t(D1)


# Simulation with 10 QTL
QTL <- 100*(1:10) #pick 10 QTL
u <- rep(0,34848) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(D2,u))
h2 <- 0.1
y <- g + rnorm(135,mean=0,sd=sqrt((1-h2)/h2*var(g)))
write.table(y, "H10Q10.txt", sep="\t")


# Simulation with 20 QTL
QTL <- 100*(1:20) #pick 20 QTL
u <- rep(0,34848) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(D2,u))
h2 <- 0.1
y <- g + rnorm(135,mean=0,sd=sqrt((1-h2)/h2*var(g)))
write.table(y, "H10Q20.txt", sep="\t")


# Simulation with 30 QTL
QTL <- 100*(1:30) #pick 30 QTL
u <- rep(0,34848) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(D2,u))
h2 <- 0.1
y <- g + rnorm(135,mean=0,sd=sqrt((1-h2)/h2*var(g)))
write.table(y, "H10Q30.txt", sep="\t")


# Simulation with 40 QTL
QTL <- 100*(1:40) #pick 40 QTL
u <- rep(0,34848) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(D2,u))
h2 <- 0.1
y <- g + rnorm(135,mean=0,sd=sqrt((1-h2)/h2*var(g)))
write.table(y, "H10Q40.txt", sep="\t")


# Simulation with 50 QTL
QTL <- 100*(1:50) #pick 50 QTL
u <- rep(0,34848) #marker effects
u[QTL] <- 1
g <- as.vector(crossprod(D2,u))
h2 <- 0.1
y <- g + rnorm(135,mean=0,sd=sqrt((1-h2)/h2*var(g)))
write.table(y, "H10Q50.txt", sep="\t")