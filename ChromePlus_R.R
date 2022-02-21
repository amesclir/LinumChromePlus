#install.packages("Rtools")
#install.packages("devtools")
library(devtools)
#install_github('coleoguy/chromePlus')
#force=TRUE
library(chromePlus)
library(diversitree)

##Chromosome number and stylar states##
#Read the tree

tree <- read.tree("Linum.tree")
tree
tree <- drop.tip(tree, "L_oligophyllum")

#Read the data

dat1 <- read.csv("Linum_StylarDataset.csv", sep=",")


#Convert data appropiate for analysis. 6 y 36 son el min y el max n de cromosomas

setdiff(tree$tip.label, dat1$Name)
setdiff(dat1$Name,tree$tip.label)

#?datatoMatrix
dat.mat <- datatoMatrix(x=dat1, range= c(6,36), hyper=T)


# make the basic likelihood function
?make.mkn
lik <- make.mkn(tree = tree, states= dat.mat, k=ncol(dat.mat), strict = F, control = list(method="ode"))

# constrain the likelihood function to a biologically realistic design
#?constrainMkn
con.lik <- constrainMkn(data = dat.mat,
                        lik = lik,
                        polyploidy = F,
                        hyper = T,
                        constrain = list(drop.demi = T,
                                         drop.poly= F))

# lets make sure we have the parameters we expect
#?argnames
argnames(lik)
argnames(con.lik)

#Una vez que tenemos lik y con.lik, sigo utilizando diversitree
#Crear modelos ML con "lik" y con "con.lik"

fit <- find.mle(lik, c(0.473, 0.54), method="optim")










##Chromosome number and distribution##

#Read the tree

tree <- read.tree("Linum.tree")
tree
tree <- drop.tip(tree, "L_oligophyllum")

#Read the data

dat2 <- read.csv("Linum_distributionDataset.csv", sep=",")


#Convert data appropiate for analysis. 6 y 36 son el min y el max n de cromosomas

setdiff(tree$tip.label, dat1$Name)
setdiff(dat1$Name,tree$tip.label)

?datatoMatrix
dat.mat2 <- datatoMatrix(x=dat2, range= c(6,36), hyper=T)


# make the basic likelihood function

?is.ultrametric
is.ultrametric(tree)

?make.mkn

lik <- make.mkn(tree = tree, states= dat.mat, k=ncol(dat.mat), strict = F, control = list(method="ode"))


# constrain the likelihood function to a biologically realistic design

?constrainMkn

con.lik <- constrainMkn(data = dat.mat,
                        lik = lik,
                        polyploidy = F,
                        hyper = T,
                        constrain = list(drop.demi = T,
                                         drop.poly= F))


# lets make sure we have the parameters we expect

?argnames
argnames(con.lik)



