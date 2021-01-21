######################################################
### PLS regression with weighted pivot coordinates ###
######################################################

library(compositions)
library(pls)

source("weightPivotCoord.R")


#load pre-processed data with no 0
load("Data.Rdata")
head(Data)
X1 = as.matrix(Data[, 2:128]) # Choose the appropriate columns.
y1 = log(Data[, 1]) # Choose the appropriate column.
n = length(y1)
D = ncol(X1)

set.seed(1) # for reproducibility

#centering the variables
y = as.vector(scale(y1, scale = F))
X = as.matrix(scale(acomp(X1), scale = F))

# WPC
wpc = weightPivotCoord(X, option = "cor", yvar = y)
Z = wpc$WPC
Zlc = wpc$w
W = Z
Wlc = Zlc
for (i in 2:D){
  wpc = weightPivotCoord(X, pivotvar=i, option="cor", yvar=y)
  Z = wpc$WPC
  W = cbind(W, Z)
  Zlc = wpc$w
  Wlc = cbind(Wlc, Zlc)
}


# CV RMSEP and R^2
misclass=mvr(y~Z, ncomp = 10, method="simpls", validation="CV")
pls::RMSEP(misclass, estimate="CV", intercept=FALSE)
pls::R2(misclass, estimate="CV", intercept=FALSE)

# optimal number of PLS components based on the randomization test approach
comp = selectNcomp(misclass)


# bootstrap
kk = 1000
BE = matrix(rep(NA, kk*D), ncol = D)
for(k in 1:kk){
  ind = sample(1:n, n, replace = TRUE)  
  testX1 = X1[ind, ]
  testy1 = y1[ind]
  testX = as.matrix(scale(acomp(testX1), scale = FALSE))
  testy = as.vector(scale(testy1, scale = FALSE)) 
  for (i in 1:D){
    testx1 = cbind(testX[,i], testX[,-i])
    WW = Wlc[, ((i-1)*(D-1)+1):(i*(D-1))]
    testZ = log(testx1)%*%WW
    resst = mvr(testy~testZ, ncomp = comp, method = "simpls")
    tbst = as.vector(coef(resst))
    BE[k, i] = tbst[1]
  }
  print(k)
}

# significance of bootstrap standardized coefficients
beta2 = apply(BE, 2, mean)
sdbeta = apply(BE, 2, sd)
alpha2 = -qnorm(0.025/D)
a = beta2/sdbeta

signary = rep(0, D)
for(i in 1:D){
  if((a[i]) > alpha2){signary[i] = 1}
  if((a[i]) < -alpha2){signary[i] = -1}
}
signary
names(signary) = colnames(X1)
signary



#############################################################################################
### PLS biplot

Z = W[, (1:(D-1))]
resst = mvr(y~Z, ncomp = comp, method = "simpls") 

G = matrix(as.numeric(resst$scores), n, comp)[, 1:2]
H = matrix(as.numeric(resst$loadings), D-1, comp)[, 1:2]
H = rbind(H, c(0,0))

for (i in 2:D){
  Z = W[, ((i-1)*(D-1)+1):(i*(D-1))]
  resst = mvr(y~Z, ncomp = comp, method = "simpls") 
  h = as.numeric(resst$loadings)[c(1, D)]
  H[i,] = h
}
H = 1*H # Scale loadings so that the arrows (points, respectively) are more visible in the biplot?

colnames(G) = c("Comp1", "Comp2")
colnames(H) = c("Comp1", "Comp2")

G = as.data.frame(G)
H = as.data.frame(H)
H$Variable = colnames(X1)
H$Angle = ((180/pi)*atan(H$Comp2/H$Comp1))
H$Adj = (1-1.5*sign(H$Comp1))/2

pdf("PLSbiplot.pdf")
g = ggplot(G, aes(x = Comp1, y = Comp2)) +
  geom_point() +
  geom_segment(data = H, aes(x = 0, y = 0, xend = Comp1, yend = Comp2),
               arrow = arrow(length = unit(1/2, "picas")), colour = c("blue","grey","red")[factor(signary)]) +
  geom_text(data = H, aes(label = Variable, x = Comp1, y = Comp2, angle = Angle, hjust = Adj), 
            size = 3, colour = c("blue","grey","red")[factor(signary)]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  xlab("PLS comp. 1") + ylab("PLS comp. 2") +
  theme_classic() + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_text(size = 15, vjust = 2, face = "bold"))
g + coord_fixed(ratio = 1) # Change ratio between y-axis and x-axis scale?
dev.off()
