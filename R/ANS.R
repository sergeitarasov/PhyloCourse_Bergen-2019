library("corHMM")
library("phytools")

setwd("~/Documents/Courses/Bergen_2019/PhyloCourse_Bergen_2019")
source('write_nexus_RB.R')

#   ____________________________________________________________________________
#   Simulate Tree and Traits                                               ####

# simulate tree using pure birth process
tree<-pbtree(n=200, scale=100, b=1, d=0)
plot(tree)

# simulate character evolution

# make rate matrix Q
Q <- matrix(
  c(
    -0.03, 0.03,
    0.1, -0.1
    ), 2,2, byrow = T)

Q

# simulate character evolution on tree using Q
hist <- sim.history(tree, Q, nsim=1)
plot(hist)


#   ____________________________________________________________________________
#   Ancestral Character State Recon                                         ####





##  ............................................................................
##  Inference using Symmetrical Model                                       ####

# one-rate symmetrical model

Q.Sym <- matrix(
  c(
    NA, 1,
    1, NA
  ), 2,2, byrow = T)

Q.Sym

# Inference
taxa <- cbind(hist$tip.label, hist$states)
Recon_Q.Sim <- rayDISC(hist, taxa, rate.mat=Q.Sym, node.states="marginal", 
        model="ARD", root.p="maddfitz")

# infered rate matrix
Recon_Q.Sim

# plot ANCE
plotRECON(tree, Recon_Q.Sim$states, piecolors=c('black', 'red'), title="1-rate Model")


##  ............................................................................
##  Inference using Asymmetrical Model                                     ####

# 2-rate symmetrical model

Q.Asym <- matrix(
  c(
    NA, 1,
    2, NA
  ), 2,2, byrow = T)

Q.Asym

# Inference
taxa <- cbind(hist$tip.label, hist$states)
Recon_Q.Asim <- rayDISC(hist, taxa, rate.mat=Q.Asym, node.states="marginal", 
                       model="ARD", root.p="maddfitz")

# infered rate matrix
Recon_Q.Asim

# plot ANCE
plotRECON(tree, Recon_Q.Asim$states, piecolors=c('black', 'red'), title="1-rate Model")



#   ____________________________________________________________________________
#   Compare Results                                                         ####

Q
Recon_Q.Sim
Recon_Q.Asim

Recon_Q.Sim$AIC-Recon_Q.Asim$AIC



#   ____________________________________________________________________________
#   Save simulated results to Run RevBayes                                  ####

# write traits
otu <- dplyr::as_tibble(taxa)
write_nexus_RevBayes(otu, file = "Sim_traits.nex")

#write tree
write.tree(tree, 'Sim_tree.tre')


