#Diagnostic check
library(rwty)


my.trees <- load.multi("/Users/bguinet/Desktop/Cynipoidea_paper/Mr_bayes/", format = "mb")


colnames(my.trees$ptable)
makeplot.param(my.trees, burnin = 0, "LnL")

my.trees.rwty <- analyze.rwty(my.trees, burnin=5000)
makeplot.all.params(my.trees, burnin=0) #
approx.ess <- topological.approx.ess(my.trees, burnin = 50)



#PLot tree datation 

library(MCMCtreeR)

directory.mb <- '/Users/bguinet/Desktop/Cynipoidea_paper/Mr_bayes/Concatenated_sequences_main_datation.nxs4.con.tre'
MCMC.tree.plot(analysis.type='mrbayes',
    directory.files=directory.mb, cex.tips=0.33,
    plot.type='phylogram', lwd.bar=2, add.time.scale=FALSE,
    node.method='bar', col.age='navy')

MCMC.tree.plot(analysis.type='mrbayes',
               directory.files=directory.mb,
               time.correction = 100, plot.type = "distributions", cex.age = 0.4, 
                cex.labels = 0.5, relative.height = 0.08, col.tree = "grey40", 
                scale.res = c("Eon", "Period"), no.margin = TRUE, label.offset = 4, 
                 density.col = "#00000050", density.border.col = "#00000080")

phy <- readMCMCtree(MCMCtree.phy, from.file = FALSE)


MCMC.tree.plot(phy, cex.tips = 0.2, time.correction = 100, 
               scale.res = c("Eon", "Period"), plot.type = "phylogram", cex.age = 0.6, cex.labels = 0.6, 
               relative.height = 0.08, col.tree = "grey40", label.offset = 4, 
               node.method = "none", no.margin = TRUE)




