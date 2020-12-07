library(circlize)
library(ape)
library(dendextend)
library(phytools)
library(RColorBrewer)

#### Data import and processing ####
# quickly get expression levels for our rosette leaf 1
rosetteExpr = read.table("RNAseq.transcript.abundance.txt",
                         stringsAsFactors = F, header = T)
rosetteExpr = rosetteExpr[,c(1,46)] # we know the data we need is only in col 46

# import out tree we got from Clustal
clustTree = read.tree(file = "mouse_outgroup_tree.txt")
# reroot the tree to the mouse outgroup
clustTree = root(clustTree, 
     outgroup = "sp|Q8CJ12|AGRG2_MOUSE",resolve.root = TRUE)
# make it ultrametric, all leaves are equidistant from root
hc = as.hclust(force.ultrametric(clustTree, method = "extend"))

# from out tree, extract labels and length, convert to R dendrogram object
labels = hc$labels
ct = cutree(hc, 6)
n = length(labels)
dend = as.dendrogram(hc)

# create matrix of our expression abundances
mat = matrix(rep(0, 98), nrow= n, ncol=1) # we have 98 values to go through
for (i in 2:98) { # R is cursed and indexes start at 1
  # match our label name with our abundance data
  mat[i] = rosetteExpr[which(rosetteExpr$transcript_id == labels[i]), 2]
}
mat[1] = 0 # set our mouse to empty
# the function that generates the various shades of black -> green
col_fun = colorRamp2(c(0, max(rosetteExpr$DV.LF.1)), c("black", "green"))

#### plotting ####

# create path to save in:
png(file="A.thaliana_circlized.png", width=1200, height=1200) # hi-def

# initialize the circle plot
circos.par(cell.padding = c(0,0,0,0))
circos.initialize("a", xlim = c(0, n))
# circlize plots work by plotting on track at a time, outside in. this is labels
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.3,
             # this function iterates through all of the labels, and uses
             # circos methods to place them according to the below parameters
             # all labels are set to black
             panel.fun = function(x, y){
               for(i in seq_len(n)) {
                 circos.text(i - 0.5, 0, labels[i], adj = c(0, 0.5),
                             facing = "clockwise", niceFacing = TRUE,
                             col = 1, cex = 0.85)}})

# this plots the heatmap, which is essentially a row of 98 colored rectangles
# luckily for us, circlize wraps it into a circle
circos.track(ylim = c(0, 0.5), bg.border = NA, 
             panel.fun = function(x, y) {
               # array of colors for each rectangle
               col_mat = col_fun(mat)
               # change first rectangle to white (mouse)
               col_mat[1] = "#FFFFFFFF"
               # first four arguments are xleft, ybottom, xright, ytop
               # since our tree is evenly spaced, we can just pass four arrays,
               # and the colors will get plotted in the order of the labels
               circos.rect(0:97, rep(0, 98), 1:98, rep(0.5, 98), 
                          border = col_mat[1:98], col = col_mat[1:98])})

# finally, plot the dendrogram using the built in function
dend_height = attr(dend, "height")
circos.track(ylim = c(0, dend_height), bg.border = NA, track.height = 0.4,
             panel.fun = function(x, y) {
               circos.dendrogram(dend)})

# add title and legend, italics for biological consistency ;)
title(expression("RNASeq Transcript Abundance of 100 "*italic(A.~thaliana)*" genes"),
      cex.main = 3, line = -2)
legend(x="bottomright", legend=c(0, max(rosetteExpr$DV.LF.1)), fill=c("black", "green"),cex = 2)

circos.clear()
dev.off()
