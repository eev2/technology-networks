folder = 'classes' # 'classes' or 'subclasses'

if (!require(igraph)) {
    install.packages("igraph", repos = "http://www.stats.bris.ac.uk/R/")
    library(igraph)
}

folder0 = paste0('../data/', folder)
alllb = scan(paste0(folder0, '/ipc_class_description_4.tsv'),
             character())
lb = alllb[!grepl("^[A-Z]99.*", alllb)]
N = length(lb)
nyr = 2011 - 1979

sedgefile = list.files(pattern = paste0(folder, "_sedge\\.csv.*"))
if(grepl(".*\\.zip", sedgefile))
  sedgefile = unz(sedgefile, sub(".zip", "", sedgefile, fixed = TRUE))
Adj = array(scan(sedgefile), c(N, N, nyr))
nodecore = matrix(scan(paste0(folder, "_nodecore.csv"), sep = ","),
                  N, byrow = TRUE)
typecols = c("white", "blue", "red")

secs = sub("^([A-Z])[0-9A-Z]+$", "\\1", lb)
seclb = unique(secs)
secidx = lapply(seclb, function(a) which(secs == a))
seccnt = table(secs)

seccols = c(A = "#0000FF", B = "#0080FF", C = "#00FFFF", D = "#80FF80", E = "#FFFF00", F = "#FF8000", G = "#FF0000", H = "#800000")

tt = seq(1,32,6)

pdf(paste0(folder, "_network_plot.pdf"), width = 8, height = 16)
oldpar = par(mar = c(.1,.1,2,.1), mfrow=c(4,2))
for (j in 1:nyr) {
    net = graph_from_adjacency_matrix(Adj[, , j], mode = "directed",
                    weighted = NULL, diag = FALSE)
    plot(net, vertex.size = 4, edge.arrow.size = .2,
         edge.arrow.width = .8, edge.curved = 0.1,
         vertex.color = adjustcolor(typecols,.5)[nodecore[,j]+1],
         #vertex.frame.color = seccols[secs],
         vertex.label = NA,
         margin = 0, main = j + 1979,
         layout = layout_with_fr)
         # layout = layout_with_lgl(net,root = 121)
         # mark.groups = secidx
         # mark.col=adjustcolor(seccols,.05)
}
dev.off()

pdf(paste0(folder, "_network_plot_sel.pdf"), width = 6, height = 9)
oldpar = par(mar = c(.1,.1,2,.1), mfrow=c(3,2))
for (j in tt) {
    AA = Adj[, , j]
    AA = crossprod(AA)*AA+AA
    net = graph_from_adjacency_matrix(AA, mode = "directed",
                    weighted = NULL, diag = FALSE)
    plot(net, vertex.size = 4, edge.arrow.size = .2,
         edge.arrow.width = .8, edge.curved = 0.1,
         vertex.color = adjustcolor(typecols,.5)[nodecore[,j]+1],
         #vertex.frame.color = seccols[secs],
         vertex.label = NA,
         margin = 0, main = j + 1979,
         layout = layout_with_fr)
         # layout = layout_with_lgl(net,root = 121)
         # mark.groups = secidx
         # mark.col=adjustcolor(seccols,.05)
}
dev.off()

pdf(paste0(folder, "_network_plot_acsonly_sel.pdf"), width = 6, height = 9)
oldpar = par(mar = c(.1,.1,2,.1), mfrow=c(3,2))
for (j in tt) {
    AA = Adj[, , j]
    jsel = nodecore[,j] > 0
    AA = AA[jsel, jsel]
    net = graph_from_adjacency_matrix(AA, mode = "directed",
                    weighted = NULL, diag = FALSE)
    plot(net, vertex.size = 4, edge.arrow.size = .2,
         edge.arrow.width = .8, edge.curved = 0.1,
         vertex.color = adjustcolor(typecols,.5)[nodecore[jsel,j]+1],
         #vertex.frame.color = seccols[secs],
         ##vertex.label = lb[jsel], vertex.label.cex = .7,
         vertex.label = NA,
         margin = 0, main = j + 1979,
         layout = layout_with_graphopt)
         # layout = layout_with_lgl(net,root = 121)
         # mark.groups = secidx
         # mark.col=adjustcolor(seccols,.05)
}
dev.off()


pdf(paste0(folder, "_network_heat_sel.pdf"), width = 7, height = 7)
for (j in tt) {
  AA = Adj[,,j]
  diag(AA) = 0
  id <- 0
  for (i in 1:8) {
    ic = id + 1; id <- id + seccnt[i]
    ib = 0
    for (k in 1:8) {
      ia = ib + 1; ib = ib + seccnt[k]
      if (!((i-k)%%2)) AA[ia:ib, ic:id] = AA[ia:ib, ic:id] + .5
    }
  }
  AA[AA > 1] <- 1
  heatmap(AA, Rowv = NA, Colv = NA, col = c("honeydew", "lavender", "black"),
          scale="none",
          margins = c(3,3), ColSideColors = rep(seccols, seccnt),
          RowSideColors = rep(seccols, seccnt), revC = TRUE,
          zlim = range(unlist(Adj)), main = j + 1979)
}
dev.off()
