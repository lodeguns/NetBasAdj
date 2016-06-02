#########################################################################################
#  Novel algorithms to detect oscillatory patterns in multi omic metabolic networks.
#
#  Network Based Adjacency (NBA)  - Fun. Extraction
#
#  Francesco Bardozzo, Pietro Liò, Roberto Tagliaferri
#
#########################################################################################

library(KEGGgraph)
library(igraph)

#########################################################################################
#
# function CountOccurrences - Algorithm 1
#
# In the paper this is the function countOccurences of the Algorithm 1.
# Substantially the relevant parameters present the same names described in the paper. 
# "tabellafordtx" is the name used in the list "osc.paths.list.RData" of the k-th omic
# pattern z_p.

getWeightsHOSC<-function(DG_k, locseq, org, d_NBA, psi, tabellafordtx, inout="out")
{
  
  
  graphIG<-igraph.from.graphNEL(DG_k, name = TRUE, weight = TRUE,
                                unlist.attrs = TRUE)
  
  
  count <- c(rep(0, length(locseq)))
  names(count)<- as.character((locseq))
  
  
  for(i in 1:(length(locseq)-(d_NBA-1))){
    
    fromOsc <- paste(org, as.character(locseq[i]), sep=":")
    toOsc <- paste(org, as.character(locseq[i+(d_NBA-1)]), sep=":")
    
    print(paste("", fromOsc, " ", toOsc))
    
    res0=NULL
    tryCatch( res0 <- res0 <-get.shortest.paths(graphIG, fromOsc, toOsc, mode = c(inout), weights = NULL, output=c("vpath", "epath", "both"),
                                                     predecessors = FALSE, inbound.edges = FALSE) ,
              warning = function(w) {print(paste("warning "));},
              error   = function(e) {print(paste("error ")); NaN}
              )
    
    

    
    #Vertex Sequence
    #vpath means that the vertices along the paths are reported.
    #
    if(length(res0)!=0){
    for( k in 1:length(res0$vpath)){
      
      if(length(V(graphIG)[res0$vpath[[k]]]$name)==psi){
        pr = paste("The genes on the pattern are at a distance d_NBA = ", 
                   d_NBA-1 ," and on the network close of psi = ", psi-1, " edges.", sep =" ")
        pr = paste(pr, ". Shortest path's end nodes: ", i, " and ", d_NBA-1, sep=" ")
        print(pr)
        #print(V(graphIG)[res0$vpath[[k]]]$name)
        
        
        
        #print(length(V(graphIG)[res0$vpath[[k]]]$name))
        namenodes<- V(graphIG)[res0$vpath[[k]]]$name
        
        for( z in 1:length(namenodes)){
          zx<-grep(gsub(paste(org, ":", sep=""), "", namenodes[z]), names(count))
          count[[zx]] <-count[[zx]]+1
          #print(paste(tabellafordtx[grep(gsub(paste(org, ":", sep=""), "", namenodes[z]), tabellafordtx$locus),]$CAI))
          
        }
        
      }
    }
    
    
    } 
    
    
  
  }
  
  return(count)
  
}


# This function do the extraction of DG_k
osc.get.graph <- function(pathid, org){
  g<-NULL
  print(pathid)
  gKGML <- getKGMLurl(pathid,organism=org)
  tryCatch( g <- parseKGML2Graph(gKGML, genesOnly=TRUE),
            error = function(c) {
              c$message <- paste(c$message, " pathway:",gKGML, ")", sep=" ")
              print(c)
              return(NULL)
            })
  
  return(g)
}

# This function reads the patterns from "osc.path.list.RData" and starts the
# computation of the occurences after the extraction of the metabolic networks' reactions.
osc.get.list.w <- function(osc.paths.list, oscill.org, inout)
{
  names.paths <- names(osc.paths.list)
  osc.w.list <- list()
  for(j in 1:length(osc.paths.list))
  {
    graph<- osc.get.graph(names.paths[j], oscill.org)
    if(!is.null(graph)){
      osc.w.list[[length(osc.w.list)+1]] <-  data.frame(osc.get.path.w(graph,
                                                                       oscill.org,
                                                                       osc.paths.list[[j]], 
                                                                       inout))
    } else 
      
    {
      osc.w.list[[length(osc.w.list)+1]] <- data.frame(N="NULL", name = names(osc.paths.list[j]), org = oscill.org)
      
      
    }
    
  }
  
  names(osc.w.list)<-names.paths
  
  return(osc.w.list)
  
}



# Function that setup d_NBA progressively bigger till the length
# of the pattern "tabellafordtx".
osc.get.path.w <- function(DG_k, org, tabellafordtx, inout)  {
  
  # Once fixed the distance d_NBA = 2 start the iteration 
  # on psi with the function: getWLOsc().
  
  outp <-  getWLOsc(DG_k, tabellafordtx$locus, "eco", 2, length(tabellafordtx[,1]), tabellafordtx, inout)
  
  c2 <- rbind(NULL, outp[-which(rowSums(outp[,c(1:length(tabellafordtx[,1]))]) == 0),])
  # Setup d_NBA progressively bigger till the length of the pattern "tabellafordtx".
  for(j in 3:length(tabellafordtx[,1]))
  {
    outp <- getWLOsc(DG_k, tabellafordtx$locus, "eco",j, length(tabellafordtx[,1]), tabellafordtx, inout)
    
    c2<-rbind(c2, outp[-which(rowSums(outp[,c(1:length(tabellafordtx[,1]))]) == 0),])
    
  }
  
  return(c2)
}


# Given d_NBA fixed to a value, search for a shortest path
# of a length of psi edges of distance. Psi progressively increase of 2 to the
# length  of the pattern "tabellafordtx".

getWLOsc<- function(DG_k, locseq, org, d_NBA, psi, tabellafordtx, inout)
{
  count =  getWeightsHOSC(DG_k, locseq, org, d_NBA, 2, tabellafordtx, inout)
  
  c1 <- rbind(NULL, c(count, d_NBA=d_NBA, psi=2))
  
  cx = 3
  while(cx<=psi){
    
    count = getWeightsHOSC(DG_k, locseq, org, d_NBA, cx, tabellafordtx, inout)
    c1 <- rbind(c1, c(count, d_NBA=d_NBA, psi=cx))
    cx = cx+1
    
  }
  return(c1)
}


# Get the weights for a particular value of d_NBA related to the number of edges (psi).
osc.select.w <- function(osc.list.w, d_NBA_i)
{
  osc.sum.up.w = NULL
  if(length(osc.list.w[[1]])!=3){
    osc.list.w.sm<-osc.list.w[[1]][osc.list.w[[1]]$d_NBA == d_NBA_i,]
    
    osc.sum.up.w <-osc.list.w.sm[1,-c(length(osc.list.w.sm[1,])-1, length(osc.list.w.sm[1,]))]
    osc.sum.up.w <- osc.sum.up.w / (osc.list.w.sm[1,]$psi -1 )
    
    for(j in 2 :length(osc.sum.up.w[,1]))
    {
      osc.sum.up.w1 <-osc.list.w.sm[j,-c(length(osc.list.w.sm[j,])-1, length(osc.list.w.sm[j,]))]
      osc.sum.up.w <- osc.sum.up.w + osc.sum.up.w1 / (osc.list.w.sm[j,]$psi -1 )
      
    }
  }
  return(osc.sum.up.w)
  
}

#####################################################################################
# function ComputeNBAWeights - Algorithm 1
#
# This is the function ComputeNBAWeights described in the Algorithm 1 in the paper.
# Select all the distinct values of the distance d_NBA and get the NBA weights on all
# the pattern.

osc.select.w.all <- function(osc.list.w, d_NBA_l = NULL)
{
  
  osc.sum.up.w = NULL
  if(length(osc.list.w[[1]])!=3){
    if(is.null(d_NBA_l)){
      d_NBA_l <- unique(osc.list.w[[1]]$d_NBA)
    }
    
    osc.sum.up.w <- osc.select.w(osc.list.w[1], d_NBA_l[1])/(d_NBA_l[1] - 1)
    
    
    for(j in 2:length(d_NBA_l))
    {
      osc.list.w1 <- osc.select.w(osc.list.w[1], d_NBA_l[j])
      
      if(!is.na(osc.list.w1[[1]])){
        osc.sum.up.w <- osc.sum.up.w +  osc.list.w1/(d_NBA_l[j] - 1)
        
      }
    }
  }
  return(osc.sum.up.w)
  
}







#########################################################################################

#########################################################################################

#Only an idea, not used in this work.
getWHOsc<- function(DG_k, locseq, org, d_NBA, psi, tabellafordtx)
{
  
  count =  getWeightsHOSC(DG_k, locseq, org, 2, psi, tabellafordtx)
  
  cx = 3
  while(cx<=d_NBA)
  {
    
    for(i in 1:length(count))
    {
      if(count[[i]] != 0 )
      { 
        count[[i]] = count[[i]] + count[[i]]
      }
    }
    
    count = count + getWeightsHOSC(DG_k, locseq, org,cx,psi, tabellafordtx)
    
    cx = cx+1
    
  }
  return(count)
}
