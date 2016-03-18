
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

osc.get.path.w <- function(toyGraph, org, tabellafordtx, inout)  {
  
  outp <-  getWLOsc(toyGraph, tabellafordtx$locus, "eco", 2, length(tabellafordtx[,1]), tabellafordtx, inout)
  
  c2 <- rbind(NULL, outp[-which(rowSums(outp[,c(1:length(tabellafordtx[,1]))]) == 0),])
  
  for(j in 3:length(tabellafordtx[,1]))
  {
    outp <- getWLOsc(toyGraph, tabellafordtx$locus, "eco",j, length(tabellafordtx[,1]), tabellafordtx, inout)
    
    c2<-rbind(c2, outp[-which(rowSums(outp[,c(1:length(tabellafordtx[,1]))]) == 0),])
    
  }
  
  return(c2)
}



getWLOsc<- function(toyGraph, locseq, org, hm, dd, tabellafordtx, inout)
{
  count =  getWeightsHOSC(toyGraph, locseq, org, hm, 2, tabellafordtx, inout)
  
  c1 <- rbind(NULL, c(count, hm=hm, cx=2))
  
  cx = 3
  while(cx<=dd){

    count = getWeightsHOSC(toyGraph, locseq, org, hm, cx, tabellafordtx, inout)
    
    
    c1 <- rbind(c1, c(count, hm=hm, cx=cx))
    
    cx = cx+1
    
  }
  return(c1)
}



getWHOsc<- function(toyGraph, locseq, org, hm, dd, tabellafordtx)
{
  
  count =  getWeightsHOSC(toyGraph, locseq, org, 2, dd, tabellafordtx)
  
  cx = 3
  while(cx<=hm)
  {
    
    for(i in 1:length(count))
    {
      if(count[[i]] != 0 )
      { 
        count[[i]] = count[[i]] + count[[i]]
      }
    }
    
    count = count + getWeightsHOSC(toyGraph, locseq, org,cx,dd, tabellafordtx)
    
    cx = cx+1
    
  }
  return(count)
}






getWeightsHOSC<-function(toyGraph, locseq, org, h, dd, tabellafordtx, inout="out")
{
  
  
  graphIG<-igraph.from.graphNEL(toyGraph, name = TRUE, weight = TRUE,
                                unlist.attrs = TRUE)
  
  
  count <- c(rep(0, length(locseq)))
  names(count)<- as.character((locseq))
  
  
  for(i in 1:(length(locseq)-(h-1))){
    
    fromOsc <- paste(org, as.character(locseq[i]), sep=":")
    toOsc <- paste(org, as.character(locseq[i+(h-1)]), sep=":")
    
    res0 <- get.shortest.paths(graphIG, from=fromOsc , to=toOsc, #graphIG, v=V(graphIG), to=V(graphIG),
                               mode = inout, 
                               weights = NULL, 
                               output=c("vpath"),
                               predecessors = TRUE,  #potrebbero servirmi.
                               inbound.edges = TRUE) 
    #Vertex Sequence
    #vpath means that the vertices along the paths are reported.
    #
    for( k in 1:length(res0$vpath)){
      
      if(length(V(graphIG)[res0$vpath[[k]]]$name)==dd){
        pr = paste("Genes on an ordered sequence near of positions n = ", 
                   h-1 ," and on the network near of ", dd-1, " edges.", sep =" ")
        pr = paste(pr, " indexes: ", i, " e ", h-1, sep=" ")
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
  
  
  return(count)
  
}


#########################################################################################
#  Multi Omic Oscillations Networks
#
#  Network Based Adjacency (NBA)  - Function for short memory selection weights. 
#
#  Francesco Bardozzo, Pietro Liò, Roberto Tagliaferri
#
#########################################################################################


osc.select.w <- function(osc.list.w, short.memory)
{
  osc.sum.up.w = NULL
  if(length(osc.list.w[[1]])!=3){
    osc.list.w.sm<-osc.list.w[[1]][osc.list.w[[1]]$hm == short.memory,]
    
    osc.sum.up.w <-osc.list.w.sm[1,-c(length(osc.list.w.sm[1,])-1, length(osc.list.w.sm[1,]))]
    osc.sum.up.w <- osc.sum.up.w / (osc.list.w.sm[1,]$cx -1 )
    
    for(j in 2 :length(osc.sum.up.w[,1]))
    {
      osc.sum.up.w1 <-osc.list.w.sm[j,-c(length(osc.list.w.sm[j,])-1, length(osc.list.w.sm[j,]))]
      osc.sum.up.w <- osc.sum.up.w + osc.sum.up.w1 / (osc.list.w.sm[j,]$cx -1 )
      
    }
  }
  return(osc.sum.up.w)
  
}

osc.select.w.all <- function(osc.list.w, short.memory = NULL)
{
  
  osc.sum.up.w = NULL
  if(length(osc.list.w[[1]])!=3){
    if(is.null(short.memory)){
      short.memory <- unique(osc.list.w[[1]]$hm)
    }
    
    osc.sum.up.w <- osc.select.w(osc.list.w[1], short.memory[1])/(short.memory[1] - 1)
    
    
    for(j in 2:length(short.memory))
    {
      osc.list.w1 <- osc.select.w(osc.list.w[1], short.memory[j])
      
      if(!is.na(osc.list.w1[[1]])){
        osc.sum.up.w <- osc.sum.up.w +  osc.list.w1/(short.memory[j] - 1)
        
      }
    }
  }
  return(osc.sum.up.w)
  
}




#########################################################################################
#########################################################################################

