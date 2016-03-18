#########################################################################################
#  Novel algorithms to detect oscillatory patterns in multi omic metabolic networks.
#
#  Network Based Adjacency (NBA)  - Fun. Extraction
#
#  Francesco Bardozzo, Pietro Liò, Roberto Tagliaferri
#
#########################################################################################

#Load DataSet
#https://github.com/lodeguns/NetBasAdj/blob/master/osc.paths.list.RData


load('./osc.paths.list.Rdata')

#Load Functions
#https://tonybreyal.wordpress.com/2011/11/24/source_https-sourcing-an-r-script-from-github/
source_https <- function(url, ...) {
  require(RCurl)
  sapply(c(url, ...), function(u) {
    eval(parse(text = getURL(u, followlocation = TRUE, 
    cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), envir = .GlobalEnv)
  })
}

source_https("https://github.com/lodeguns/NetBasAdj/raw/master/NetBasAdj.R")


########################################################################
#Run the example for ingoing and outgoing degree. At the index 1 of the
#dataset "osc.paths.list" are tabulated the omics values for the metabolic
#network Glycolysis of Escherichia coli K-12 MG1655.

oscill.org = "eco"

glyco.w.in  <- osc.get.list.w(osc.paths.list[1], oscill.org, "in")
osc.select.w.all(glyco.w.in)

glyco.w.out <- osc.get.list.w(osc.paths.list[1], oscill.org, "out")

#Selected Adj.s
osc.select.w.all(glyco.w.out, c(2,3,4))

#All
osc.select.w.all(glyco.w.out)

#All the datas are publicly available and comes from KEGG, NCBI, PaXDB.
########################################################################



