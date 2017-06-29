#=============================================================================#
# topGO.R           													      #
#																			  #
# Version: 1.0   															  #
# Date: 6 Dec 2016											             	  #
# Author: Michiel Adriaens, PhD; MaCSBio, Maastricht University               #
# History:																	  #
#  1.0: Creation															  #
#                                                                             #
#=============================================================================#

genesOfInterest <- ... # vector of Ensembl Gene IDs, e.g. your set of diff. expressed genes			
library(topGO)
library(xlsx)
library(org.Hs.eg.db) # Homo Sapiens, also available for e.g. rat: org.Rn.eg.db
for (j in c("BP", "CC", "MF")) {			
	message(j)		
	xx <- annFUN.org(j, mapping = "org.Rn.eg.db", ID = "ensembl") # Also "entrez", "symbol" and other options
	allGenes <- unique(unlist(xx))		
	allGenesFactor <- factor(as.integer(allGenes %in% genesOfInterest))		
	names(allGenesFactor) <- allGenes		
	GOdata <- new("topGOdata",		
			ontology = j,
			allGenes = allGenesFactor,
			nodeSize = 5, # recommended setting
			annot = annFUN.org,
			mapping = "org.Rn.eg.db",
			ID = "ensembl")
	goRes <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher") # recommended setting		
	allRes <- GenTable(GOdata, Pvalue = goRes, 		
			topNodes = length(goRes@score),
			orderBy = "Pvalue", ranksOf = "Pvalue",
			numChar = 1E9)
	allRes <- allRes[1:max(which(as.numeric(allRes[,6]) < 0.05)), ]	# only save with p < 0.05	
	write.xlsx(allRes,		
			file = "topGOres.xlsx", 
			sheet = j,
			row.names = FALSE,
			append = TRUE)
}			

