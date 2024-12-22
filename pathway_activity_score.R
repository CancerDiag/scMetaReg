library(scde)
library('pagoda2')

PAS_pagoda <- function(counts, gSets, log.scale = F, nrandom = 100, proj_name){
    nPcs = min(round(ncol(counts)/5),5)
    p2 = Pagoda2$new(counts, log.scale=F)
    p2$adjustVariance(plot=F)
    p2$calculatePcaReduction(nPcs = nPcs,use.odgenes=F,fastpath=F)
    
    path_names = c()
    env = new.env(parent=globalenv())
    invisible(lapply(1:length(gSets),function(i) {
      genes = intersect(gSets[[i]],rownames(counts))
      #name = paste0(names(gSets[i]),i)
      name = names(gSets[i])
      if(length(genes)>=3){
        assign(name, genes, envir = env)
        path_names = c(path_names, name)
      }
    }))
    

    p2$testPathwayOverdispersion(setenv = env, verbose = T,
                                 recalculate.pca = T,
                                 min.pathway.size = 1, n.randomizations = nrandom)

    path_names = names(p2$misc$pwpca)
    score = matrix(NA,nrow=length(path_names),ncol=ncol(counts))
    rownames(score) = path_names
    colnames(score) = colnames(counts)
    for(i in 1:length(p2$misc$pwpca)){
      if(!is.null(p2$misc$pwpca[[i]][1]$xp$scores)){
        score[i,] = as.numeric(p2$misc$pwpca[[i]][1]$xp$scores)
        }
    }
    score <- score[complete.cases(score),]
    saveRDS(p2, file.path(outDir,paste0(proj_name, 'pagoda2obj.rds')))
    return(score)
   }