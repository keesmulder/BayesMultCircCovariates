library(parallel)

source('Simulation/simulationStudyCircGLM.R')

# DOES NOT WORK
parallelSimStudyCircGLM <- function(predDesigns = predDesigns, ...) {

  n_cores <- detectCores() - 1

  # Register a parallel computing cluster.
  cl <- makePSOCKcluster(n_cores)

  params <- list(...)

  browser()

  # Export function from the session.
  clusterExport(cl, params)

  allSimRes <- parLapply(cl, predDesigns, function(pddes) {

        thisPdDesSimRes <- do.call(simStudyCircGLM,
                                   c(list(predDesigns = pddes), params))
        thisPdDesSimRes
  })

  stopCluster(cl)

  return(allSimRes)

}
