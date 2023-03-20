#' Convenience Function to set up a cluster used for multiprocessing
#'
#' @param ncores Numbers of processes the cluster should use
#' @param cluster_type "FORK" or "PSOCK". Fork is faster but doesn't work on
#'                     Windows
#' @param debug TRUE if the cluster log should be written to file
#' @param outfile Path where the cluster log is written to in debug mode.
#'
#' @return
#' @export
#'
#' @examples
cluster_setup <- function(ncores = 2,
                          cluster_type = "FORK",
                          debug = FALSE,
                          outfile = NULL) {
  if (debug) {
    cluster <- parallel::makeCluster(ncores,
                                     type = cluster_type,
                                     outfile = outfile)
  } else {
    cluster <- parallel::makeCluster(ncores,
                                     type = cluster_type
    )
  }
  # Multiprocessing setup
  doParallel::registerDoParallel(cluster)
  parallel::clusterEvalQ(cl = cluster, {
    options(mc.cores = 1)
  })
}
