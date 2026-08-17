#' performs SFS to identify combination of input variables maximizing a
#' criterion
#'
#' @param input_raster SpatRaster or list of SpatRaster
#' @param obs_vect SpatVector or SpatVectorCollection
#' @param obs2optimize numeric .list of ground obs diversity metrics
#' corresponding to obs_vect.
#' Expected values: richness, shannon, simpson, hill
#' @param obs_criterion character. richness, shannon, simpson, hill
#' @param corrMethod character. select between method available for cor.test
#' @param input_mask SpatRaster corresponding to mask
#' @param Hill_order numeric.
#' @param nb_clusters numeric.
#' @param min_sun numeric.
#' @param nb_pix numeric.
#' @param nb_iter numeric.
#' @param pcelim numeric.
#' @param verbose boolean.
#' @param nbWorkers numeric.
#' @param algorithm character.
#' @param nbCPU numeric.
#'
#' @return list including performances (correlation) of SFS with additional
#' features and assessed diversity metrics corresponding to each step
#' @importFrom doFuture registerDoFuture
#' @importFrom future plan multisession sequential
#' @importFrom foreach foreach %dopar%
#' @importFrom dplyr group_split
#' @importFrom vegan mantel
#' @importFrom progressr progressor handlers with_progress
#' @importFrom stats cor.test
#' @importFrom parallel makeCluster stopCluster
#'
#' @export

biodivMapR_sfs <- function(input_raster, obs_vect, obs2optimize,
                           obs_criterion = 'shannon', corrMethod = 'spearman',
                           input_mask = NULL, Hill_order = 1, nb_clusters = 50,
                           min_sun = 0.25, nb_pix = 1e5, nb_iter = 10,
                           pcelim = 0.02, verbose = TRUE, nbWorkers = 1,
                           algorithm = "Hartigan-Wong", nbCPU = 1){

  FullListIndices <- c('richness', 'shannon', 'simpson', 'hill',
                       'bray', 'brayturn', 'jaccard', 'jaccardturn', 'sorensen', 'simpson_diss',
                       'FRic', 'FEve', 'FDiv', 'FDis', 'FRaoq')
  #### Which diversity metrics should be computed?
  alphamet <- c('richness', 'shannon', 'simpson', 'hill')
  betamet <- c('bray', 'brayturn', 'jaccard', 'jaccardturn', 'sorensen', 'simpson_diss')
  fmet <- c('FRic', 'FEve', 'FDiv', 'FDis', 'FRaoq')
  # if computation of functional metrics required
  alpha_metrics <- alphamet[which(alphamet %in% names(obs2optimize))]
  if (length(alpha_metrics)==0)
    alpha_metrics <- NULL
  # computation of beta diversity required?
  beta_metrics <- betamet[which(betamet %in% names(obs2optimize))]
  if (length(beta_metrics)==0)
    getBeta <- FALSE
  if (length(beta_metrics)>0)
    getBeta <- TRUE
  functionalmetrics <- fmet[which(fmet %in% names(obs2optimize))]
  if (length(functionalmetrics)==0)
    fd_metrics <- NULL
  if (length(functionalmetrics)>0)
    fd_metrics <- functionalmetrics

  # 1- prepare for kmeans over the full spatial extent
  message('sampling pixels to compute spectral species')
  rast_sample <- rast_sample_sfs(input_raster = input_raster,
                                 input_mask = input_mask, nb_pix = nb_pix,
                                 nb_iter = nb_iter)

  message('plot extraction')
  extracted_sfs <- extract_from_rast_sfs(input_raster = input_raster,
                                         input_mask = input_mask,
                                         min_sun = min_sun,
                                         obs_vect = obs_vect)
  rast_val <- extracted_sfs$rast_val
  Attributes <- extracted_sfs$Attributes
  nbPlots_total <- extracted_sfs$nbPlots_total
  IDplot <- extracted_sfs$IDplot

  # 3- perform SFS
  message('perform SFS')
  # initialize lists and vars
  Corr_criterion <- EvolCorr <-   CorrSFS <- AssessSFS <- list()
  SelectedVars <- c()
  for (idx in FullListIndices)
    EvolCorr[[idx]] <- c()

  AllVars <- names(rast_sample)
  nb_pcs_to_keep <- length(AllVars)

  # initialize parallel computing
  registerDoFuture()
  cl <- parallel::makeCluster(nbWorkers)
  with(plan('cluster', workers = cl), local = TRUE)
  progressr::handlers(global = TRUE)
  progressr::handlers('cli')
  progressr::with_progress({
    p <- progressr::progressor(steps = nb_pcs_to_keep)
    # increment number of variables
    for (nbvars2select in seq_len(nb_pcs_to_keep)){
      NumVar_list <- as.list(seq_len(length(AllVars)))
      numvar <- win_ID <- NULL
      subfeatures_SFS <- function() {
        foreach(numvar = NumVar_list,  .packages='biodivMapR') %dopar% {
          CorrVal <- Assess <- list()
          SelFeat_tmp <- c(SelectedVars,AllVars[[numvar]])
          # compute kmeans over selected features
          Kmeans_info <- get_kmeans(rast_sample = rast_sample[SelFeat_tmp],
                                    nb_iter = nb_iter, algorithm = algorithm,
                                    nb_clusters = nb_clusters,
                                    nbCPU = 1, progressbar = FALSE)
          SSValid <- get_spectralSpecies(inputdata = rast_val[SelFeat_tmp],
                                         Kmeans_info = Kmeans_info)
          SSValid$win_ID <- IDplot

          # alpha diversity metrics
          alpha_sfs <- alpha_metrics_sfs(alpha_metrics = alpha_metrics,
                                         obs2optimize = obs2optimize,
                                         SSValid = SSValid, corrMethod = corrMethod,
                                         nb_clusters = nb_clusters, Hill_order = Hill_order,
                                         pcelim = pcelim, nbPlots_total = nbPlots_total)

          # beta diversity metrics
          beta_sfs <- beta_metrics_sfs(obs2optimize = obs2optimize,
                                       SSValid = SSValid,
                                       beta_metrics = beta_metrics,
                                       nb_clusters = nb_clusters,
                                       pcelim = pcelim,
                                       nb_iter = nb_iter,
                                       nbPlots_total = nbPlots_total,
                                       Attributes = Attributes,
                                       corrMethod = corrMethod,
                                       IDwindow  = alpha_sfs$IDwindow)

          # functional diversity metrics
          fd_sfs <- fd_metrics_sfs(fd_metrics = fd_metrics,
                                   obs2optimize = obs2optimize,
                                   rast_val = rast_val[SelFeat_tmp],
                                   Kmeans_info = Kmeans_info, IDplot = IDplot,
                                   corrMethod = corrMethod)

          for (alpha in alpha_metrics){
            Assess[[alpha]] <- alpha_sfs$Assess[[alpha]]
            CorrVal[[alpha]] <- alpha_sfs$CorrVal[[alpha]]
          }
          for (beta in beta_metrics){
            Assess[[beta]] <- beta_sfs$Assess[[beta]]
            CorrVal[[beta]] <- beta_sfs$CorrVal[[beta]]

          }
          for (fd in fd_metrics){
            Assess[[fd]] <- fd_sfs$Assess[[fd]]
            CorrVal[[fd]] <- fd_sfs$CorrVal[[fd]]
          }
          return(list('crit2Opt' = CorrVal, 'AssessedVal' = Assess))
        }
      }
      subSFS <- subfeatures_SFS()
      p(message = sprintf('Perform feature selection %g', nbvars2select))
      CorrSFS[[nbvars2select]] <- list()
      for (ind in FullListIndices)
        CorrSFS[[nbvars2select]][[ind]] <- unlist(lapply(lapply(subSFS,
                                                                '[[','crit2Opt'),
                                                         '[[',ind))
      CorrSFS[[nbvars2select]] <- data.frame(CorrSFS[[nbvars2select]])
      rownames(CorrSFS[[nbvars2select]]) <- AllVars

      # define variable maximizing a criterion
      SelVar <- maximize_sfs_criterion(CorrSFS, obs_criterion, nbvars2select)
      AssessSFS[[nbvars2select]] <- list()
      for (ind in FullListIndices)
        AssessSFS[[nbvars2select]][[ind]] <- unlist((lapply(lapply(subSFS,'[[',
                                                                   'AssessedVal'),
                                                            '[[',ind))[SelVar[1]])
      # which criterion to maximize with SFS?
      WhichVar <- rownames(CorrSFS[[nbvars2select]])[SelVar[1]]
      # add selected component to selected vars
      SelectedVars <- c(SelectedVars,WhichVar)
      # delete selected component from AllVars
      AllVars <- AllVars[-which(AllVars==WhichVar)]
      for (idx in FullListIndices)
        EvolCorr[[idx]] <- c(EvolCorr[[idx]], CorrSFS[[nbvars2select]][[idx]][SelVar[1]])

      # EvolCorr$richness <- c(EvolCorr$richness, CorrSFS[[nbvars2select]]$richness[SelVar[1]])
      # EvolCorr$shannon <- c(EvolCorr$shannon, CorrSFS[[nbvars2select]]$shannon[SelVar[1]])
      # EvolCorr$simpson <- c(EvolCorr$simpson, CorrSFS[[nbvars2select]]$simpson[SelVar[1]])
      # EvolCorr$hill <- c(EvolCorr$hill, CorrSFS[[nbvars2select]]$hill[SelVar[1]])
      # EvolCorr$BC <- c(EvolCorr$BC, CorrSFS[[nbvars2select]]$BC[SelVar[1]])
      # EvolCorr$FRic <- c(EvolCorr$FRic, CorrSFS[[nbvars2select]]$FRic[SelVar[1]])
      # EvolCorr$FEve <- c(EvolCorr$FEve, CorrSFS[[nbvars2select]]$FEve[SelVar[1]])
      # EvolCorr$FDiv <- c(EvolCorr$FDiv, CorrSFS[[nbvars2select]]$FDiv[SelVar[1]])
    }
  })
  parallel::stopCluster(cl)
  plan(sequential)

  EvolCorr <- data.frame(EvolCorr)
  rownames(EvolCorr) <- SelectedVars
  AssessSFS_comp <- list()
  for (ind in FullListIndices) {
    if (!is.null(AssessSFS[[1]][[ind]])) {
      AssessSFS_comp[[ind]] <- data.frame(lapply(AssessSFS,'[[',ind))
      names(AssessSFS_comp[[ind]]) <- SelectedVars
      AssessSFS_comp[[ind]]$obs <- c(obs2optimize[[ind]])
    } else {
      AssessSFS_comp[[ind]] <- NA
    }
  }
  return(list('SFS_perf' = EvolCorr, 'AssessDiv' = AssessSFS_comp))
}
