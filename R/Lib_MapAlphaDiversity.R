# ==============================================================================
# DiversityMappR
# Lib_MapAlphaDiversity.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2018/07 Jean-Baptiste FERET
# ==============================================================================
# This Library produces maps of alpha diversity indicators (Shannon, Simpson,
# Fischer...) based on spectral species file
# ==============================================================================

#' maps alpha diversity indicators  based on prior selection of PCs
#'
#' @param Input.Image.File Path and name of the image to be processed
#' @param Output.Dir output directory
#' @param Spatial.Unit dimensions of the spatial unit
#' @param TypePCA Type of PCA (PCA, SPCA, NLPCA...)
#' @param nbclusters number of clusters defined in k-Means
#' @param MinSun minimum proportion of sunlit pixels required to consider plot
#' @param pcelim minimum contribution (in %) required for a spectral species
#' @param FullRes
#' @param LowRes
#' @param nbCPU
#' @param MaxRAM
#' @param Index.Alpha 'Shannon', 'Simpson, 'Fisher'
#'
#' @return
#' @export
Map.Alpha.Diversity = function(Input.Image.File,Output.Dir,Spatial.Unit,TypePCA='SPCA',nbclusters=50,MinSun= 0.25,pcelim=0.02,Index.Alpha = 'Shannon',FullRes=TRUE,LowRes = FALSE,MapSTD = FALSE,nbCPU=FALSE,MaxRAM=FALSE){
  Output.Dir.SS         = Define.Output.SubDir(Output.Dir,Input.Image.File,TypePCA,'SpectralSpecies')
  Output.Dir.PCA        = Define.Output.SubDir(Output.Dir,Input.Image.File,TypePCA,'PCA')
  Spectral.Species.Path = paste(Output.Dir.SS,'SpectralSpecies',sep="")
  # 1- COMPUTE ALPHA DIVERSITY
  ALPHA         = Compute.Alpha.Diversity(Spectral.Species.Path,Spatial.Unit,nbclusters,MinSun,pcelim,nbCPU=nbCPU,MaxRAM=MaxRAM,Index.Alpha = Index.Alpha)
  # 2- SAVE ALPHA DIVERSITY MAPS
  print('Write alpha diversity maps')
  # which spectral indices will be computed
  Shannon = Simpson = Fisher = FALSE
  if (length((grep('Shannon',Index.Alpha)))>0)    Shannon = TRUE
  if (length((grep('Simpson',Index.Alpha)))>0)    Simpson = TRUE
  if (length((grep('Fisher',Index.Alpha)))>0)     Fisher = TRUE

  Output.Dir.Alpha  = Define.Output.SubDir(Output.Dir,Input.Image.File,TypePCA,'ALPHA')
  if (Shannon==TRUE){
    Index       = 'Shannon'
    Alpha.Path  = paste(Output.Dir.Alpha,Index,'_',Spatial.Unit,sep='')
    Write.Image.Alpha(ALPHA$Shannon,ALPHA$HDR,Alpha.Path,Spatial.Unit,Index,FullRes = FullRes, LowRes = LowRes)
    if (MapSTD==TRUE){
      Index       = 'Shannon_SD'
      Alpha.Path  = paste(Output.Dir.Alpha,Index,'_',Spatial.Unit,sep='')
      Write.Image.Alpha(ALPHA$Shannon.SD,ALPHA$HDR,Alpha.Path,Spatial.Unit,Index,FullRes = FullRes, LowRes = LowRes)
    }
  }

  if (Fisher==TRUE){
    Index       = 'Fisher'
    Alpha.Path  = paste(Output.Dir.Alpha,Index,'_',Spatial.Unit,sep='')
    Write.Image.Alpha(ALPHA$Fisher,ALPHA$HDR,Alpha.Path,Spatial.Unit,Index,FullRes = FullRes, LowRes = LowRes)
    if (MapSTD==TRUE){
      Index       = 'Fisher_SD'
      Alpha.Path  = paste(Output.Dir.Alpha,Index,'_',Spatial.Unit,sep='')
      Write.Image.Alpha(ALPHA$Fisher.SD,ALPHA$HDR,Alpha.Path,Spatial.Unit,Index,FullRes = FullRes, LowRes = LowRes)
    }
  }

  if (Simpson==TRUE){
    Index       = 'Simpson'
    Alpha.Path  = paste(Output.Dir.Alpha,Index,'_',Spatial.Unit,sep='')
    Write.Image.Alpha(ALPHA$Simpson,ALPHA$HDR,Alpha.Path,Spatial.Unit,Index,FullRes = FullRes, LowRes = LowRes)
    if (MapSTD==TRUE){
      Index       = 'Simpson_SD'
      Alpha.Path  = paste(Output.Dir.Alpha,Index,'_',Spatial.Unit,sep='')
      Write.Image.Alpha(ALPHA$Simpson.SD,ALPHA$HDR,Alpha.Path,Spatial.Unit,Index,FullRes = FullRes, LowRes = LowRes)
    }
  }
  return()
}

#' maps alpha diversity indicators  based on prior selection of PCs
#'
#' @param ImNames Path and name of the images to be processed
#' @param Output.Dir output directory
#' @param Spatial.Unit dimensions of the spatial unit
#' @param TypePCA Type of PCA (PCA, SPCA, NLPCA...)
#' @param nbclusters number of clusters defined in k-Means
#' @param MinSun minimum proportion of sunlit pixels required to consider plot
#' @param pcelim minimum contribution (in %) required for a spectral species
#' @param FullRes
#' @param LowRes
#' @param nbCPU
#' @param MaxRAM
#' @param Index.Alpha 'Shannon', 'Simpson, 'Fisher'
#'
#' @return
Map.Alpha.Diversity.TestnbCluster = function(ImNames,Output.Dir,Spatial.Unit,TypePCA='SPCA',nbclusters=50,MinSun= 0.25,pcelim=0.02,Index.Alpha = 'Shannon',FullRes=TRUE,LowRes = FALSE,nbCPU=FALSE,MaxRAM=FALSE){
  Output.Dir.SS         = Define.Output.SubDir(Output.Dir,ImNames$Input.Image,TypePCA,paste('SpectralSpecies_',nbclusters,sep=''))
  Output.Dir.PCA        = Define.Output.SubDir(Output.Dir,ImNames$Input.Image,TypePCA,'PCA')
  Spectral.Species.Path = paste(Output.Dir.SS,'SpectralSpecies',sep="")
  # 1- COMPUTE ALPHA DIVERSITY
  ALPHA         = Compute.Alpha.Diversity(Spectral.Species.Path,Spatial.Unit,nbclusters,MinSun,pcelim,nbCPU=nbCPU,MaxRAM=MaxRAM,Index.Alpha = Index.Alpha)
  # 2- SAVE ALPHA DIVERSITY MAPS
  print('Write alpha diversity maps')
  # which spectral indices will be computed
  Shannon = Simpson = Fisher = FALSE
  if (length((grep('Shannon',Index.Alpha)))>0)    Shannon = TRUE
  if (length((grep('Simpson',Index.Alpha)))>0)    Simpson = TRUE
  if (length((grep('Fisher',Index.Alpha)))>0)     Fisher = TRUE

  Output.Dir.Alpha  = Define.Output.SubDir(Output.Dir,ImNames$Input.Image,TypePCA,paste('ALPHA_',nbclusters,sep=''))
  if (Shannon==TRUE){
    Index       = 'Shannon'
    Alpha.Path  = paste(Output.Dir.Alpha,Index,'_',Spatial.Unit,sep='')
    Write.Image.Alpha(ALPHA$Shannon,ALPHA$HDR,Alpha.Path,Spatial.Unit,Index,FullRes = FullRes, LowRes = LowRes)
    Index       = 'Shannon_SD'
    Alpha.Path  = paste(Output.Dir.Alpha,Index,'_',Spatial.Unit,sep='')
    Write.Image.Alpha(ALPHA$Shannon.SD,ALPHA$HDR,Alpha.Path,Spatial.Unit,Index,FullRes = FullRes, LowRes = LowRes)
  }

  if (Fisher==TRUE){
    Index       = 'Fisher'
    Alpha.Path  = paste(Output.Dir.Alpha,Index,'_',Spatial.Unit,sep='')
    Write.Image.Alpha(ALPHA$Fisher,ALPHA$HDR,Alpha.Path,Spatial.Unit,Index,FullRes = FullRes, LowRes = LowRes)
    Index       = 'Fisher_SD'
    Alpha.Path  = paste(Output.Dir.Alpha,Index,'_',Spatial.Unit,sep='')
    Write.Image.Alpha(ALPHA$Fisher.SD,ALPHA$HDR,Alpha.Path,Spatial.Unit,Index,FullRes = FullRes, LowRes = LowRes)
  }

  if (Simpson==TRUE){
    Index       = 'Simpson'
    Alpha.Path  = paste(Output.Dir.Alpha,Index,'_',Spatial.Unit,sep='')
    Write.Image.Alpha(ALPHA$Simpson,ALPHA$HDR,Alpha.Path,Spatial.Unit,Index,FullRes = FullRes, LowRes = LowRes)
    Index       = 'Simpson_SD'
    Alpha.Path  = paste(Output.Dir.Alpha,Index,'_',Spatial.Unit,sep='')
    Write.Image.Alpha(ALPHA$Simpson.SD,ALPHA$HDR,Alpha.Path,Spatial.Unit,Index,FullRes = FullRes, LowRes = LowRes)
  }
  return()
}

#' Map alpha diversity metrics based on spectral species
#'
#' @param Spectral.Species.Path path for spectral species file to be written
#' @param Spatial.Unit size of spatial units (in pixels) to compute diversity
#' @param nbclusters number of clusters defined in k-Means
#' @param pcelim
#' @param nbCPU
#' @param MaxRAM
#' @param Index.Alpha
#' @param MinSun minimum proportion of sunlit pixels required to consider plot
#'
#' @return list of mean and SD of alpha diversity metrics
Compute.Alpha.Diversity = function(Spectral.Species.Path,Spatial.Unit,nbclusters,MinSun,pcelim,nbCPU=FALSE,MaxRAM=FALSE,Index.Alpha = 'Shannon'){
  ##      read SpectralSpecies file and write distribution per spatial unit   ##
  SS.HDR          = Get.HDR.Name(Spectral.Species.Path)
  HDR.SS          = read.ENVI.header(SS.HDR)
  if (MaxRAM == FALSE){
    MaxRAM = 0.25
  }
  nbPieces.Min    = Split.Image(HDR.SS,MaxRAM)
  if (nbCPU==FALSE){
    nbCPU = detectCores()
  }
  if (nbPieces.Min<nbCPU){
    nbPieces.Min  = nbCPU
  }
  SeqRead.SS      = Where.To.Read.Kernel(HDR.SS,nbPieces.Min,Spatial.Unit)

  ##          prepare SS distribution map and corresponding sunlit map        ##
  # prepare SS distribution map
  SSD.Path        = paste(Spectral.Species.Path,'_Distribution',sep='')
  HDR.SSD         = HDR.SS
  # define number of bands
  HDR.SSD$bands   = HDR.SS$bands*nbclusters
  # define image size
  HDR.SSD$samples = floor(HDR.SS$samples/Spatial.Unit)
  HDR.SSD$lines   = floor(HDR.SS$lines/Spatial.Unit)
  # change resolution
  HDR.SSD         = Change.Resolution.HDR(HDR.SSD,Spatial.Unit)
  HDR.SSD$`band names`= NULL
  # create SSD file
  fidSSD          = file(description = SSD.Path, open = "wb", blocking = TRUE,
                         encoding = getOption("encoding"), raw = FALSE)
  close(fidSSD)
  headerFpath     = paste(SSD.Path,'.hdr',sep='')
  write.ENVI.header(HDR.SSD, headerFpath)
  SeqWrite.SSD    = Where.To.Write.Kernel(HDR.SS,HDR.SSD,nbPieces.Min,Spatial.Unit)

  # prepare proportion of sunlit pixels from each spatial unit
  Sunlit.Path     = paste(SSD.Path,'_Sunlit',sep='')
  HDR.Sunlit      = HDR.SSD
  # define number of bands
  HDR.Sunlit$bands= 1
  # define number of bands
  HDR.Sunlit$`data type` = 4
  # create SSD Sunlit mask
  fidSunlit       = file(description = Sunlit.Path, open = "wb", blocking = TRUE,
                         encoding = getOption("encoding"), raw = FALSE)
  close(fidSunlit)
  headerFpath     = paste(Sunlit.Path,'.hdr',sep='')
  write.ENVI.header(HDR.Sunlit, headerFpath)
  SeqWrite.Sunlit = Where.To.Write.Kernel(HDR.SS,HDR.Sunlit,nbPieces.Min,Spatial.Unit)

  # for each piece of image
  ReadWrite = list()
  for (i in 1:nbPieces.Min){
    ReadWrite[[i]] = list()
    ReadWrite[[i]]$RW.SS  = ReadWrite[[i]]$RW.SSD = ReadWrite[[i]]$RW.Sunlit = list()
    ReadWrite[[i]]$RW.SS$Byte.Start     = SeqRead.SS$ReadByte.Start[i]
    ReadWrite[[i]]$RW.SS$nbLines        = SeqRead.SS$Lines.Per.Chunk[i]
    ReadWrite[[i]]$RW.SS$lenBin         = SeqRead.SS$ReadByte.End[i]-SeqRead.SS$ReadByte.Start[i]+1

    ReadWrite[[i]]$RW.SSD$Byte.Start    = SeqWrite.SSD$ReadByte.Start[i]
    ReadWrite[[i]]$RW.SSD$nbLines       = SeqWrite.SSD$Lines.Per.Chunk[i]
    ReadWrite[[i]]$RW.SSD$lenBin        = SeqWrite.SSD$ReadByte.End[i]-SeqWrite.SSD$ReadByte.Start[i]+1

    ReadWrite[[i]]$RW.Sunlit$Byte.Start = SeqWrite.Sunlit$ReadByte.Start[i]
    ReadWrite[[i]]$RW.Sunlit$nbLines    = SeqWrite.Sunlit$Lines.Per.Chunk[i]
    ReadWrite[[i]]$RW.Sunlit$lenBin     = SeqWrite.Sunlit$ReadByte.End[i]-SeqWrite.Sunlit$ReadByte.Start[i]+1
  }
  ImgFormat     = '3D'
  SSD.Format    = ENVI.Type2Bytes(HDR.SSD)
  SS.Format     = ENVI.Type2Bytes(HDR.SS)
  Sunlit.Format = ENVI.Type2Bytes(HDR.Sunlit)

  # multiprocess of spectral species distribution and alpha diversity metrics
  plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
  Schedule.Per.Thread = ceiling(nbPieces.Min/nbCPU)
  ALPHA               = future_lapply(ReadWrite, FUN = Convert.PCA.to.SSD, Spectral.Species.Path = Spectral.Species.Path,
                                      HDR.SS = HDR.SS, HDR.SSD = HDR.SSD, SS.Format = SS.Format, SSD.Format = SSD.Format,
                                      ImgFormat = ImgFormat, Spatial.Unit = Spatial.Unit, nbclusters = nbclusters, MinSun = MinSun,
                                      pcelim = pcelim, Index.Alpha = Index.Alpha, SSD.Path = SSD.Path, Sunlit.Path= Sunlit.Path,
                                      Sunlit.Format= Sunlit.Format, future.scheduling = Schedule.Per.Thread)
  plan(sequential)
  # create ful map from chunks
  Shannon.Mean.Chunk = Fisher.Mean.Chunk = Simpson.Mean.Chunk = list()
  Shannon.SD.Chunk = Fisher.SD.Chunk = Simpson.SD.Chunk = list()
  for (i in 1:length(ALPHA)){
    Shannon.Mean.Chunk[[i]] = ALPHA[[i]]$Shannon
    Fisher.Mean.Chunk[[i]]  = ALPHA[[i]]$Fisher
    Simpson.Mean.Chunk[[i]] = ALPHA[[i]]$Simpson
    Shannon.SD.Chunk[[i]]   = ALPHA[[i]]$Shannon.SD
    Fisher.SD.Chunk[[i]]    = ALPHA[[i]]$Fisher.SD
    Simpson.SD.Chunk[[i]]   = ALPHA[[i]]$Simpson.SD
  }
  Shannon.Mean  = do.call(rbind,Shannon.Mean.Chunk)
  Fisher.Mean   = do.call(rbind,Fisher.Mean.Chunk)
  Simpson.Mean  = do.call(rbind,Simpson.Mean.Chunk)
  Shannon.SD    = do.call(rbind,Shannon.SD.Chunk)
  Fisher.SD     = do.call(rbind,Fisher.SD.Chunk)
  Simpson.SD    = do.call(rbind,Simpson.SD.Chunk)
  my_list <- list("Shannon"=Shannon.Mean,"Fisher"=Fisher.Mean,"Simpson"=Simpson.Mean,
                  "Shannon.SD"=Shannon.SD,"Fisher.SD"=Fisher.SD,"Simpson.SD"=Simpson.SD,"HDR"=HDR.SSD)
  return(my_list)
}

#' Convert PCA into SSD based on previous clustering
#'
#' @param ReadWrite
#' @param Spectral.Species.Path
#' @param HDR.SS
#' @param HDR.SSD
#' @param SS.Format
#' @param SSD.Format
#' @param ImgFormat
#' @param Spatial.Unit
#' @param nbclusters
#' @param MinSun
#' @param pcelim
#' @param Index.Alpha
#' @param SSD.Path
#' @param Sunlit.Path
#' @param Sunlit.Format
#'
#' @param
#' @param
#' @return
Convert.PCA.to.SSD = function(ReadWrite, Spectral.Species.Path, HDR.SS, HDR.SSD,
                              SS.Format, SSD.Format, ImgFormat, Spatial.Unit, nbclusters,
                              MinSun, pcelim, Index.Alpha, SSD.Path, Sunlit.Path, Sunlit.Format){
  SS.Chunk      = Read.Image.Subset(Spectral.Species.Path,HDR.SS,
                                    ReadWrite$RW.SS$Byte.Start,ReadWrite$RW.SS$lenBin,
                                    ReadWrite$RW.SS$nbLines,SS.Format,ImgFormat)
  SSD.Alpha     = Compute.SSD(SS.Chunk,Spatial.Unit,nbclusters,MinSun,pcelim,Index.Alpha=Index.Alpha)
  # write spectral Species ditribution file
  fidSSD        = file(description = SSD.Path, open = "r+b", blocking = TRUE,
                       encoding = getOption("encoding"), raw = FALSE)
  if (!ReadWrite$RW.SSD$Byte.Start==1){
    seek(fidSSD, where = ReadWrite$RW.SSD$Byte.Start-1, origin = "start", rw = "write")
  }
  SSD.Chunk     = aperm(array(SSD.Alpha$SSD,c(ReadWrite$RW.SSD$nbLines,HDR.SSD$samples,HDR.SSD$bands)),c(2,3,1))
  writeBin(c(SSD.Chunk), fidSSD, size = SSD.Format$Bytes,endian = .Platform$endian, useBytes = FALSE)
  close(fidSSD)

  # write PCsunlit pixels corresponding to SSD file
  fidSunlit     = file(description = Sunlit.Path, open = "r+b", blocking = TRUE,
                       encoding = getOption("encoding"), raw = FALSE)
  if (!ReadWrite$RW.Sunlit$Byte.Start==1){
    seek(fidSunlit, where = ReadWrite$RW.Sunlit$Byte.Start-1, origin = "start", rw = "write")
  }
  Sunlit.Chunk  = t(SSD.Alpha$PCsun)
  writeBin(c(Sunlit.Chunk), fidSunlit, size = Sunlit.Format$Bytes,endian = .Platform$endian, useBytes = FALSE)
  close(fidSunlit)

  rm(SSD.Chunk)
  rm(Sunlit.Chunk)
  gc()
  Shannon.Mean.Chunk  = apply(SSD.Alpha$Shannon,1:2,mean)
  Fisher.Mean.Chunk   = apply(SSD.Alpha$Fisher,1:2,mean)
  Simpson.Mean.Chunk  = apply(SSD.Alpha$Simpson,1:2,mean)
  Shannon.SD.Chunk    = apply(SSD.Alpha$Shannon,1:2,sd)
  Fisher.SD.Chunk     = apply(SSD.Alpha$Fisher,1:2,sd)
  Simpson.SD.Chunk    = apply(SSD.Alpha$Simpson,1:2,sd)
  rm(SSD.Alpha)
  gc()
  my_list <- list("Shannon"=Shannon.Mean.Chunk,"Fisher"=Fisher.Mean.Chunk,"Simpson"=Simpson.Mean.Chunk,
                  "Shannon.SD"=Shannon.SD.Chunk,"Fisher.SD"=Fisher.SD.Chunk,"Simpson.SD"=Simpson.SD.Chunk)
  return(my_list)
}

#' compute spectral species distribution from original spectral species map
#'
#' @param Image.Chunk 3D image chunk of spectral species
#' @param Spatial.Unit size of spatial units (in pixels) to compute diversity
#' @param nbclusters number of clusters defined in k-Means
#' @param MinSun minimum proportion of sunlit pixels required to consider plot
#' @param Index.Alpha
#' @param pcelim minimum proportion for a spectral species to be included
#'
#' @return list of alpha diversity metrics for each iteration
Compute.SSD = function(Image.Chunk,Spatial.Unit,nbclusters,MinSun,pcelim,Index.Alpha='Shannon'){

  nbi         = floor(size(Image.Chunk)[1]/Spatial.Unit)
  nbj         = floor(size(Image.Chunk)[2]/Spatial.Unit)
  nbIter      = size(Image.Chunk)[3]
  SSDMap      = array(NA,c(nbi,nbj,nbIter*nbclusters))
  shannonIter = FisherAlpha = SimpsonAlpha = array(NA,dim=c(nbi,nbj,nbIter))
  PCsun       = matrix(NA,nrow=nbi,ncol=nbj)

  # which spectral indices will be computed
  Shannon = Simpson = Fisher = FALSE
  if (length((grep('Shannon',Index.Alpha)))>0)    Shannon = TRUE
  if (length((grep('Simpson',Index.Alpha)))>0)    Simpson = TRUE
  if (length((grep('Fisher',Index.Alpha)))>0)     Fisher = TRUE

  # for each kernel in the line
  for (ii in 1:nbi){
    # for each kernel in the column
    for (jj in 1:nbj){
      li      = ((ii-1)*Spatial.Unit)+1
      ui      = ii*Spatial.Unit
      lj      = ((jj-1)*Spatial.Unit)+1
      uj      = jj*Spatial.Unit
      # put all iterations in a 2D matrix shape
      ijit    = t(matrix(Image.Chunk[li:ui,lj:uj,],ncol=nbIter))
      # keep non zero values
      ijit          = matrix(ijit[,which(!ijit[1,]==0)],nrow = nbIter)
      nb.Pix.Sunlit = size(ijit)[2]
      PCsun[ii,jj]= nb.Pix.Sunlit/Spatial.Unit**2
      if (PCsun[ii,jj]>MinSun){
        # for each iteration
        for (it in 1:nbIter){
          lbk         = (it-1)*nbclusters
          SSD         = as.vector(table(ijit[it,]))
          ClusterID   = sort(unique(ijit[it,]))
          if (pcelim>0){
            KeepSS      = which(SSD>=pcelim*nb.Pix.Sunlit)
            ClusterID   = ClusterID[KeepSS]
            SSD         = SSD[KeepSS]
          }
          SSDMap[ii,jj,(lbk+ClusterID)] = SSD
          if (Shannon==TRUE){
            shannonIter[ii,jj,it]  = Get.Shannon(SSD)
          }
          if (Simpson==TRUE){
            SimpsonAlpha[ii,jj,it] = Get.Simpson(SSD)
          }
          if (Fisher==TRUE){
            if (length(SSD)>2){
              FisherAlpha[ii,jj,it]  = fisher.alpha(SSD)
            } else {
              FisherAlpha[ii,jj,it]  = 0
            }
          }
        }
      } else {
        shannonIter[ii,jj,]  = NA
        FisherAlpha[ii,jj,]  = NA
        SimpsonAlpha[ii,jj,] = NA
      }
    }
  }
  my_list <- list("Shannon"=shannonIter,"Fisher"=FisherAlpha,"Simpson"=SimpsonAlpha,"SSD"=SSDMap,"PCsun"=PCsun)
  return(my_list)
}

#' computes shannon index from a distribution
#'
#' @param Distrib Distribution
#'
#' @return Shannon index correspnding to the distribution
Get.Shannon=function(Distrib){

  Distrib     = Distrib/sum(Distrib,na.rm=TRUE)
  Distrib     = Distrib[which(!Distrib==0)]
  shannon     = -1*sum(Distrib*log(Distrib),na.rm=TRUE)
  return(shannon)
}

#' computes Simpson index from a distribution
#'
#' @param Distrib Distribution
#'
#' @return Simpson index correspnding to the distribution
Get.Simpson=function(Distrib){

  Distrib   = Distrib/sum(Distrib,na.rm=TRUE)
  Simpson   = 1-sum(Distrib*Distrib,na.rm=TRUE)
  return(Simpson)
}

#' Writes image of alpha diversity indicator (1 band) and smoothed alpha diversity
#'
#' @param Image 2D matrix of image to be written
#' @param HDR.SSD hdr template derived from SSD to modify
#' @param ImagePath path of image file to be written
#' @param Spatial.Unit spatial units use dto compute diversiy (in pixels)
#' @param Index name of the index (eg. Shannon)
#' @param FullRes should full resolution image be written (original pixel size)
#' @param LowRes should low resolution image be written (one value per spatial unit)
#'
#' @return
Write.Image.Alpha=function(Image,HDR.SSD,ImagePath,Spatial.Unit,Index,FullRes = TRUE, LowRes = FALSE){

  # Write image with resolution corresponding to Spatial.Unit
  HDR.Alpha             = HDR.SSD
  HDR.Alpha$bands       = 1
  HDR.Alpha$`data type` = 4
  HDR.Alpha$`band names`= Index
  Image.Format  = ENVI.Type2Bytes(HDR.Alpha)
  if (LowRes ==TRUE){
    headerFpath = paste(ImagePath,'.hdr',sep='')
    write.ENVI.header(HDR.Alpha, headerFpath)
    ImgWrite    = array(Image,c(HDR.Alpha$lines,HDR.Alpha$samples,1))
    ImgWrite    = aperm(ImgWrite,c(2,3,1))
    fidOUT      = file(description = ImagePath, open = "wb", blocking = TRUE,
                       encoding = getOption("encoding"), raw = FALSE)
    writeBin(c(ImgWrite), fidOUT, size = Image.Format$Bytes,endian = .Platform$endian, useBytes = FALSE)
    close(fidOUT)
  }
  if (FullRes ==TRUE){
    # Write image with Full native resolution
    HDR.Full          = HDR.Alpha
    HDR.Full$samples  = HDR.Alpha$samples*Spatial.Unit
    HDR.Full$lines    = HDR.Alpha$lines*Spatial.Unit
    HDR.Full          = Revert.Resolution.HDR(HDR.Full,Spatial.Unit)
    ImagePath.FullRes = paste(ImagePath,'_Fullres',sep='')
    headerFpath       = paste(ImagePath.FullRes,'.hdr',sep='')
    write.ENVI.header(HDR.Full, headerFpath)
    Image.Format  = ENVI.Type2Bytes(HDR.Full)
    Image.FullRes = matrix(NA,ncol=HDR.Full$samples,nrow=HDR.Full$lines)
    for (i in 1:HDR.SSD$lines){
      for (j in 1:HDR.SSD$samples){
        Image.FullRes[((i-1)*Spatial.Unit+1):(i*Spatial.Unit),((j-1)*Spatial.Unit+1):(j*Spatial.Unit)]=Image[i,j]
      }
    }
    ImgWrite    = array(Image.FullRes,c(HDR.Full$lines,HDR.Full$samples,1))
    ImgWrite    = aperm(ImgWrite,c(2,3,1))
    fidOUT      = file(description = ImagePath.FullRes, open = "wb", blocking = TRUE,
                       encoding = getOption("encoding"), raw = FALSE)
    writeBin(c(ImgWrite), fidOUT, size = Image.Format$Bytes,endian = .Platform$endian, useBytes = FALSE)
    close(fidOUT)
    # zip resulting file
    ZipFile(ImagePath.FullRes)
  }
  # write smoothed map
  SizeFilt  = 1
  nbi       = dim(Image)[1]
  nbj       = dim(Image)[2]
  if (min(c(nbi,nbj))>(2*SizeFilt+1)){
    Image.Smooth      = Mean.Filter(Image,nbi,nbj,SizeFilt)
    Image.Smooth[which(is.na(Image))] = NA
    ImagePath.Smooth  = paste(ImagePath,'_MeanFilter',sep='')
    headerFpath       = paste(ImagePath.Smooth,'.hdr',sep='')
    Image.Format      = ENVI.Type2Bytes(HDR.Alpha)
    if (LowRes ==TRUE){
      write.ENVI.header(HDR.Alpha, headerFpath)
      ImgWrite          = array(Image.Smooth,c(HDR.Alpha$lines,HDR.Alpha$samples,1))
      ImgWrite          = aperm(ImgWrite,c(2,3,1))
      fidOUT            = file(description = ImagePath.Smooth, open = "wb", blocking = TRUE,
                               encoding = getOption("encoding"), raw = FALSE)
      writeBin(c(ImgWrite), fidOUT, size = Image.Format$Bytes,endian = .Platform$endian, useBytes = FALSE)
      close(fidOUT)
    }
    if (FullRes ==TRUE){
      # Write image with Full native resolution
      ImagePath.FullRes = paste(ImagePath.Smooth,'_Fullres',sep='')
      headerFpath       = paste(ImagePath.FullRes,'.hdr',sep='')
      write.ENVI.header(HDR.Full, headerFpath)
      Image.Format  = ENVI.Type2Bytes(HDR.Full)
      Image.FullRes = matrix(NA,ncol=HDR.Full$samples,nrow=HDR.Full$lines)
      for (i in 1:HDR.SSD$lines){
        for (j in 1:HDR.SSD$samples){
          Image.FullRes[((i-1)*Spatial.Unit+1):(i*Spatial.Unit),((j-1)*Spatial.Unit+1):(j*Spatial.Unit)]=Image.Smooth[i,j]
        }
      }
      ImgWrite    = array(Image.FullRes,c(HDR.Full$lines,HDR.Full$samples,1))
      ImgWrite    = aperm(ImgWrite,c(2,3,1))
      fidOUT      = file(description = ImagePath.FullRes, open = "wb", blocking = TRUE,
                         encoding = getOption("encoding"), raw = FALSE)
      writeBin(c(ImgWrite), fidOUT, size = Image.Format$Bytes,endian = .Platform$endian, useBytes = FALSE)
      close(fidOUT)
      # zip resulting file
      ZipFile(ImagePath.FullRes)
    }
  }
  return("")
}
