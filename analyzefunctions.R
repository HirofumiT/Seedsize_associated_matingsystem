
initializing = function(){
  setwd('source')
  library("Matrix")
  library("lme4")
  library("stringr")
  library('conflicted')
  library("lmerTest") 
  library("MuMIn")
  #library("maptools")
  
  glmm_chart <<- data.frame(array(rep(0,7),dim = c(0,7)))
  aicscore_chart <<- data.frame(array(rep(0,3),dim = c(0,4)))
  glmm.pcn_chart <<- data.frame(array(rep(0,7),dim = c(0,7)))
  a_bra_tmp <<- read.csv('all/bra_analyzingset.csv')
  a_ast_tmp <<- read.csv('all/ast_analyzingset.csv')
  a_sol_tmp <<- read.csv('all/sol_analyzingset.csv')
  a_all_tmp <<- rbind(rbind(a_bra_tmp,a_ast_tmp),a_sol_tmp)
  
  h_bra_tmp <<- read.csv('herbonly/bra_analyzingset.csv')
  h_ast_tmp <<- read.csv('herbonly/ast_analyzingset.csv')
  h_sol_tmp <<- read.csv('herbonly/sol_analyzingset.csv')
  
  
  a_sol_foml <<- Y ~ XS  + log(srad)  + TP + wind  + (1|Z)
  a_ast_foml <<- Y ~ XS  + rf + tmin + varpr  + (1|Z)
  a_bra_foml <<- Y ~ XS  + rf + varpr + (1|Z)
  h_sol_foml <<- Y ~ XS  + log(XC) + tmin + wind  + (1|Z)
  h_ast_foml <<- Y ~ XS  + rf + tmin + varpr  + (1|Z)
  h_bra_foml <<- Y ~ XS  +log(srad) +  (1|Z)
  datanames <<- c('sol','bra','ast')
  setwd('/Users/hiro.t/Documents/GitHub/Seedsize_associated_matingsystem/source')
  
  
}

savesummarys = function(){
  names(glmm_chart) <<- c('famname','estimate-effect','Std.Error','df','tvalue','Pr','forml')
  write.csv(glmm_chart,file = 'analyzesummary/output/glmm-summary.csv',row.names = F)
  
  names(aicscore_chart) <<- c('famname','AICbestmodel','AICscore')
  write.csv(aicscore_chart,file = 'analyzesummary/output/AIC-PCn.csv',row.names = F)
  
  names(glmm.pcn_chart) <<- c('famname','estimate-effect','Std. Error','df','t value','Pr','forml')
  write.csv(glmm.pcn_chart,file = 'analyzesummary/output/PCnAICglmm-summary.csv',row.names = F)
  
}

d1tod4 = function(data,datanum){
 
  data = data[is.na(data$temp) != TRUE |is.na(data$cx) != TRUE,]
  data = data[is.na(data$rf) != TRUE,]
  
  #data$rf = log(data$rf)
  data = data[is.infinite(data$rf) != TRUE,]
  
  data$SCSI = as.character(data$SCSI)
  data$SCSI = gsub(pattern = ("2|3|5"),replacement = "SI",data$SCSI)
  data$SCSI = gsub(pattern = ("4"),replacement = "SC",data$SCSI)
  data$SCSI[is.na(data$SCSI) == TRUE]  = c("ND")  
  
  
  SCSIlist = c("SC","SI","ND")
  data_2 <- data[data$SCSI %in% SCSIlist,]
  fixGenus = list(rep,1:length(data_2$Genus),character(0))
  for(a in 1:length(data_2$Genus)){
   if(data_2$Genus[a] =="Solanum"){
      if(is.na(data_2$class5[a]) == FALSE & data_2$class5[a] != "n"){
        fixGenus[a] =c(as.character(paste("Solanum_",as.character(data_2$class5[a]))))
        print(a)
        print(data_2$Genus[a])}else{fixGenus[a] = c(as.character(data_2$Genus[a]))}
        }else{
      fixGenus[a] = c(as.character(data_2$Genus[a]))
        }
     }
  fixGenus=t(as.data.frame(fixGenus))
  rownames(fixGenus) = NULL
  attach(data_2)
  data_3=data.frame(Genus,Species_Epithet,SCSI,Weight_list,ID,fixGenus,cx,temp,tmax,tmin,varpr,srad,wind,rf,growth)
  names(data_3) <- c("Genus", "Species_Epither", "SCSI","Weight_list","ID","fixGenus","cx","temp","tmax","tmin","varpr","srad","wind","rf",'growth')
  rownames(data_3) = NULL 
  detach(data_2)
  
  numofdata <- sapply(1:length( unique(data_3$fixGenus) ), function(i){length(data_3[data_3$fixGenus == unique(data_3$fixGenus)[i], 1])})
  #print(numofdata)
  Genuslist=unique(data_3$fixGenus)
  fixGenuslist=list()
  for ( a in 1:length(Genuslist) ){
    b=Genuslist[a]
    if(is.na(b) == TRUE){
      print(paste("seeyou",b))
    }else{
      print(b)
      SClist = c(data_3$Weight_list[data_3$SCSI =="SC" & data_3$fixGenus == b])
      SIlist = c(data_3$Weight_list[data_3$SCSI =="SI" & data_3$fixGenus == b])
      #SClist = SClist[SClist != "NULL"]
      #SIlist = SClist[SIlist != "NULL"]
      SClength = length(SClist)
      SIlength = length(SIlist)
      if((SClength < 2)|(SIlength < 2)){
        #print(paste("seeyou",b))
      }else{
        if(b == "Solanum_ "){ ##
          print(paste("seeyou",b))
        }else{
          tempSClist = c(data_3$temp[data_3$SCSI =="SC" & data_3$fixGenus == b])
          tempSIlist = c(data_3$temp[data_3$SCSI =="SI" & data_3$fixGenus == b])
          tempSClength = length(tempSClist)
          tempSIlength =  length(tempSIlist)
          if((tempSClength < 2)|(tempSIlength < 2)){
            print(paste("tempseeyou",b))
          }else{
            fixGenuslist=c(fixGenuslist,as.character(b))
            #print(paste("welcome",b))
          }
        }
      }
    }
  }

  if(datanum == '2'){
    rejectlist = c("Solanum_","Solanum_ nouse")
    fixGenuslist =  fixGenuslist[-which(fixGenuslist %in% rejectlist)]
  }else{
    print("No reject")
  }
  
  data_4 = data_3[data_3$fixGenus %in% fixGenuslist,]
  data_4$fixGenus = droplevels(data_4$fixGenus)
  data_4$Weight_list = log(data_4$Weight_list)#############attach log here!!!!!!!!!!!!
  
  data_4 = data_4[is.na(data_4$Weight_list) == FALSE,]
  data_4 = data_4[is.na(data_4$rf) == FALSE,]
  return(data_4)
}########################atatch log in this function!!!!!!!!!!!!!!

d_4toAnaly = function(d4){
  
  d4 = d4[d4$SCSI != 'ND',]
  XS = d4$SCSI
  XS = lapply(XS, gsub, pattern="SC", replacement = 1)
  XS = lapply(XS, gsub, pattern="SI", replacement = 0)
  #XS = lapply(XS, gsub, pattern="ND", replacement = NA)
  
  XS = as.numeric(XS)
  XC = d4$cx
  Z = as.factor(d4$fixGenus)
  genuscount = unique(Z)
  tmp <- data.frame(d4$Weight_list,XS,XC,d4$temp,d4$rf,d4$tmax,d4$tmin,d4$varpr,d4$srad,d4$wind,d4$growth,Z)
  names(tmp) = c("Y","XS","XC","TP","rf","tmax","tmin","varpr","srad","wind","growth","Z")
  return(tmp)
}

doAIC = function(analy,fam){
  options(na.action = "na.fail")
  analy_deg = dredge(analy,rank = "AIC", fixed = NULL, m.lim = NULL,trace = TRUE)
  return(analy_deg)
}


OutputSummary_AICed = function(analy,fam){
  
  
  
  #write(paste("data:",length(tmp$Y),"\n"),file = str_c('output/',fam,'_a_AICed.txt',sep =""),append = T)
  write(capture.output(summary(analy)),file = str_c('analyzesummary/output/',fam,'_a_AICed.txt',sep =""),append = T)
  write(capture.output(ranef(analy)),file = str_c('analyzesummary/output/',fam,'_a_AICed.txt',sep =""),append = T)
}

fPCA =function(templ,fam){
  templ$srad = log(templ$srad)
  templ.pca = prcomp(templ[5:10],scale = TRUE)
  save(templ.pca, file = paste('analyzesummary/rdata/',fam,'.rda',sep = ''))
  tmp_2 = transform(templ,PC1  = templ.pca$x[,1],PC2  = templ.pca$x[,2],PC3 = templ.pca$x[,3],PC4 = templ.pca$x[,4],PC5 = templ.pca$x[,5],PC6 = templ.pca$x[,6])
  #Panaly = lmerTest::lmer(Y ~ XS + PC1 + PC2 + (1|Z),data=tmp_2 ,REML = FALSE)
  #print(summary(templ.pca))
  #biplot(templ.pca,scale = FALSE,main=fam)
  # write.csv(tmp_2,paste('fromAnalyzer/',fam,'_analyzingset_pca`.csv',sep=''),row.names = F)
  #write(capture.output(summary(Panaly)),file = str_c('output/',fam,'_a_PCA.txt',sep =""),append = T)
  #write(capture.output(ranef(Panaly)),file = str_c('output/',fam,'_a_PCA.txt',sep =""),append = T)
  write.csv(templ.pca$x,file = str_c('analyzesummary/output/pcasummary_',fam,'.csv',sep = ''),row.names = F)
  
  return(tmp_2)
}


savetxtglmm = function(dataset,forml,famname,analy){
  txtsavedir = str_c('analyzesummary/output/GLMM_',famname,'.txt',sep ="")
  print(forml)
  write(paste(forml,sep = ''),file = txtsavedir)
  write(capture.output(summary(analy)),file = txtsavedir, append = T)
  write(capture.output(ranef(analy)),file = txtsavedir,append = T)
}
doandsaveglmm = function(dataset,forml,famname){
  analy = lmerTest::lmer(forml,dataset, REML = FALSE)
  savetxtglmm(dataset,forml,famname,analy)
  save(analy,file = sprintf('analyzesummary/rdata/GLMM_%s.rda',famname))
  return(analy)
}

pcaAICchart = function(aic){
  best = as.data.frame(aic)[1,]
  print(best)
  print(best[2:6])
  bestforml = ''
  for(i in 1:5){
    if(colnames(best[i+1]) == 'df'){
      break
    }else{
      if(is.na(best[i+1]) == FALSE){
        #print(colnames(best[i+2]))
        if(bestforml != ''){
          bestforml = paste(bestforml,colnames(best[i+1]),sep = ' + ')
          #print(bestforml)
        }else{
          #print(bestforml)
          bestforml = colnames(best[i+1])
        }
        
      }else{}
    }
  }
  bestforml = sprintf('Y ~ %s + (1|Z)',bestforml)
  
  return(bestforml)
}

