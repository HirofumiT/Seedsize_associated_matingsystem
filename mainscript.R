##build 20191009

##DoMainGLMM_Analysis ##
## Run code below
setwd('YOUR_FILE_ADDRESS/Seedsize_associated_matingsystem')
source( "analyzefunctions.R" )
initializing()
for (type in c('a','h')){
  for (fam in datanames){
    
    famname = paste(type,fam,sep = '_')
    dataset = get(paste(famname,'tmp',sep = '_'))
    #forml = get(paste(famname,'foml',sep = '_'))
    name_analyfile = paste(famname,'analy',sep = '_')
    name_analyfilepca = paste('pca',famname,'tmp',sep = '_')
    print(famname)
    #analy_AICmodel = doandsaveglmm(dataset,forml,paste(famname,'pc1to3_AICmodel',sep = ''))
    #a_AIC.coef = as.data.frame(summary(analy_AICmodel)[["coefficients"]])
    #condition = rownames(a_AIC.coef) == 'XS'
    #score = data.frame(array(c(paste(famname,'pc1to3_AICmodel',sep=''),a_AIC.coef$Estimate[condition],a_AIC.coef$`Std. Error`[condition],a_AIC.coef$df[condition],a_AIC.coef$`t value`[condition],a_AIC.coef$`Pr(>|t|)`[condition],as.character(forml[3])),dim = c(1,7)))
    #glmm_chart =　rbind(glmm_chart,score) 
    
    pcadataset = fPCA(dataset,famname)
    write.csv(pcadataset,file = paste(name_analyfilepca,'.csv',sep = ""),row.names = F)
    pcaanaly = lmerTest::lmer(Y ~ XS + PC1 + PC2 + PC3  + (1|Z),data=pcadataset ,REML = FALSE)
    save(pcaanaly, file = sprintf('analyzesummary/rdata/glmm_pc1to3_%s.rda',famname))
    aic = doAIC(pcaanaly,famname)
    write.csv(aic,file=sprintf('analyzesummary/output/AIC_pc1to3_%s.csv',famname))
    
    aicbest = pcaAICchart(aic)
    print(paste(aicbest,famname))
    aicscore = data.frame(array(c(famname,aicbest,aic$AIC[1]),dim = c(1,3)))
    aicscore_chart = rbind(aicscore_chart,aicscore)
    
    aicedpca_analy = doandsaveglmm(pcadataset,aicbest,paste(famname,'_PCnAICmodel',sep = ''))
    pcnaic.coef = as.data.frame(summary(aicedpca_analy)[["coefficients"]])
    condition = rownames(pcnaic.coef) == 'XS'
    pcnscore = data.frame(array(c(paste(famname,'_PCnAICmodel',sep=''),pcnaic.coef$Estimate[condition],pcnaic.coef$`Std. Error`[condition],pcnaic.coef$df[condition],pcnaic.coef$`t value`[condition],pcnaic.coef$`Pr(>|t|)`[condition],as.character(aicbest)),dim = c(1,7)))
    
    glmm.pcn_chart =　rbind(glmm.pcn_chart,pcnscore)
    #rm(famname,dataset,forml,pcadataset,aic)
     
    
  }
}
savesummarys()
##End##
#summary can see in 'analyzesummary/output/PCnAICglmm-summary.csv'




##codes for GLMM Y~XS
for (type in c('a')){
  for (fam in datanames){
    param = 'XS'
    famname = paste(type,fam,sep = '_')
    dataset = get(paste(famname,'tmp',sep = '_'))
    forml = sprintf('Y ~ %s  + (1|Z)','XS')
    name_analyfile = paste(famname,param,'analy',sep = '_')
    print(famname)
    analy_AICmodel = doandsaveglmm(dataset,forml,paste(famname,param,'_param',sep = ''))
    a_AIC.coef = as.data.frame(summary(analy_AICmodel)[["coefficients"]])
    condition = rownames(a_AIC.coef) == param
    score = data.frame(array(c(paste(famname,'_',param,sep=''),a_AIC.coef$Estimate[condition],a_AIC.coef$`Std. Error`[condition],a_AIC.coef$df[condition],a_AIC.coef$`t value`[condition],a_AIC.coef$`Pr(>|t|)`[condition],as.character(forml)),dim = c(1,4)))
    #glmm_paramchart =　rbind(glmm_paramchart,score)
      
    
  }
}
