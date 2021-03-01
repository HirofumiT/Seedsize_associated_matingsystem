**# Seedsize_associated_matingsystem**

Detail
- Analysis script and original full data of the paper.
source/{Asteraceae|Brassicaceae|Solanaceae} .csv  

- Analysis script
mainscript.R

How to run code

1.  Download this repository
`$git clone https://github.com/HirofumiT/Seedsize_associated_matingsystem.git`
or
press “Code”, then choose “Download ZIP”
2. Open ‘mainscript.R’ , then chenge working dir
```r
setwd('YOUR_FILE_ADDRESS/Seedsize_associated_matingsystem')   
```
3. Loading analysis functions and initializing  
```r
source( "analyzefunctions.R" )
initializing()
```

“all/XX_analyzingset.csv” or “herbonly/XX _analyzingset.csv” is loaded as analysis dataset.
4. Run analysis (PCA,AIC, generalized linear mixed model(GLMM))
```r
for (type in c('a','h')){
....
  }
savesummarys()
```
5. Analysis summary is in “Seedsize_associated_matingsystem/source/analyzesummary/output”
	- “PCnAICglmm-summary.csv” : summary of GLMM with bestmodel based on AIC (Raw data of Tables 3,4 in paper) . “a” or “h” attached on “famname” column means “all growth-form” or “only herbaceous species”
6. To run Single regression analysis by a GLMM
```r
##codes for GLMM Y~XS
for (type in c('a')){
  ...      
  }
}
```

	
