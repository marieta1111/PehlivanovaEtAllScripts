###########################################################################################################
# GAM  model function to run models for multiple networks and extract stats for logk and other covariates #
#           expanded version of regDataGam.R which only extracts stats for the first independent variable #
# 									                        	  #
# ms:       statements of models in the form that lm() accepts		                                  #
# dataset:  which dataset to use					                                  #
# varNames: row names of individual networks	                 				          #
# s: 	    1=main independent variable of interest is linear 		                                  #
#	    2=main independent variable is modeled with splines		                                  #
# nVar:     number of independent variables to output; default is 1		   			  #		
###########################################################################################################

regDataGamMultOut <- function(ms,dataset,varNames,s,nVar=1) {
        # s is a splines variable, if 1, no splines on logk
        # if 2, splines on logk
        nn<-length(varNames)
        model.results<-lapply(ms, function(x) {
        foo<-summary(gam(as.formula(x), data=dataset, method="REML"))
        # takes stats for first covariate, specifically
        if (s==1) {
        	if (nVar==1){
	        return(c(foo$n,foo$p.table[2,]))
	}
		else if (nVar!=1){
		return(c(foo$n,c(t(foo$p.table[2:(nVar+1),]))))
	}		
        } else if (s==2) {
                if (nVar==1){
		return(c(foo$n,foo$s.table[1,]))
	}
		else if (nVar!=1){
		return(c(foo$n,foo$s.table[1,],c(t(foo$p.table[2:nVar,]))))
	}
        }
        })
        # convert results to data frame
        gamOut1<-data.frame(matrix(unlist(model.results), nrow=nn, byrow=T))
        if (s==1) {
                 names(gamOut1) <- c("n",rep(c("coef","se","tval","p"),nVar))
        } else if (s==2) {
                 names(gamOut1) <-c("n","edf","rfDF","F","p",rep(c("coef","se","tval","p"),nVar-1))
        }
        row.names(gamOut1)<-varNames
        # multiple comparison testing
        gamOut1[,dim(gamOut1)[2]+1]<-p.adjust(gamOut1[,5], method = "fdr", n = nn)
        gamOut1[,dim(gamOut1)[2]+1]<-ms
        names(gamOut1)[(dim(gamOut1)[2]-1):dim(gamOut1)[2]]<-c("fdr.p","model")

        return(gamOut1)
}


