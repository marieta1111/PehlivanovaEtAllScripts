#######################################################################################
# GAM  model function to run models over multiple networks and extract stats for logk #
# 									              #
# ms:       statements of models in the form that lm() accepts		              #
# dataset:  which dataset to use					  	      #
# varNames: row names of individual networks				              #
# s: 	    1=main variable of interest is linear 			              #
#	    2=main variable is modeled with splines			              #	
#######################################################################################

regDataGam <- function(ms,dataset,varNames,s) {
        # s is a splines variable, if 1, no splines on logk
        # if 2, splines on logk
        nn<-length(varNames)
        model.results<-lapply(ms, function(x) {
        foo<-summary(gam(as.formula(x), data=dataset, method="REML"))
        # takes stats for first covariate, specifically
        if (s==1) {
                return(c(foo$n,foo$p.table[2,]))
        } else if (s==2) {
                return(c(foo$n,foo$s.table[1,]))
        }
        })
        # convert results to data frame
        gamOut1<-data.frame(matrix(unlist(model.results), nrow=nn, byrow=T))
        if (s==1) {
                 names(gamOut1) <-c("n","coef","se","tval","p")
        } else if (s==2) {
                 names(gamOut1) <-c("n","edf","rfDF","F","p")
        }
        row.names(gamOut1)<-varNames
        # multiple comparison testing
        gamOut1[,6]<-p.adjust(gamOut1[,5], method = "fdr", n = nn)
        gamOut1[,7]<-ms
        names(gamOut1)[6:7]<-c("fdr.p","model")
        return(gamOut1)
}


