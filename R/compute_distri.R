compute_distri<-function(model, dfpred, time_var="year",
                         factor_name=NA, nreplicates=1000,...)
{


    ## Check that the model was run with parameters estimated through REML
    if(!model$method %in% c("REML", "ML"))
        print("Predictions computed through this function are using the unconditionnal=T option of the vcov function: this requires models to be fitted with method=REML or ML")

    ## Some information on the variables to properly index the resulting matrix
    ## Need to check if the factor_name isn't included in the model parameters -> this will
    ## mess up the results matrix (this could happen when the objective is to predict
    ## over a grid for example (though X-Y), as one needs to keep track of the cell id
    Xp<-predict.gam(model, newdata=dfpred, type="lpmatrix") #
    betas_mod<-coef(model)
    vcov_mod<-vcov(model, unconditionnal=T) #unconditionnal = T allows
                                            #parameter estimates to be
                                            #treated as uncertain at
                                            #their estimates, but this
                                            #requires the model to be
                                            #estimated through REML or
                                            #ML

    ## sample MVN values of beta
    betas_rep<-rmvn(n=nreplicates, betas_mod, vcov_mod)

    ## Product
    lp<-t(exp(Xp %*% t(betas_rep))) #

    ## If there is a structuring factor, store the results in a list
    if (!is.na(factor_name)) {

        ## A list to store the results
        list_distri<-list()

        ## Convert grouping factor if it is not
        if(!is.factor(dfpred[,factor_name])){
            dfpred[,factor_name]<-factor(dfpred[,factor_name])
        }

        ## return the ids of each
        index_fact<-lapply(levels(dfpred[,factor_name]), function(x) which(dfpred[,factor_name]==x))

        for (i in 1:length(index_fact)) {
            list_distri[[i]] <- lp[,index_fact[[i]]]
            colnames(list_distri[[i]])<-dfpred[,time_var][index_fact[[i]]]
        }

        ## names for each element of the list
        names(list_distri)<-levels(dfpred[,factor_name])

        ## this list is returned as the basis to compute mean trend and confidene intervals
        return(list_distri)
    }  else{ #if no structuring factor, then just name lp (with the year index) and return lp
        colnames(lp)<-dfpred[,time_var]
        return(lp)
    }

}
