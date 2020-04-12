setMethod("vareff",
          signature(object = "bild"),
          function (object) 
          {
            coef <- object@coefficients
            nas <- is.na(coef[, 1])
            cnames <- names(coef[, 1][!nas])

            # for dependence="indR" 
            if (all(is.na(match(cnames, "log.psi1"))) && !all(is.na(match(cnames, "omega"))))
            {
              coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
             # cat("\nRandom effects:\t\n") 
              table.r1<-data.frame(coef[dim(coef)[[1]], 1])
              dimnames(table.r1)<-list("(Intercept)", "Variance")
              return(table.r1)
            }
            

            # for dependence="MC1R" 
            else if (!all(is.na(match(cnames, "log.psi1"))) && 
                     all(is.na(match(cnames, "log.psi2")))  &&
                     !all(is.na(match(cnames, "omega"))))
            {
              coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
             # cat("\nRandom effects:\t\n") 
              table.r1<-data.frame(coef[dim(coef)[[1]], 1])
              dimnames(table.r1)<-list("(Intercept)", "Variance")
              return(table.r1)
            }
            
            # for dependence="MC2R" 
            else if (!all(is.na(match(cnames, "log.psi1"))) && 
                     !all(is.na(match(cnames, "log.psi2")))  &&
                     !all(is.na(match(cnames, "omega"))))
            {
              coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
              cat("\nRandom effect:\t\n") 
              table.r1<-data.frame(coef[dim(coef)[[1]], 1])
              dimnames(table.r1)<-list("(Intercept)", "Variance")
              print(table.r1)
            }
              
            else 
              warning("\nOnly to Random effects model") 
              
          }
        )