setMethod("fixeff", 
          signature(object="bild"),
          function(object) 
          {
            coef <- object@coefficients
            nas <- is.na(coef[, 1])
            cnames <- names(coef[, 1][!nas])
         
            
            # for dependence="ind" 
            if(all(is.na(match(cnames, "log.psi1"))) && 
               all(is.na(match(cnames, "log.psi2")))  &&
               all(is.na(match(cnames, "omega"))))
            {  
              return(coef[ ,  ])
            }
            

            # for dependence="MC1" 
            else if(!all(is.na(match(cnames, "log.psi1"))) && 
                    all(is.na(match(cnames, "log.psi2")))  &&
                    all(is.na(match(cnames, "omega"))))
            {
              return(coef[1:(dim(coef)[[1]]-1) ,  ])
            }
            
            
            # for dependence="MC2" 
            else if(!all(is.na(match(cnames, "log.psi1"))) && 
                    !all(is.na(match(cnames, "log.psi2")))  &&
                    all(is.na(match(cnames, "omega"))))
            {
              return(coef[1:(dim(coef)[[1]]-2) ,  ])
            }
            
 
            # for dependence="indR" 
            else if (all(is.na(match(cnames, "log.psi1"))) && 
                     all(is.na(match(cnames, "log.psi2"))) &&
                     !all(is.na(match(cnames, "omega"))))
            {
              return(coef[-dim(coef)[[1]],  ])
            }
            
            # for dependence="MC1R" 
            else if (!all(is.na(match(cnames, "log.psi1"))) && 
                      all(is.na(match(cnames, "log.psi2")))  &&
                     !all(is.na(match(cnames, "omega"))))
            {
              return(coef[1:(dim(coef)[[1]]-2),  ])
            }

            
            # for dependence="MC2R" 
            else if (!all(is.na(match(cnames, "log.psi1"))) && 
                     !all(is.na(match(cnames, "log.psi2")))  &&
                     !all(is.na(match(cnames, "omega"))))
            {
              return(coef[1:(dim(coef)[[1]]-3),  ])
            } 
            
          }
        )