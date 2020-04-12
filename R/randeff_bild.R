setMethod("randeff",
          signature(object = "bild"),
          function (object) 
          {
            coef <- object@coefficients
            bi.est<- object@bi.estimate
            nas <- is.na(coef[, 1])
            cnames <- names(coef[, 1][!nas])
            names(coef)<-cnames
            X<-object@model.matrix
 
            # for one random effect
            if ( !all(is.na(match(cnames, "omega"))))
              
            {  n.cases<- length(bi.est)
            id<-c(1:n.cases)
            tabela2<-as.data.frame(matrix(as.double(bi.est),nrow=n.cases,ncol=1))
            dimnames(tabela2) <- list(c(1:n.cases), "(Intercept)")
            return(tabela2)
            } 
            
          else 
              warning("\nOnly to Random effects model") 
            
        })