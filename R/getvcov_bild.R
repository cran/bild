setMethod("getvcov",
    signature(object = "bild"),
    function (object) 
    {
      cov<-object@covariance
      cnames <- rownames(cov)
      r1<-nrow(cov)
      c1<-ncol(cov)
 
      # for dependence="ind" or dependence="AR1"
       if(all(is.na(match(cnames, "omega")))  )
        {     cov
          }
      
      
      # for dependence="indR" or dependence="MC1R" or dependence="MC2R"
      else if (!all(is.na(match(cnames, "omega"))))
        {  cov.aux<-cov[1:(r1-1),1:(c1-1)]
           cov.aux
        }
        
} )