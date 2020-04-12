setMethod(f="getcoef",
          signature=c(object = "bild"),
          function (object) 
          {
          return(object@coefficients)
          }
        )

