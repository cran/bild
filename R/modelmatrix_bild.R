setMethod("model.mat",
    signature(object = "bild"),
    function (object) 
    {return(object@model.matrix)})