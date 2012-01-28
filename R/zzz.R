
.onLoad <- .First.lib <- function(libname,pkgname)
{
library.dynam("bild",pkgname,libname)
invisible()
}
