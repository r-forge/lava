###{{{ Misc
Col <- function (col, alpha = 0.2) 
{
    sapply(col, function(x) do.call(rgb, as.list(c(col2rgb(x)/255, 
        alpha))))
}
###}}} Misc
