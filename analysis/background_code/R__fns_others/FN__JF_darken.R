## https://gist.github.com/Jfortin1/72ef064469d1703c6b30)
# JFortin1
darken <- function(color, factor=1.4){
	col <- col2rgb(color)
	col <- col/factor
	col <- rgb(t(col), maxColorValue=255)
	col
       }

