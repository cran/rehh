distribplot <- function(data,lty = 1,lwd = 1.5,col = c("blue","red"),main = "Genome-wide distribution",xlab = "",cex.main = 1.5,cex.lab = 1.25,qqplot = TRUE) {
	if (!is.null(names(dev.list())) && ((names(dev.cur()) == "windows") | (names(dev.cur()) == "X11") | (names(dev.cur()) == "quartz"))) {
   		dev.new()
    }
	par(mar = c(5,5,4,2) + 0.1)
	plot(density(data,na.rm = TRUE),main = main,xlab = xlab,col = col[1],lty = lty,lwd = lwd,cex.main = cex.main,cex.lab = cex.lab)
	curve(dnorm,col = col[2],add = TRUE)
	legend("topright",c("Observed","Gaussian"),bty = "n",col = col,lty = lty,lwd = lwd)
	if (qqplot) {
	    if (!is.null(names(dev.list())) && ((names(dev.cur()) == "windows") | (names(dev.cur()) == "X11") | (names(dev.cur()) == "quartz"))) {
   	    		dev.new()
    	}
		par(mar = c(5,5,4,2) + 0.1)
		qqnorm(data[!is.na(data)],cex.main = cex.main,cex.lab = cex.lab,pch = 16,cex = 0.75)
		abline(a=0,b=1,lty=2)
	}
	par(mar = c(5,4,4,2) + 0.1)
}
