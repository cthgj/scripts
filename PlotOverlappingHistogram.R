library(ggplot2)
nolim<-data.frame("./no_lim.clasnum")
oneh<-data.frame("./100.clasnum")
twoh<-data.frame("./200.clasnum")
nolim$aa<-'noLimit'



plotOverlappingHist <- function(a, b, colors=c("blue", "blue"),
                            breaks=NULL, xlim=NULL, ylim=NULL){

  ahist=NULL
  bhist=NULL
 
  if(!(is.null(breaks))){
    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
  } else {
    ahist=hist(a,plot=F)
    bhist=hist(b,plot=F)
    dist = ahist$breaks[2]-ahist$breaks[1]
    breaks = seq(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks),dist)

    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
  }

  if(is.null(xlim)){
    xlim = c(min(ahist$breaks,bhist$breaks),max(ahist$breaks,bhist$breaks))
  }

  if(is.null(ylim)){
    ylim = c(0,max(ahist$counts,bhist$counts))
  }

  overlap = ahist
  for(i in 1:length(overlap$counts)){
    if(ahist$counts[i] > 0 & bhist$counts[i] > 0){
      overlap$counts[i] = min(ahist$counts[i],bhist$counts[i])
    } else {
      overlap$counts[i] = 0
    }
  }

  plot(ahist, xlim=xlim, ylim=ylim, col=colors[1])
  plot(bhist, xlim=xlim, ylim=ylim, col=colors[2], add=T)

}


plotOverlappingHist3 <- function(a, b, c, colors=c("green","blue", "red", "white"),
                            breaks=NULL, xlim=NULL, ylim=NULL){

  ahist=NULL
  bhist=NULL
  chist=NULL
  if(!(is.null(breaks))){
    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
    chist=hist(c,breaks=breaks,plot=F)
  } else {
    ahist=hist(a,plot=F)
    bhist=hist(b,plot=F)
    chist=hist(c,plot=F)
    dist = ahist$breaks[2]-ahist$breaks[1]
    breaks = seq(min(ahist$breaks,bhist$breaks, chist$breaks),max(ahist$breaks,bhist$breaks, chist$breaks),dist)

    ahist=hist(a,breaks=breaks,plot=F)
    bhist=hist(b,breaks=breaks,plot=F)
    chist=hist(c,breaks=breaks,plot=F)
  }

  if(is.null(xlim)){
    xlim = c(min(ahist$breaks,bhist$breaks, chist$breaks),max(ahist$breaks,bhist$breaks, chist$breaks))
  }

  if(is.null(ylim)){
    ylim = c(0,max(ahist$counts,bhist$counts,chist$counts))
  }

  overlap = ahist
  for(i in 1:length(overlap$counts)){
    if(ahist$counts[i] > 0 & bhist$counts[i] > 0  & chist$counts[i] > 0){
      overlap$counts[i] = min(ahist$counts[i],bhist$counts[i],chist$counts[i])
    } else {
      overlap$counts[i] = 0
    }
  }

  plot(ahist, xlim=xlim, ylim=ylim, col=colors[1])
  plot(bhist, xlim=xlim, ylim=ylim, col=colors[2], add=T)
  plot(chist, xlim=xlim, ylim=ylim, col=colors[3], add=T)


}
