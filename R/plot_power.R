#'  plot_power
#'
#' This function draws the power graph, with curves sorted by the mean power and smoothed for easier reading.
#' @param  pwr  a matrix of power values, usually from the twosample_power command 
#' @param  Smooth =TRUE lines are smoothed for easier reading
#' @param  xname Name of variable on x axis
#' @return plt, an object of class ggplot.
#' @export 

plot_power=function(pwr, xname=" ", Smooth=TRUE) {

# sort methods by their average power
  mu=apply(pwr,2,mean) 
  lvls = colnames(pwr)[order(mu, decreasing = TRUE)]
# create data frame  
  df=data.frame(x=rep(as.numeric(rownames(pwr)), ncol(pwr)),
       y=100*c(pwr),
       Method=factor(rep(colnames(pwr), each=nrow(pwr)),  
                     levels=lvls,
                     ordered = TRUE))
# create ggplot graphics                     
plt=ggplot2::ggplot(data=df, ggplot2::aes(df$x, df$y, color=df$Method))+
     ggplot2::xlab(xname)+
     ggplot2::ylab("Power")+
     ggplot2::scale_color_manual(labels=lvls,
                                 values=1:12,
                                 name='Method')
 if(Smooth) plt=plt+
    ggplot2::geom_smooth(formula = y ~ x, method="loess", se=FALSE, span=0.25)
 else  plt=plt+ggplot2::geom_line()
 plt
}