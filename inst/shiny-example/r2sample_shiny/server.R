red <- paste("<font color=\"red\">")
library(ggplot2)
library(Rcpp)
library(grid)
library("R2sample")
#methods <- c("chi large", "chi small", "t test", "KS", "Kuiper", "AD", "CdM", "LR", "ZK", "ZA", "ZC", "Wassp1")

shinyServer(function(input, output,session) {

   indata <- reactive({  
      if(input$datatype=="Continuous") {
          inFile1 <- input$fileData1
          if (is.null(inFile1)) return(NULL)
          inFile2 <- input$fileData2
          if (is.null(inFile2)) return(NULL)
          return(list(
            x=scan(inFile1$datapath),
            y=scan(inFile2$datapath)
          ))
      }          
      if(input$datatype=="Discrete") {
          inFile3 <- input$fileData3
          if (is.null(inFile3)) return(NULL)
          tmp=read.table(inFile3$datapath)
          return(list(x=tmp[,1],y=tmp[,2],vals=tmp[,3]))
      }

    })


    
    

    test=eventReactive(input$gobutton, {
        z=indata() 
        nbins=as.numeric(strsplit(input$nbins,",")[[1]]) 
        if(input$datatype=="Continuous") {
            methods=input$cmethods
            out=R2sample::twosample_test(z$x, z$y, 
                      B=as.numeric(input$B),
                      nbins=nbins, 
                      maxProcessor=as.numeric(input$maxProcessor),
                      doMethod=methods)                     
            methods=names(out$statistics)          
            out=cbind(out$statistics, out$p.values)
            colnames(out)=c("Statistics", "p - value")
            rownames(out)=methods
        }    
        if(input$datatype=="Discrete") {
            out=R2sample::twosample_test(z$x, z$y, z$vals, 
                     B=as.numeric(input$B),
                     nbins=nbins,
                     maxProcessor=as.numeric(input$maxProcessor),
                     doMethod=input$dmethods)
            out=cbind(out$statistics, out$p.values)
            colnames(out)=c("Statistics", "p - value")
            rownames(out)=input$dmethods
        }                 
        out       
    })
    
    output$tblTest <- renderTable({ test()  }, striped=TRUE,rownames=TRUE,colnames=TRUE)
    
    
})     
