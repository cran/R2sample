allmethods <- list(cont=c("KS", "Kuiper", "CvM", "AD",  "LR", "ZA", "ZK", "ZC", "Wassp1", "ES large", "ES small", "EP large", "EP small"),
                   disc= c("KS", "Kuiper", "CvM", "AD",  "LR", "ZA", "Wassp1", "large", "small"))

shinyUI(fluidPage(
  titlePanel("Twosample Tests"),
  HTML("<h4>Enter all the information required and then hit Go.</h4>"), 
  fluidRow(
      column(2, actionButton("gobutton",HTML("<h4><font color=\"red\">Go<font color=\"black\"></h4>"))), 
      column(2, textInput("B", "# of Simulations", value="5000", width="85%")),
      column(3, textInput("nbins", "Number of Bins", value="50, 10", width="85%")),
      column(3, textInput("maxProcessor", "Number of Processors to Use", value="10", width="85%"))
  ), 
  HTML("<hr>"), 
  fluidRow(
    conditionalPanel( condition = "input.datatype == 'Continuous'",
           column(10, checkboxGroupInput("cmethods", HTML("<h4>Methods</h4>"), 
                allmethods$cont, inline=TRUE, selected=allmethods$cont))
    ),       
    conditionalPanel( condition = "input.datatype == 'Discrete'",
           column(10, checkboxGroupInput("dmethods", HTML("<h4>Methods</h4>"), 
                allmethods$disc, inline=TRUE, selected=allmethods$disc))
    )    
  ),  
  HTML("<hr>"),  
  fluidRow(
     column(3, selectInput("datatype", "Data is ...",choices = c("Continuous", "Discrete"))),
     conditionalPanel( condition = "input.datatype == 'Continuous'",
        HTML("<h4>Two files for the two data sets</h4>"),
        column(3, fileInput('fileData1', 'Upload file with data set 1')),
        column(3, fileInput('fileData2', 'Upload file with data set 2'))   
     ),
     conditionalPanel( condition = "input.datatype == 'Discrete'",
        HTML("<h4>File should have three columns for the x counts,  y counts and values, separated by an empty spaces</h4>"),
        column(3, fileInput('fileData3', 'Upload file'))
     )
  ),
  column(5, tableOutput("tblTest"))  
))