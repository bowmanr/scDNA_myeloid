tabPanel(title = "Clonograph", 
         sidebarLayout(
           sidebarPanel(
             selectInput("clonoInput", "Sample Names:", names(sample_list),multiple = TRUE,selected="MSK45")
             
           ),
           mainPanel(
             fluidRow(
               uiOutput("plot.ui"),
               #downloadButton(outputId = "downloadClonograph", label = "Download your plot")
             )
           )
         )
         
         
)