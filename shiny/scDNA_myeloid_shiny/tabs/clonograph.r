tabPanel(title = "Clonograph", 
         sidebarLayout(
           sidebarPanel(
             selectInput("clonoInput", "Select Gene of Interest:", colnames(sample_list)[1:9],multiple = TRUE,selected="MSK45"),
             selectInput("clonoInput", "Select Patients With Gene:", names(sample_list),multiple = TRUE,selected="MSK45")
             
             
           ),
           mainPanel(
             fluidRow(
               uiOutput("plot.ui"),
               #downloadButton(outputId = "downloadClonograph", label = "Download your plot")
             )
           )
         )
         
         
)