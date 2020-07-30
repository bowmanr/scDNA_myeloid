tabPanel(title = "Clonograph", 
         sidebarLayout(
           sidebarPanel(
             h2("customization options...")
             
           ),
           mainPanel(
             fluidRow(
               column(plotOutput(network2e)),
               column(plotOutput(network2f))
             )
           )
         )
         
         
)