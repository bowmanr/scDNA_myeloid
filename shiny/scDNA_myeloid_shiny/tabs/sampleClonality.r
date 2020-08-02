tabPanel(title = "Sample Clonality", 

       sidebarLayout(
         sidebarPanel(
           #choose the graph
           radioButtons("sc", "Premade Graphs:", 
                        c("Figure 1C" = "oneC",
                          "Figure 1E" = "oneE",
                          "Figure 2A" = "twoA",
                          "Figure 2B" = "twoB",
                          "Figure 3A" = "threeA"
                        )),
           #choose the x-axis
           h2("Customization Tools:"),
           selectInput("selected_group", "Customization", available_groups,multiple = TRUE,selected=c("DTAI"))
           
         ),
         mainPanel(
           tabsetPanel(type = "tabs", 
                       tabPanel("Premade", plotOutput("sampleClonPlotP",width="50%")),
                       tabPanel("Custom", plotOutput("sampleClonPlotC",width="50%"))
             
           )
         )
       )
    
)

