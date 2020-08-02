tabPanel(title = "Sample Clonality", 

       sidebarLayout(
         sidebarPanel(
           #choose the graph
           h3("Premade Graphs:"),
           radioButtons("sc", "", 
                        c("Figure 1C" = "oneC",
                          "Figure 1E" = "oneE",
                          "Figure 2A" = "twoA",
                          "Figure 2B" = "twoB",
                          "Figure 3A" = "threeA"
                        )),
           #choose the x-axis
           h3("Customization Tools:"),
           selectInput("sampleClonGroups", "Select Groups", available_groups, multiple = TRUE,selected=c("DTAI"))
           
         ),
         mainPanel(
           tabsetPanel(type = "tabs", 
                       tabPanel("Premade", plotOutput("sampleClonPlotP",width="50%")),
                       tabPanel("Custom", plotOutput("sampleClonPlotC",width="50%"))
             
           )
         )
       )
    
)

