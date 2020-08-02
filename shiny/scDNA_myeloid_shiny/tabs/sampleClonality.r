
tabPanel(title = "Sample Clonality", 

       sidebarLayout(
         sidebarPanel(
           #choose the graph
           radioButtons("sc", "Sample Clonality Graphs:", 
                        c("Figure 1C" = "oneC",
                          "Figure 1E" = "oneE",
                          "Figure 2A" = "twoA",
                          "Figure 2B" = "twoB",
                          "Figure 3A" = "threeA"
                        )),
           #choose the x-axis
           selectInput("selected_group", "Grouping",c("Final_group","Dx","Group") ,selected="Final_group")
           
         ),
         mainPanel(
           fluidRow(
             plotOutput("sampleClonPlot",width="50%")
             #downloadButton(outputId = "downloadSampleClonPlot", label = "Download your plot")
           )
         )
       )
    
)

