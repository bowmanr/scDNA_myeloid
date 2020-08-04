
<<<<<<< Updated upstream
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
               plotOutput("sampleClonPlot",width="50%"),
               #downloadButton(outputId = "downloadSampleClonPlot", label = "Download your plot")
             )
           )
         )
         
    
)
=======
       sidebarLayout(
         sidebarPanel(
           h3("Graphs:"),
           radioButtons("sc", "", 
                        c("Number of Mutations" = "oneC",
                          "Number of Clones" = "oneE",
                          "Shannon Diversity" = "twoA",
                          "Mutations in \n Dominant Clone" = "twoB",
                          "Dominant Clone Size" = "threeA"
                        )),
           conditionalPanel(
             condition = "input.tabselected==1",
             selectInput("sampleClonGroups", "Grouping",c("Final_group","Dx","Group") ,selected="Final_group")
           ),
           conditionalPanel(
             condition = "input.tabselected==2",
             splitLayout( 
               textInput("Group1",label="Group Name",value="DNMT3A_only"),
               selectInput("Group1_include", label="Include", 
                           colnames(clone_mutations)[6:33],
                           multiple = TRUE,
                           selected=c("DNMT3A")),
               selectInput("Group1_exclude", label="Exclude", 
                           colnames(clone_mutations)[6:33],
                           multiple = TRUE,
                           selected=c("RAS","FLT3"))),
             splitLayout( 
               textInput("Group2",label=NULL,value="DNMT3A_NRAS"),
               selectInput("Group2_include", label=NULL,
                           colnames(clone_mutations)[6:33],
                           multiple = TRUE,
                           selected=c("DNMT3A","NRAS")),
               selectInput("Group2_exclude", label=NULL,
                           colnames(clone_mutations)[6:33],
                           multiple = TRUE,
                           selected=c("FLT3"))),
             splitLayout( 
               textInput("Group3",label=NULL,value="DNMT3A_FLT3"),
               selectInput("Group3_include", label=NULL,
                           colnames(clone_mutations)[6:33],
                           multiple = TRUE,
                           selected=c("DNMT3A","FLT3")),
               selectInput("Group3_exclude", label=NULL,
                           colnames(clone_mutations)[6:33],
                           multiple = TRUE,
                           selected=c("NRAS"))),
         )),
         mainPanel(
           tabsetPanel(type = "tabs", id = "tabselected",
                       tabPanel("Premade", value=1, plotOutput("sampleClonPlotP",width="50%")),
                       tabPanel("Custom",  value=2, plotOutput("sampleClonPlotC",width="50%"))
           )
         )
       )
)

>>>>>>> Stashed changes
