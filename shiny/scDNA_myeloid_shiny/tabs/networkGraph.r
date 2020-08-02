tabPanel(title = "Network Graph", 
         sidebarLayout(
           sidebarPanel(
             selectInput("networkInput", "Gene", colnames(clone_mutations)[6:33],multiple = TRUE,selected=c("DNMT3A","TET2","IDH1","IDH2","ASXL1"))
           ),
           mainPanel(
          #   fluidRow(
               plotOutput("networkPlot",width="100%",height="500px"),
               #downloadButton(outputId = "downloadNetworkGraph", label = "Download your plot")
            # )
           )
         )
         
)