tabPanel(title = "Clonograph", 
         sidebarLayout(
<<<<<<< Updated upstream
           sidebarPanel(
             selectInput("clonoInput", "Sample Names:", names(sample_list),multiple = TRUE,selected="MSK45")
             
           ),
=======
            sidebarPanel(
              radioButtons("selection_Feature", "Plot Type",
                           c(Sample = "Sample",
                             Gene = "Gene")),
                   conditionalPanel(
                     condition = "input.selection_Feature == 'Gene'",
                         selectInput("clonoInputGene", 
                                     "Select Gene of Interest:", 
                                     colnames(sample_mutations)[c(2:27,30)],
                                      multiple = TRUE,selected=c("IDH2","DNMT3A")) ,
                      ) ,
                     conditionalPanel(
                       condition = "input.selection_Feature == 'Sample'",
                         selectInput("clonoInputSample", 
                                     "Select Patients:",
                                     names(sample_list),
                                     multiple = TRUE,selected="MSK45"),
                     )
               ),
>>>>>>> Stashed changes
           mainPanel(
             fluidRow(
               uiOutput("plot.ui"),
               #downloadButton(outputId = "downloadClonograph", label = "Download your plot")
             )
           )
         )
)