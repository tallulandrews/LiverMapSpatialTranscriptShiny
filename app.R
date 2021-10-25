#### Packages

library(shiny)

library(Seurat)
slice_obj <- readRDS("data/C73_C1_RESEQ_raw.rds")
all_genes <- sort(unique(rownames(slice_obj)))
all_features <- colnames(slice_obj@meta.data)[2:ncol(slice_obj@meta.data)]

not_gene_token <- "(None)"

#### UI
ui <- fluidPage(
	# App title ----
	titlePanel("Human Liver Atlas Spatial Transcriptomics"),

	# Sidebar Layout with input and output definitions ----
	sidebarLayout (

		# Sidebar panel for inputs ----
		sidebarPanel(
		
			# Input: Slider for transparency of the gene expression
			sliderInput("alpha",
				label = "Plot Opacity",
                       	min = 0, 
				max = 100, 
				value = 100),
	
			# Input: Select Gene to display
			selectInput(inputId = "gene",
	                  label = "Choose a gene:",
                  	choices = c(not_gene_token, all_genes), 
				selected = not_gene_token),

			# Input: Select Gene to display
			selectInput(inputId = "feature",
	                  label = paste('Choose a feature (Set gene to:"', not_gene_token, '")', sep =""),
                  	choices = all_features, 
				selected = "zonation_score"),

			# Input: Set resolution of image to download
			sliderInput("res",
				label = "Download Image Resolution (dpi)",
                       	min = 50, 
				max = 500, 
				step = 50,
				value = 300),

			# Download Button
			downloadButton("downloadPlot", "Download Plot")
			
		),
	
		# Main panel for displaying outputs ----
		mainPanel(

			# Output: SpatialFeaturePlot
			plotOutput(outputId = "sptFeaturePlot")

		)
	)
)

#### Server
server <- function(input, output) {

	# Seurat Spatial Feature Plot of the specified gene and slice.
	plotImage <- function() {
		toplot <- input$gene
		if (input$gene == not_gene_token) {
			toplot <- input$feature
		}
		print(suppressWarnings(Seurat::SpatialFeaturePlot( slice_obj, features=toplot, alpha=input$alpha/100)))
	}


	# renderPlot = react to input & output is a plot
	output$sptFeaturePlot <- renderPlot({
		plotImage()
	})

	output$downloadPlot <- downloadHandler(
		filename = paste("Spatial", input$gene, input$slice, "plot.png", sep="_"),
		content = function(file) {
		png(file, width=8, height=8, units="in", res=input$res)
		plotImage()
		dev.off()
    })  

}

#### shinyApp Call
shinyApp(ui = ui, server = server)