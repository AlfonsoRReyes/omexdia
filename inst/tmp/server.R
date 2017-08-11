
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyServer(function(input, output) {

  output$distPlot <- renderPlot({

    # generate bins based on input$bins from ui.R
    # x    <- faithful[, 2]
    # theBins <- seq(min(x), max(x), length.out = input$theBins + 1)

    # draw the histogram with the specified number of bins
    # hist(x, breaks = theBins, col = 'darkgray', border = 'white')
      omexdia_solved <- OMEXDIAsteady()
      
      

  })

})
