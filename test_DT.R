library(DT)    #  for data tables

library(shiny)

shinyApp(
  ui = fluidPage(fluidRow(column(12, DTOutput('tbl')))),
  server = function(input, output) {
    output$tbl = renderDT(
      iris, options = list(lengthChange = FALSE)
    )
  }
)