shinyUI(pageWithSidebar(
  headerPanel('Item model explorer'),
  sidebarPanel(
    checkboxInput("showParameters", label = "Show/Edit parameters", value = TRUE),
    conditionalPanel(
      condition = "input.showParameters",
      hr(),
      selectInput('editPar', 'Edit Parameter:', 'discrimination'),
      sliderInput('editParValue', "Parameter Value:",
                  min=-5, max=5, value=1.749, step=.01, ticks=FALSE),
      hr(),
      tableOutput("parView")
    )
  ),
  mainPanel(
    plotOutput('plot1')
  )
))
