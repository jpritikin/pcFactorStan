default <- itemModelExplorer.default

shinyUI(pageWithSidebar(
  headerPanel('Item model explorer'),
  sidebarPanel(
    sliderInput("numThresholds", "Thresholds:",
                min = 1, max = 10,
                value = length(default$thresholds), ticks = FALSE),
    checkboxInput("showParameters", label = "Show/Edit parameters", value = TRUE),
    conditionalPanel(
      condition = "input.showParameters",
      hr(),
      selectInput('editPar', 'Edit Parameter:', 'discrimination'),
      sliderInput('editParValue', "Parameter Value:",
                  min=-5, max=5, value=default$discrimination,
                  step=.01, ticks=FALSE),
      hr(),
      tableOutput("parView")
    )
  ),
  mainPanel(
    plotOutput('plot1')
  )
))
