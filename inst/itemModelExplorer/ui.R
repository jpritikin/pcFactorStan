default <- itemModelExplorer.default

shinyUI(pageWithSidebar(
  headerPanel('Item model explorer'),
  sidebarPanel(
    conditionalPanel(
      "0==1",  # debugging only
      checkboxInput("newModel", label = "New Model", value = TRUE)),
    sliderInput("numThresholds", "Thresholds:",
                min = 1, max = 10,
                value = length(default$thresholds), ticks = FALSE),
    checkboxInput("showParameters", label = "Show/Edit parameters", value = TRUE),
    conditionalPanel(
      condition = "input.showParameters",
      hr(),
      selectInput('editPar', 'Edit Parameter:',
                  c('discrimination', 'scale', paste0('th',1:2))),
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
