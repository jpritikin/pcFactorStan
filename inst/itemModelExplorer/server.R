library(ggplot2)
library(reshape2)

#options(shiny.reactlog=TRUE)

verbose <- FALSE
#verbose <- TRUE

moveSomeParameter <- function(input, state, whichPar, newValue) {
  if (verbose) cat("moveParameter", whichPar,"to",newValue, fill=T)
  oldValue <- isolate(state$par[whichPar])
  doit <- newValue != oldValue
  if (doit) {
    state$par[whichPar] <- newValue
  }
  return(doit)
}

moveParameter <- function(input, state, newValue) {
  whichPar <- isolate(input$editPar)
  moveSomeParameter(input, state, whichPar, newValue)
}

softmax <- function(y) exp(y) / sum(exp(y))

cmp_probs <- function(alpha, scale, rawDiff, thRaw) {
  th <- cumsum(thRaw)
  diff = scale * rawDiff
  unsummed <- c(0, c(diff + rev(th)), c(diff - th), use.names = FALSE)
  cumsum(unsummed * alpha)
}

calcProb <- function(par, theta) {
  sapply(theta, function(x) softmax(cmp_probs(par['discrimination'], par['scale'], x,
    c(par['th1'], par['th2']))), USE.NAMES=FALSE)
}

shinyServer(function(input, output, session) {
  par <- c(discrimination=1.749, scale=1, th1=.8, th2=1.7)
  state <- reactiveValues(par = par)
  
  parNames <- names(par)
  updateSelectInput(session, "editPar", choices = parNames, selected = parNames[1])

  observe({
     newVal <- input$editParValue
     if (!is.numeric(newVal)) return()
     moveParameter(input, state, newVal)
   })
  
  observe({
    whichPar <- input$editPar
    val <- isolate(state$par[whichPar])
    if (verbose) cat("switch editing parameter to", whichPar, val, fill = TRUE)
    updateSliderInput(session, "editParValue",
                      label = paste("Parameter", whichPar), value = as.numeric(val))
  })

  output$parView <- renderTable({
    data.frame(par=state$par)
  })
  
  output$plot1 <- renderPlot({
    width <- 4
    grid <- expand.grid(theta=seq(-width,width,.1))
    
    trace <- try(calcProb(state$par, grid$theta), silent = TRUE)
    if (inherits(trace, "try-error") || any(is.na(trace))) {
      pl <- ggplot(grid, aes(theta, 0)) + geom_line() + ylim(0,1) +
        geom_text(label="Invalid parameters", y=.5, x=0, size=14, color="red")
      return(pl)
    }
    grid <- cbind(grid, t(trace))
    grid2 <- melt(grid, id.vars=c("theta"), variable.name="category", value.name="p")
    
    ggplot(grid2, aes(theta, p, color=category)) + geom_line() +
      ylim(0,1) + xlim(-width, width)  + theme(legend.position="none")
  })
})

# runApp('.',display.mode="showcase")

