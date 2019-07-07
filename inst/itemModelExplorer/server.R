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
  unsummed <- c(0, diff + rev(th), diff - th, use.names = FALSE)
  cumsum(unsummed * alpha)
}

calcProb <- function(par, theta) {
  sapply(theta, function(x) softmax(cmp_probs(par['discrimination'], par['scale'], x,
    par[-(1:2)])), USE.NAMES=FALSE)
}

shinyServer(function(input, output, session) {
  state <- reactiveValues(
    par = c(discrimination=2, scale=1, th1=.8, th2=1.6))

  observe({
    numThr <- input$numThresholds
    if (numThr + 2 == length(par)) return()
    oldPar <- isolate(state$par)
    par <- oldPar[1:3]
    if (numThr > 1) par <- c(par, rep(2*par[3], numThr-1))
    names(par)[3:length(par)] <- paste0('th',1:numThr)
    state$par <- par
    updateSelectInput(session, "editPar", choices = names(par),
                      selected = 'discrimination')
  })

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
    data.frame(name=names(state$par), par=state$par)
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

