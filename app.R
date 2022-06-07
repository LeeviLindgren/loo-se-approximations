library(purrr)
library(shiny)
library(ggplot2)

results <- readRDS('simres.Rds')
params_df <- map_df(results, 'sim_params')


ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectInput('n', 'Sample size', choices = unique(params_df$n)),
      selectInput('p', 'Number of predictors', choices = unique(params_df$p)),
      selectInput('sigma', 'Observation noise', choices = unique(params_df$sigma))
    ),
    mainPanel(plotOutput('plot'))
  )
)

server <- function(input, output, session) {
  output$plot <- renderPlot({
    n <- input$n
    p <- input$p
    sigma <- input$sigma
    
    tmp <- results %>%
      keep(function(x) {
        x['sim_params'][[1]]$n == n &
        x['sim_params'][[1]]$p == p &
        x['sim_params'][[1]]$sigma == sigma
      })
    
    
    ggplot(tmp[[1]]$df, aes(sd_taylor, sd_bb, color = R2)) + 
      geom_point() + 
      geom_abline() + 
      xlab('standard error taylor approximation') + ylab('standard error bayesian bootstrap')
  })
}

shinyApp(ui, server)




