library(shiny)
library(ggplot2)

ui <- tagList(
  # Full width intro text outside fluidPage container
  tags$div(
    style = "width:100%; padding: 20px 40px; background:#f9f9f9; box-sizing: border-box;",
    h3("Introduction"),
    p(HTML("Michaelis–Menten kinetics describes how the rate of an enzyme-catalyzed reaction depends on substrate concentration, and two key parameters: the maximum reaction rate (V<sub>max</sub>) and the substrate concentration at which the reaction rate is half of V<sub>max</sub> (K<sub>m</sub>). This framework provides a useful analogy for understanding how Rubisco activity is influenced by CO<sub>2</sub> availability in C<sub>3</sub> plants.")), 
    br(),
    p(HTML("This app demonstrates the Michaelis-Menten enzyme kinetics model. A model is a simplified mathematical description of a process - in this case, how an enzyme's reaction rate depends on substrate concentration.")),
    br(),
    p(HTML("Adjust the kinetic parameters using the sliders below. The model used is: V = (V<sub>max</sub> × [S]) / (K<sub>m</sub> + [S])")),
    br()
  ),
  
  fluidPage(
   # titlePanel("Michaelis-Menten Kinetics"),
    
    fluidRow(
      column(6, offset = 3,
             plotOutput("mmPlot"),
             br(),
             fluidRow(
               column(6,
                      sliderInput("Vmax", "Vmax (μmol/min/mg enzyme)",
                                  min = 0.5, max = 5, value = 1.5, step = 0.1)
               ),
               column(6,
                      sliderInput("Km", "Km (mM)",
                                  min = 0.01, max = 0.5, value = 0.05, step = 0.01)
               )
             )
      )
    )
  )
)

server <- function(input, output) {
  michaelis_menten <- function(S, Vmax, Km) {
    Vmax * S / (Km + S)
  }
  
  default_Vmax <- 1.5
  default_Km <- 0.05
  substrate <- seq(0, 0.5, by = 0.01)
  default_velocities <- michaelis_menten(substrate, default_Vmax, default_Km)
  default_df <- data.frame(S = substrate, V = default_velocities)
  
  output$mmPlot <- renderPlot({
    current_velocities <- michaelis_menten(substrate, input$Vmax, input$Km)
    current_df <- data.frame(S = substrate, V = current_velocities)
    
    ggplot() +
      geom_line(data = default_df, aes(x = S, y = V),
                color = "lightblue", size = 1.2, alpha = 0.9) +
      geom_line(data = current_df, aes(x = S, y = V),
                color = "seagreen", size = 1.5) +
      labs(
        x = "Substrate concentration [S] (mM)",
        y = "Reaction rate, V (μmol/min/mg enzyme)"
      ) +
      coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 2)) +
      theme_classic()
  })
}

shinyApp(ui, server)
