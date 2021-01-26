#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)


# Define UI for application that draws a histogram
shinyUI(fluidPage(
    theme = shinytheme("superhero"),
    # Application title
    titlePanel("Models For Sorption Isotherms"),
    h4("Piyadi G. J. Lakshika"),
    h4("Statistical Consultanct Unit of USJP, 2020"),
    br(),
    br(),
    br(),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            #sliderInput("bins",
            #            "Number of bins:",
            #            min = 1,
            #            max = 50,
            #            value = 30),
            fileInput("file","Upload the file"), #fileInput() function is used to upload a file
            #helpText("Default file size is 5MB"),
            helpText("Select the parameters below"),
            checkboxInput(inputId = 'header', label = "Header", value = TRUE),
            radioButtons(inputId = 'sep', label = 'Separator', choices = c(Comma=',',Semicolon=';',Tab='\t', Space = ''), selected = ','),
            uiOutput("selectModel")
            
        ),
        
        
        

        # Show a plot of the generated distribution
        mainPanel(
           # plotOutput("distPlot")
            uiOutput("tb")
        )
    )
))
