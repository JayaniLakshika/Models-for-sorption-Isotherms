#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(ggplot2) #To draw plots
library(boot) #To run bootstrap regression



shinyServer(function(input, output) {
    
    data <- reactive({
        if(is.null(input$file)){return()} 
        read.csv(input$file$datapath, header=input$header, sep=input$sep)
    })
    
    # Output a data table for the upload tab page
    output$contents <- renderTable({
        if(is.null(input$file)){return()}
        read.csv(input$file$datapath, header=input$header, sep=input$sep)
        
    })
    
    output$selectModel <- renderUI({
        if(is.null(input$file)){return()}
        list(hr(),
             helpText("Select the model"),
             selectInput("Select","Select", choices = c("All","Brunauer Emmett Teller (BET) Model","Oswin Model","Smith Model","Halsey Model","Guggenheim Anderson de Boer (GAB) Model"))
        )
        
    })
    
    
    # Function to obtain regression weights
    bs1 <- function(formula, data, indices) {
        d <- data[indices,] # allows boot to select sample
        fit <- lm(formula, data=d)
        return(c(coef(fit),summary(fit)$r.square))
    }
    
    
    # Function to obtain regression weights
   
    Smith<-function(){
        
        if(is.null(input$file)){return()}
        data <- read.csv(input$file$datapath, header=input$header, sep=input$sep)
        data$x <- log10(1-data$aw) #To add small value to neglect the infinity case
        results <- boot(data=data, statistic=bs1,
                        R=1000, formula=Mw ~ x)# Bootstrapping with 1000 replications
        coef <- as.data.frame(results[2]) # To get coefficients of the bootstrap results
        names(coef) <- c("Beta_0", "Beta_1","R_square")
        
        c1 <- median(coef$Beta_0[!is.na(coef$Beta_0)])# Estimate for c1
        #print(paste0("Estimate for C1 is ",round(c1,2)))
        c2 <- median(coef$Beta_1[!is.na(coef$Beta_1)]) * (-1) # Estimate for c2
        #print(paste0("Estimate for C2 is ",round(c2,2)))
        
        data$y_hat <- c1 - c2 * log10(1-data$aw)
        data$dif <- data$Mw -data$y_hat#To get the difference between actual and predicted values(Mw)
        MSE <- sum((data$dif)^2)/nrow(data)#To get the Mean Squared Error
        #print(paste0("Mean Squared Error(MSE) of the Smith model is ",round(MSE,2)))
        
        rsq <- median(coef$R_square) * 100#To get the R squared
        #print(paste0("R squared value is ",round(rsq,2), "%"))
        output <- data.frame(c1 = c1, 
                             c2=c2, 
                             R_square=rsq,
                             MSE = MSE)
        return(output)
        
    }
    
    
    Halsey<-function(){
        
        if(is.null(input$file)){return()}
        data <- read.csv(input$file$datapath, header=input$header, sep=input$sep)
    
        #add the column Y to dataset by transforming Mw as given in linearized equation
        data$Y<-log(data$Mw)
        
        
        #add the column X to dataset by transforming aw as given in linearized equation
    
        data$X<-log(-log(data$aw))
        
        #repeat 1000 times
        results<-boot(data=data,statistic=bs1,R=1000,formula=Y~X)
        
        coef <- as.data.frame(results[2]) # To get coefficients of the bootstrap results
        names(coef) <- c("C_hat", "m_hat","R_square")
        
        #median as the point estimate for m_hat(slope)
        m<-median(coef$m_hat[!is.na(coef$m_hat)])
        
        #median as the point estimate for C_hat(intercept)
        C<-median(coef$C_hat[!is.na(coef$C_hat)])
        
        #model constants
        n<-(-1/m)
        c<-exp(n*C)
        
        data$y_hat <- m*data$X+C
        data$dif <- data$Y -data$y_hat
        
        MSE <- sum((data$dif)^2)/nrow(data)#To get the Mean Squared Error
        #print(paste0("Mean Squared Error(MSE) of the Smith model is ",round(MSE,2)))
        
        rsq <- median(coef$R_square) * 100#To get the R squared
        #print(paste0("R squared value is ",round(rsq,2), "%"))
        output <- data.frame(n = n, 
                             c = c, 
                             R_square=rsq,
                             MSE = MSE)
        return(output)
        
        
        
    }
    
    
    Bet <- function(){
        
        if(is.null(input$file)){return()}
        data <- read.csv(input$file$datapath, header=input$header, sep=input$sep)
        data$x <- data$aw
        data$y <- data$aw /((1-data$aw)*data$Mw)
        
        results <- boot(data=data, statistic=bs1,
                        R=1000, formula=y~x)
        coef <- as.data.frame(results[2]) # To get coefficients of the bootstrap results
        names(coef) <- c("Beta_0", "Beta_1","R_square")
        
        
        C <- (median(coef$Beta_1[!is.na(coef$Beta_1)])/median(coef$Beta_0[!is.na(coef$Beta_0)])) + 1
        #print(paste0("Estimate for C is ",round(C,2)))
        M0 <- 1/(median(coef$Beta_0[!is.na(coef$Beta_0)])*C)
        #print(paste0("Estimate for M0 is ",round(M0,2)))
        
        #bet_fun <- function(x) M0*C*x/((1-x)*(1+x*(C-1)))
        
        data$y_hat <- M0*C*data$aw/((1-data$aw)*(1+data$aw*(C-1)))
        data$dif <- data$y -data$y_hat
        MSE <- sum((data$dif)^2)/nrow(data)#find the Mean Squared Error
        
        #print(paste0("Mean Squared Error(MSE) of the bet model is ",round(MSE,2)))
        
        rsq <- median(coef$R_square) * 100#To get the R squared
        
        #print(paste0("R squared value is ",round(rsq,2), "%"))
        output <- data.frame(C = C, 
                             M0 = M0, 
                             R_square=rsq,
                             MSE = MSE)
        return(output)
        
        
    }
    
    
    oswin <- function()
    {
        if(is.null(input$file)){return()}
        data <- read.csv(input$file$datapath, header=input$header, sep=input$sep)
        
        data$y<-log(data$Mw)
        data$x<-log(data$aw/(1-data$aw))
        
        
        #Number of repetition times = 1000
        results<-boot(data = data,statistic = bs1,R=1000,formula=y~x)
        coef <- as.data.frame(results[2]) # To get coefficients of the bootstrap results
        names(coef) <- c("Beta_0", "Beta_1","R_square")
        
        
        
        #Obtaining median of intercepts as point estimates
        c <- median(coef$Beta_0[!is.na(coef$Beta_0)])
        
        #Obtaining median of gradients as point estimates
        m <- median(coef$Beta_1[!is.na(coef$Beta_1)])
        
        #Calculating the coeficients
        
        n=m
        C = exp(c)
        
        
        #Fitted values
        data$y_hat<-c + m *data$x
        
        data$dif <- data$y -data$y_hat
        MSE <- sum((data$dif)^2)/nrow(data)#find the Mean Squared Error
        
        rsq <- median(coef$R_square) * 100
        output <- data.frame(n = n, 
                             C = C, 
                             R_square=rsq,
                             MSE = MSE)
        return(output)
        
    }
    
    
    GAB <- function(){
        
        if(is.null(input$file)){return()}
        data <- read.csv(input$file$datapath, header=input$header, sep=input$sep)
        data$x <- data$aw
        data$y <- data$aw/data$Mw
        
        results <- boot(data=data, statistic=bs1,
                        R=1000, formula=y~I(x^2)+x)
        coef <- as.data.frame(results[2]) # To get coefficients of the bootstrap results
        names(coef) <- c("Beta_0", "Beta_1","Beta_2","R_square")
        
        #print(paste0("Estimate for C is ",round(C,2)))
        M0 <- 1/sqrt(median(coef$Beta_1[!is.na(coef$Beta_1)])^2-(4*median(coef$Beta_0[!is.na(coef$Beta_0)])*median(coef$Beta_2[!is.na(coef$Beta_2)])))
        #print(paste0("Estimate for M0 is ",round(M0,2)))
        
        C <- 2/(1-median(coef$Beta_1[!is.na(coef$Beta_1)])*M0)
        
        k <- 1/(M0*C*median(coef$Beta_0[!is.na(coef$Beta_0)]))
        
        data$y_hat <- 1/(M0*C*k) + ((C-2)/M0*C)*data$aw + ((1-C)*k/M0*C)*(data$aw)^2
        data$dif <- data$y -data$y_hat
        MSE <- sum((data$dif)^2)/nrow(data)#find the Mean Squared Error
        
        #print(paste0("Mean Squared Error(MSE) of the bet model is ",round(MSE,2)))
        
        rsq <- median(coef$R_square) * 100#To get the R squared
        
        #print(paste0("R squared value is ",round(rsq,2), "%"))
        output <- data.frame(C = C, 
                             M0 = M0,
                             k = k,
                             R_square=rsq,
                             MSE = MSE)
        return(output)
        
        
    }
    
    
    
    
    output$Model <- renderTable({
        if(is.null(input$file)){return()}
        if(input$Select == "Smith Model"){
            
            return(Smith())}
        
        else if(input$Select == "Halsey Model"){
            
            return(Halsey())}
        
        else if(input$Select == "Brunauer Emmett Teller (BET) Model"){
            
            return(Bet())}
        
        else if(input$Select == "Oswin Model"){
            
            return(oswin())}
        else if(input$Select == "Guggenheim Anderson de Boer (GAB) Model"){
            
            return(GAB())}
        
        else if(input$Select == "All"){
            df <- data.frame(Model = c("Smith Model","Halsey Model","Brunauer Emmett Teller (BET) Model","Oswin Model","Guggenheim Anderson de Boer (GAB) Model"),
                             MSE = c(Smith()$MSE,Halsey()$MSE,Bet()$MSE,oswin()$MSE,GAB()$MSE),
                             R_squared = c(Smith()$R_square,Halsey()$R_square,Bet()$R_square,oswin()$R_square,GAB()$R_square))
            return(df)}
        
        
    })
    
    
    
    
    
    
    output$Plot <- renderPlot({
        
        if(is.null(input$file)){return()}
        aw_new<-seq(0,0.999,length=100)
        
        
        if(input$Select == "Smith Model"){
            data <- read.csv(input$file$datapath, header=input$header, sep=input$sep)
            dd2 <- Smith()
            MW1 <- dd2$c1 - dd2$c2 * log(1-aw_new)  #fitted Mw
            df_fitted <- data.frame(aw_new,MW1)#Fitted values dataframe
        
            ggplot() + 
                geom_line(data=df_fitted, aes(x=aw_new, y=MW1), color='green') + 
                geom_point(data=data, aes(x=aw, y=Mw))+xlab("aw")+ylab("Mw")+ggtitle( "Fitted vs Actual curves")
        }
        
        else if(input$Select == "Halsey Model"){
            data <- read.csv(input$file$datapath, header=input$header, sep=input$sep)
            dd <- Halsey()
            MW1 <-(-dd$c/log(aw_new))^(1/dd$n)  #fitted Mw
            df_fitted <- data.frame(aw_new,MW1)#Fitted values dataframe
            
            ggplot() + 
                geom_line(data=df_fitted, aes(x=aw_new, y=MW1), color='red') + 
                geom_point(data=data, aes(x=aw, y=Mw))+xlab("aw")+ylab("Mw")+ggtitle( "Fitted vs Actual curves")
            
        }
        
        else if(input$Select == "Brunauer Emmett Teller (BET) Model"){
            data <- read.csv(input$file$datapath, header=input$header, sep=input$sep)
            dd1 <- Bet()
            MW1 <- dd1$M0*dd1$C*aw_new/((1-aw_new)*(1+aw_new*(dd1$C-1)))  #fitted Mw
            df_fitted<-data.frame(aw_new,MW1)#Fitted values dataframe
            
            ggplot() + 
                geom_line(data=df_fitted, aes(x=aw_new, y=MW1), color='forestgreen') + 
                geom_point(data=data, aes(x=aw, y=Mw))+xlab("aw")+ylab("Mw")+ggtitle( "Fitted vs Actual curves")
            
        }
        
        else if(input$Select == "Oswin Model"){
            
            data <- read.csv(input$file$datapath, header=input$header, sep=input$sep)
            dd3 <- oswin()
            
            MW1 <- dd3$C*(aw_new/(1-aw_new))^dd3$n  #fitted Mw
            df_fitted<-data.frame(aw_new,MW1)#Fitted values dataframe
            
            ggplot() + 
                geom_line(data=df_fitted, aes(x=aw_new, y=MW1), color='blue') + 
                geom_point(data=data, aes(x=aw, y=Mw))+xlab("aw")+ylab("Mw")+ggtitle( "Fitted vs Actual curves")
            
            
        }
        
        else if(input$Select == "All"){
            
            data <- read.csv(input$file$datapath, header=input$header, sep=input$sep)
            
            dd2 <- Smith()
            MW1 <- dd2$c1 - dd2$c2 * log(1-aw_new)  #fitted Mw
            df_fitted2 <- data.frame(aw_new,MW1)#Fitted values dataframe
            
            dd3 <- oswin()
            
            MW11 <- dd3$C*(aw_new/(1-aw_new))^dd3$n  #fitted Mw
            df_fitted3<-data.frame(aw_new,MW11)#Fitted values dataframe
            
            dd <- Halsey()
            
            MW111 <-(-dd$c/log(aw_new))^(1/dd$n)  #fitted Mw
            df_fitted <- data.frame(aw_new,MW111)#Fitted values dataframe
            
            dd1 <- Bet()
            MW1111 <- dd1$M0*dd1$C*aw_new/((1-aw_new)*(1+aw_new*(dd1$C-1)))  #fitted Mw
            df_fitted4 <-data.frame(aw_new,MW1111)#Fitted values dataframe
            
            
            dd4 <- GAB()
            MW12 <- (dd4$M0*dd4$C*dd4$k*aw_new)/((1-dd4$k*aw_new)*(1-dd4$k*aw_new + dd4$C*dd4$k*aw_new))
            df_fitted5 <-data.frame(aw_new,MW12)#Fitted values dataframe
            
            ggplot() + 
                geom_line(data=df_fitted2, aes(x=aw_new, y=MW1), color='green') +
                geom_line(data=df_fitted3, aes(x=aw_new, y=MW11), color='red') +
                geom_line(data=df_fitted, aes(x=aw_new, y=MW111), color='forestgreen') +
                geom_line(data=df_fitted4, aes(x=aw_new, y=MW1111), color='blue') +
                geom_line(data=df_fitted5, aes(x=aw_new, y=MW12), color='orange') +
                geom_point(data=data, aes(x=aw, y=Mw))+xlab("aw")+ylab("Mw")+ggtitle( "Fitted vs Actual curves")
            
            
        }
        
        else if(input$Select == "Guggenheim Anderson de Boer (GAB) Model"){
            
            data <- read.csv(input$file$datapath, header=input$header, sep=input$sep)
            dd4 <- GAB()
            MW1 <- (dd4$M0*dd4$C*dd4$k*aw_new)/((1-dd4$k*aw_new)*(1-dd4$k*aw_new + dd4$C*dd4$k*aw_new))
            df_fitted<-data.frame(aw_new,MW1)#Fitted values dataframe
            
            ggplot() + 
                geom_line(data=df_fitted, aes(x=aw_new, y=MW1), color='orange') + 
                geom_point(data=data, aes(x=aw, y=Mw))+xlab("aw")+ylab("Mw")+ggtitle( "Fitted vs Actual curves")
            
            
        }
        
        
        
        
    })
    
    
    output$tb <- renderUI({
        if(is.null(input$file)){return()}
        else
            tabsetPanel(
                tabPanel("Dataset", tableOutput("contents")),
                tabPanel("Models for sorption isotherms",
                         fluidRow(column(10, offset = 1,
                                         tableOutput('Model'))),
                         fluidRow(column(10, offset = 1,
                                         plotOutput('Plot')))))
        
    })
    #output$distPlot <- renderPlot({

        # generate bins based on input$bins from ui.R
        #x    <- faithful[, 2]
        #bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        #hist(x, breaks = bins, col = 'darkgray', border = 'white')

    #})

})
