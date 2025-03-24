library(shiny)
library(ggplot2)
library(faux)
library(dplyr)
library(ThemePark)
library(tidyverse)
library(QuantPsyc)
library(pROC)
library(margins)


ui <- navbarPage(
  "Interactive Data Analysis",
    # Page 1: Group Distributions and Metrics
  tabPanel(
    "Between subject design",
    titlePanel("Compare Two Groups with Distributions"),
    sidebarLayout(
      sidebarPanel(
        numericInput("mean1", "Group 1 Mean:", value = 5),
        numericInput("sd1", "Group 1 SD:", value = 2),
        numericInput("n1", "Group 1 Sample Size:", value = 30),
        numericInput("mean2", "Group 2 Mean:", value = 6.6),
        numericInput("sd2", "Group 2 SD:", value = 2),
        numericInput("n2", "Group 2 Sample Size:", value = 30),
        radioButtons("plot_type", "Select Plot Type:", 
                     choices = c("Density Plot" = "density", "Dotplot" = "dotplot"),
                     selected = "density")
        ),
      mainPanel(
        plotOutput("distribution_plot"),
        h3("Metrics"),
        verbatimTextOutput("metrics_output")
      )
    )
  ),
  tabPanel(
    "Within subject design",
    titlePanel("Compare Two Groups with Distributions"),
    sidebarLayout(
      sidebarPanel(
        numericInput("mean1_w", "Group 1 Mean:", value = 5),
        numericInput("sd1_w", "Group 1 SD:", value = 2),
        numericInput("mean2_w", "Group 2 Mean:", value = 6.6),
        numericInput("sd2_w", "Group 2 SD:", value = 2),
        numericInput("r", "Correlation", value = .5, min = -.99, max= .99, step=.01),
        numericInput("N", "Sample Size:", value = 30),
        numericInput("SESOI", "Smallest Effect size of Interest", value = .2),
        radioButtons("plot_type_w", "Select Plot Type:", 
                     choices = c("Density Plot" = "density_w", "Dotplot" = "dotplot_w"),
                     selected = "density_w")
      ),
      mainPanel(
        plotOutput("distribution_plot_w"),
        h3("Metrics_w"),
        verbatimTextOutput("metrics_output_w")
      )
    )
  ),
  # Tab 3: Correlation Analysis
  tabPanel(
    "Correlation Analysis",
    titlePanel("Compare Two Groups with Distributions"),
    sidebarLayout(
      sidebarPanel(
        numericInput("mean1_c", "Group 1 Mean:", value = 5),
        numericInput("sd1_c", "Group 1 SD:", value = 2),
        numericInput("mean2_c", "Group 2 Mean:", value = 6.6),
        numericInput("sd2_c", "Group 2 SD:", value = 2),
        numericInput("r_c", "Correlation", value = .5, min = -.99, max= .99, step=.01),
        numericInput("N_c", "Sample Size:", value = 30),
        radioButtons("line_type", "Select line type:", 
                     choices = c("Linear model" = "lm", "LOESS" = "LOESS"),
                     selected = "lm")
        ),
      mainPanel(
        plotOutput("Scatterplot"),
        h3("Metrics_c"),
        verbatimTextOutput("metrics_output_c")
             )
           )
  ),
  # Tab 4: Logistic Regression
  tabPanel("Logistic Regression Analysis",
           sidebarLayout(
             sidebarPanel(
               numericInput("probgroup0", "Probability group 1:", value = .2, max = 1, min = 0, step = .01),
               numericInput("probgroup1", "Probability group 2:", value = .5, max = 1, min = 0, step = .01),
               numericInput("N_lr1", "Sample Size group 1:", value = 100, min = 10),
               numericInput("N_lr2", "Sample Size group 2:", value = 100, min = 10),
               radioButtons("graph_type_log", "Select Plot Type:", 
                            choices = c("Dotplot" = "dotplot", "Bargraph" = "bargraph"),
                            selected = "dotplot")
             ),
             mainPanel(
               plotOutput("logistic_plot"),
               verbatimTextOutput("logistic_metrics")
             )
           )
  ),
  # Tab 5: Poisson Regression
  tabPanel("Poisson Regression Analysis",
           sidebarLayout(
             sidebarPanel(
               numericInput("countgroup1", "Expected event rate group 1:", value = 2, min = 0, step = 1),
               numericInput("countgroup2", "Expected event rate group 2:", value = 3, min = 0, step = 1),
               numericInput("N_pr1", "Sample Size group 1:", value = 100, min = 10),
               numericInput("N_pr2", "Sample Size group 2:", value = 100, min = 10),
               radioButtons("graph_type", "Select Plot Type:", 
                            choices = c("Dotplot" = "dotplot", "Histogram" = "histogram"),
                            selected = "dotplot")
             ),
             mainPanel(
               plotOutput("poisson_plot"),
               verbatimTextOutput("poisson_metrics")
             )
           )
  )
)



server <- function(input, output, session) {
  ### Page 1: Between subject design
  metrics <- reactive({
    mean1 <- input$mean1
    sd1 <- input$sd1
    n1 <- input$n1
    mean2 <- input$mean2
    sd2 <- input$sd2
    n2 <- input$n2
    
    pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
    cohend <- (mean2 - mean1) / pooled_sd
    auc <- pnorm(cohend / sqrt(2))
    success_rate_diff <- 2 * auc - 1
    nnt <- 1 / (pnorm(cohend + qnorm(.2)) - .2)
    
    list("Cohen's d" = cohend, AUC = auc, SRD = success_rate_diff, NNT = nnt, Pooled_SD = pooled_sd)
  })
  
  output$distribution_plot <- renderPlot({
    mean1 <- input$mean1
    sd1 <- input$sd1
    n1 <- input$n1
    mean2 <- input$mean2
    sd2 <- input$sd2
    n2 <- input$n2
    
    if (input$plot_type == "density") {
      x <- seq(min(mean1 - 4 * sd1, mean2 - 4 * sd2), 
               max(mean1 + 4 * sd1, mean2 + 4 * sd2), length.out = 100)
      y1 <- dnorm(x, mean1, sd1)
      y2 <- dnorm(x, mean2, sd2)
      
      plot(x, y1, type = "l", col = "blue", lwd = 2, ylim = c(0, max(y1, y2)),
           xlab = "Value", ylab = "Density", main = "Group Distributions")
      lines(x, y2, col = "red", lwd = 2)
      legend("topright", legend = c("Group 1", "Group 2"), 
             col = c("blue", "red"), lwd = 2)
    } else if (input$plot_type == "dotplot") {
      set.seed(2794)
      correct2 <- sim_design(
        between = list(group = c("g1", "g2")), # Independent variables
        n = c(n1,n2),                                       # Total number of participants
        mu = c(mean1, mean2),                        # Means for each group
        sd = c(sd1, sd2),                            # Standard deviations
        rep = 1,
        id = "participants",                         # Identifier for data points
        plot = FALSE
      )
      
      # Dynamic dot size adjustment based on the number of points
      dotsize <- 1 / log10(n1 + 10)
      
      # Calculate y-axis limits
      y_min <- min(correct2$y)
      y_max <- max(correct2$y)

            # Create the dotplot
      ggplot(correct2, aes(x = group, y = y, fill = group)) +
        geom_dotplot(
          binaxis = "y", 
          stackdir = "center", 
          dotsize = dotsize
        ) +
        scale_fill_manual(values = c("blue", "red")) +
        labs(
          x = "Group", 
          y = "Value", 
          title = "Dotplot of Group Distributions"
        ) +
        ylim(y_min, y_max) + # Adjust Y-axis limits
        theme_barbie() +
        theme(
          plot.title = element_text(hjust = 0.5), # Center align the title
          axis.title.x = element_text(size = 14), 
          axis.title.y = element_text(size = 14)
        )
    }
  })
  
  output$metrics_output <- renderText({
    m <- metrics()
    
    paste(
      "Summary of Metrics:",
      "",
      sprintf("1. Cohen's d: %.2f", m[["Cohen's d"]]),
      "   Interpretation: Cohen's d represents the standardized mean difference.",
      "",
      sprintf("2. AUC: %.2f", m[["AUC"]]),
      "   Interpretation: AUC (Area Under the Curve) reflects the probability that a randomly chosen individual from Group 1 has a higher value than a randomly chosen individual from Group 2. Values closer to 1 indicate strong discrimination.",
      "",
      sprintf("3. Success Rate Difference (SRD): %.2f", m[["SRD"]]),
      "   Interpretation: SRD represents the difference in success rates between groups. Values range from -1 to 1, with higher positive values indicating better success for Group 1.",
      "",
      sprintf("4. Number Needed to Treat (NNT): %.2f", m[["NNT"]]),
      "   Interpretation: NNT indicates how many individuals need to be treated in Group 1 to achieve one additional success compared to Group 2. Lower values indicate a more effective treatment.",
      "",
      sprintf("5. Pooled Standard Deviation: %.2f", m[["Pooled_SD"]]),
      "   Interpretation: The pooled standard deviation combines the variability of both groups and is used as the denominator in the calculation of Cohen's d.",
      sep = "\n"
    )
  })
  
  ### Page 2: Within subject design
  metrics_w <- reactive({
    # Input parameters
    mean1 <- input$mean1_w
    sd1 <- input$sd1_w
    mean2 <- input$mean2_w
    sd2 <- input$sd2_w
    r <- input$r
    n <- input$N
    SESOI <- input$SESOI

    # Simulate the data
    set.seed(2794)
    correct2 <- sim_design(
      within = list(group = c("t1", "t2")), # Independent variables
      n = n,                                       # Total number of participants
      mu = c(mean1, mean2),                        # Means for each group
      sd = c(sd1, sd2),                            # Standard deviations
      r = r,                                       # Correlation between groups
      rep = 1,
      id = "participants",                         # Identifier for data points
      plot = FALSE
    )
    
    pooled_sd_w <- sqrt(sd1^2 + sd2^2 - 2 * r *sd1 * sd2)
    cohend_w <- (mean2-mean1)/pooled_sd_w
    auc_w <- pnorm(cohend_w/sqrt(2))
    success_rate_diff_w <- 2*auc_w - 1
    nnt_w <- 1 / (pnorm(cohend_w + qnorm(.2))-.2)
    correct2$direction <- ifelse(
      (correct2$t2 - correct2$t1) > SESOI, "greater",
      ifelse(
        (correct2$t2 - correct2$t1) < -SESOI, "smaller",
        "equivalent"
      )
    )
    greater <- mean(correct2$direction == "greater")
    smaller <- mean(correct2$direction == "smaller")
    equal <- mean(correct2$direction == "equivalent")
    
    list("Cohen's dz" = cohend_w, "t2 >= t1" = greater, "t2 <= t1" = smaller, "t2 = t1" = equal, AUC = auc_w, SRD = success_rate_diff_w, NNT = nnt_w, Pooled_SD = pooled_sd_w)
  })
  
  output$distribution_plot_w <- renderPlot({
    # Input parameters
    mean1 <- input$mean1_w
    sd1 <- input$sd1_w
    mean2 <- input$mean2_w
    sd2 <- input$sd2_w
    r <- input$r
    n <- input$N
    
    # Simulate the data
    if (input$plot_type_w == "density_w") {
      x <- seq(min(mean1 - 4 * sd1, mean2 - 4 * sd2), 
               max(mean1 + 4 * sd1, mean2 + 4 * sd2), length.out = 100)
      y1 <- dnorm(x, mean1, sd1)
      y2 <- dnorm(x, mean2, sd2)
      
      plot(x, y1, type = "l", col = "blue", lwd = 2, ylim = c(0, max(y1, y2)),
           xlab = "Value", ylab = "Density", main = "Group Distributions")
      lines(x, y2, col = "red", lwd = 2)
      legend("topright", legend = c("Group 1", "Group 2"), 
             col = c("blue", "red"), lwd = 2)
    } else if (input$plot_type_w == "dotplot_w") {
      set.seed(2794)
      correct2 <- sim_design(
        within = list(group = c("t1", "t2")), # Independent variables
        n = n,                                       # Total number of participants
        mu = c(mean1, mean2),                        # Means for each group
        sd = c(sd1, sd2),                            # Standard deviations
        r = r,                                       # Correlation between groups
        rep = 1,
        id = "participants",                         # Identifier for data points
        plot = FALSE
      )
      
      # Pivot the data from wide to long format
      data <- correct2 %>%
        pivot_longer(
          cols = c(t1, t2), # Columns to pivot
          names_to = "Group",       # New column for variable names
          values_to = "value"       # New column for values
        )
      
      # Dynamic dot size adjustment based on the number of points
      dotsize <- 1 / log10(n + 10)
      
      # Calculate y-axis limits
      y_min <- min(data$value)
      y_max <- max(data$value)
      
      # Create the dotplot
      ggplot(data, aes(x = Group, y = value, fill = Group)) +
        geom_dotplot(
          binaxis = "y", 
          stackdir = "center", 
          dotsize = dotsize
        ) +
        scale_fill_manual(values = c("blue", "red")) +
        labs(
          x = "Group", 
          y = "Value", 
          title = "Dotplot of Group Distributions"
        ) +
        ylim(y_min, y_max) + # Adjust Y-axis limits
        theme_barbie() +
        theme(
          plot.title = element_text(hjust = 0.5), # Center align the title
          axis.title.x = element_text(size = 14), 
          axis.title.y = element_text(size = 14)
        )
    }
  })
  
  output$metrics_output_w <- renderText({
    m <- metrics_w()
    
    paste(
      "Summary of Metrics:",
      "",
      sprintf("1. Cohen's dz: %.2f", m[["Cohen's dz"]]),
      "   Interpretation: Cohen's dz is the standardized mean difference for paired samples.",
      "",
      sprintf("2. t2 >= t1 Proportion: %.2f", m[["t2 >= t1"]]),
      "   Interpretation: The proportion of cases where t2 is greater than or equal to t1, adjusted for the SESOI.",
      "",
      sprintf("3. t2 <= t1 Proportion: %.2f", m[["t2 <= t1"]]),
      "   Interpretation: The proportion of cases where t2 is less than or equal to t1, adjusted for the SESOI.",
      "",
      sprintf("4. t2 = t1 Proportion: %.2f", m[["t2 = t1"]]),
      "   Interpretation: The proportion of cases where t2 and t1 are equivalent within the SESOI range.",
      "",
      sprintf("5. AUC: %.2f", m[["AUC"]]),
      "   Interpretation: The AUC (Area Under the Curve) indicates the probability that a randomly chosen individual has a higher value in t2 compared to t1.",
      "",
      sprintf("6. Success Rate Difference (SRD): %.2f", m[["SRD"]]),
      "   Interpretation: SRD represents the success rate difference between conditions. \n   Values closer to ±1 indicate stronger differences.",
      "",
      sprintf("7. Number Needed to Treat (NNT): %.2f", m[["NNT"]]),
      "   Interpretation: NNT reflects how many participants need to be in condition t2 to achieve one additional success compared to condition t1.",
      "",
      sprintf("8. Pooled Standard Deviation: %.2f", m[["Pooled_SD"]]),
      "   Interpretation: The pooled standard deviation accounts for the correlation between the two conditions and measures overall variability.",
      sep = "\n"
    )
  })
  
  
  ### Tab 3: Correlation Analysis
  metrics_c <- reactive({
    # Input parameters
    mean1 <- input$mean1_c
    sd1 <- input$sd1_c
    mean2 <- input$mean2_c
    sd2 <- input$sd2_c
    r <- input$r_c
    n <- input$N_c

    # Simulate the data
    set.seed(2794)
    correct2 <- sim_design(
      within = list(group = c("X", "Y")), # Independent variables
      n = n,                                       # Total number of participants
      mu = c(mean1, mean2),                        # Means for each group
      sd = c(sd1, sd2),                            # Standard deviations
      r = r,                                       # Correlation between groups
      rep = 1,
      id = "participants",                         # Identifier for data points
      plot = FALSE
    )
    
    mod <- lm(Y~X, data=correct2)
    simplES <- mod$coefficients[2]
    cor <- lm.beta(mod)
    BESD <- cor*100/2+50
    NNT <- pi / (atan((abs(cor)/sqrt(1-(cor^2)))) + asin(cor))
    r2 <- summary(mod)$r.squared

    list("Unstandardized regression coefficient" = simplES, Correlation = cor, BESD = BESD, NNT=NNT, r2 = r2)
  })
  
  output$Scatterplot <- renderPlot({
    mean1 <- input$mean1_c
    sd1 <- input$sd1_c
    mean2 <- input$mean2_c
    sd2 <- input$sd2_c
    r <- input$r_c
    n <- input$N_c
    
    set.seed(2794)
    correct2 <- sim_design(
    within = list(group = c("X", "Y")), # Independent variables
    n = n,                                       # Total number of participants
    mu = c(mean1, mean2),                        # Means for each group
    sd = c(sd1, sd2),                            # Standard deviations
    r = r,                                       # Correlation between groups
    rep = 1,
    id = "participants",                         # Identifier for data points
    plot = FALSE
    )    
    if (input$line_type == "lm") {
    ggplot(correct2, aes(x = X, y = Y)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      labs(
        title = "Simulated Data with Custom Slope",
        x = "Predictor (X)",
        y = "Outcome (Y)"
      ) +
      theme_barbie()
      } else if(input$line_type == "LOESS") {
      ggplot(correct2, aes(x = X, y = Y)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "loess", se = TRUE, color = "red") +
        labs(
          title = "Simulated Data with Custom Slope",
          x = "Predictor (X)",
          y = "Outcome (Y)"
        ) +
        theme_barbie()}
  })
  
  output$metrics_output_c <- renderText({
    m <- metrics_c()
    paste(
      "Summary of Metrics:",
      "",
      sprintf("1. Unstandardized regression coefficient: %.2f", m[["Unstandardized regression coefficient"]]),
      "   Interpretation: The unstandardized regression coefficient indicates the expected change in the outcome (Y) for a one-unit increase in the predictor (X).",
      "",
      sprintf("2. Correlation: %.2f", m[["Correlation"]]),
      "   Interpretation: The correlation measures the strength and direction of the linear relationship between X and Y.",
      "",
      sprintf("3. Rsquared: %.2f", m[["r2"]]),
      "   Interpretation: How much variation in the dependent variable is explained by the regression model.",
      "",
      sprintf("4. BESD: %.2f", m[["BESD"]]),
      sprintf(
        "   Interpretation: Based on the predictor, %.2f%% of individuals in the treatment group succeed compared to %.2f%% in the control group.",
        m[["BESD"]], 100 - m[["BESD"]]
      ),
      "",
      sprintf("5. Number Needed to Treat (NNT): %.2f", m[["NNT"]]),
      "   Interpretation: NNT reflects how many participants need to be in condition t2 to achieve one additional success compared to condition t1.",
      "",
      sep = "\n"
    )
  })
  
  # Logistic Regression Data Generation
  logistic_data <- reactive({
    set.seed(2794)
    
    N1 <- input$N_lr1  # Total sample size
    N2 <- input$N_lr2  # Total sample size
    
    p_group0 <- input$probgroup0  # Probability of success in group 1
    p_group1 <- input$probgroup1  # Probability of success in group 2
    
    group <- c(rep(0, each = N1), rep(1, each = N2))  # Binary group variable
    
    # Generate binary outcomes based on group probability
    y <- rbinom(N1 + N2, size = 1, prob = ifelse(group == 0, p_group0, p_group1))
    
    # Create dataframe
    sim_data <- data.frame(y = y, group = group)
    
    # Fit logistic regression model
    logit_model <- glm(y ~ group, family = binomial, data = sim_data)
    
    # Get predicted probabilities
    sim_data$prob <- predict(logit_model, type = "response")  # Predicted probabilities
    
    return(sim_data)
  })
  
  # Plot Logistic Regression Predictions
  output$logistic_plot <- renderPlot({
    data <- logistic_data()

    if (input$graph_type_log == "dotplot") {
      ggplot(data, aes(x = as.factor(group), y = y, color = as.factor(y))) +  # Use group instead of x
      geom_jitter(alpha = 0.6, position = position_jitter(width = 0.1, height = 0.1)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 1), labels = c("Fail", "Success")) +
      labs(
        title = "Logistic Regression: Predicted Probabilities",
        x = "Group",
        y = "Predicted Probability",
        color = "Observed Outcome"
      ) +
      theme_barbie()
    } else if (input$graph_type_log == "bargraph") {
      ggplot(data, aes(x = as.factor(group), fill = as.factor(y))) +
        geom_bar(position = "fill") +  # Stacked bars, scaled to proportions
        scale_y_continuous(labels = scales::percent_format()) +  # Convert to percentages
        labs(
          title = "Logistic Regression: Proportion of Successes and Failures",
          x = "Group",
          y = "Proportion",
          fill = "Outcome"
        ) +
        theme_minimal()}
  })
  
  # Logistic Regression Performance Metrics
  output$logistic_metrics <- renderText({
    data <- logistic_data()
    
    # Handle NULL data scenario
    if (is.null(data)) return("No data available. Adjust N or probabilities.")
    
    # Fit logistic regression model
    logit_model <- glm(y ~ group, family = binomial, data = data)
    
    # Compute Odds Ratio (OR)
    OR <- exp(coef(logit_model)["group"])
    
    # Compute Average Marginal Effect (AME)
    AME <- summary(margins(logit_model))
    
    # Compute Risk Ratio (RR)
    p1 <- mean(data$y[data$group == 0])  # Mean success rate in Group 1
    p2 <- mean(data$y[data$group == 1])  # Mean success rate in Group 2
    RR <- p1 / p2  # Risk Ratio formula
    
    # Risk difference
    risk_difference <- p1 - p2  

    # AUC calculation
    auc <- roc(data$y, data$group)$auc  
    
    # Survival ratio
    survival_ratio <- (1 - p1) / (1 - p2) 
    
    # Compute NNT / NNH
    nnt_nnh <- ifelse(abs(risk_difference) > 0, 1 / abs(risk_difference), Inf)
    nnt_nnh_label <- ifelse(risk_difference > 0, "Number Needed to Treat (NNT)", "Number Needed to Harm (NNH)")
    
    # Interpretations
    paste(
      "Logistic Regression Effect Sizes:",
      "",
      sprintf("1. Odds Ratio (OR): %.2f", OR),
      "   Interpretation: The odds ratio represents how much more (or less) likely success is in Group 1 compared to Group 2.",
      "",
      sprintf("2. Risk Ratio (RR): %.2f", RR),
      "   Interpretation: The risk ratio indicates how many times more (or less) likely success is in Group 1 compared to Group 2.",
      "",
      sprintf("3. Average Marginal Effect (AME): %.3f", AME$AME),
      sprintf("   Interpretation: On average, being in Group 1 changes the probability of success by %.2f percentage points.", abs(AME$AME) * 100),
      "",
      sprintf("4. Survival Ratio (SR): %.2f", survival_ratio),
      "   Interpretation: If SR > 1, Group 1 has a higher survival rate; if SR < 1, Group 1 has a higher failure rate compared to Group 2.",
      "",
      sprintf("5. Area Under Curve (AUC): %.2f", auc),
      "   Interpretation: AUC (Area Under the Curve) reflects the probability that a randomly chosen individual from Group 1 has a higher value than a randomly chosen individual from Group 2. Values closer to 1 indicate strong discrimination.",
      "",
      sprintf("6. %s: %.2f", nnt_nnh_label, nnt_nnh),
      "   Interpretation: NNT indicates how many individuals need to be treated in Group 1 to achieve one additional success compared to Group 2. Lower values indicate a more effective treatment.",
      "",
      sep = "\n"
    )
  })
  
  # Poisson Regression Data Generation
  poisson_data <- reactive({
    set.seed(2794)
    
    N1 <- input$N_pr1  # Sample size group 1
    N2 <- input$N_pr2  # Sample size group 1
    
    
    lambda_group1 <- input$countgroup1  # Probability of success in group 1
    lambda_group2 <- input$countgroup2  # Probability of success in group 2
    
    group <- c(rep(0, each = N1),rep(1, each=N2))  # Binary group variable
    
    # Step 3: Simulate Poisson-distributed counts
    counts <- rpois(N1 + N2, lambda = ifelse(group == 0, lambda_group1, lambda_group2))
    
    # Step 4: Create data frame
    sim_data <- data.frame(counts = counts, group = factor(group))
    
    # Step 5: Fit Poisson regression model
    poisson_model <- glm(counts ~ group, family = poisson, data = sim_data)
    
    return(sim_data)
  })
  
  # Plot Poisson Regression Predictions
  output$poisson_plot <- renderPlot({
    data <- poisson_data()
    if (input$graph_type == "dotplot"){
    ggplot(data, aes(x = group, y = counts, fill = group)) +
      geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.25, binwidth = 1, alpha = 0.6) +
      scale_x_discrete(labels = c("Group 0", "Group 1")) +
      labs(title = "Poisson Regression: Dot Plot of Observed Counts by Group",
           x = "Group", 
           y = "Event Counts") +
      theme_barbie() +
      theme(legend.position = "none")
    } else if (input$graph_type == "histogram"){
      ggplot(data, aes(counts, fill = group)) +
        geom_histogram(binwidth=.5, position="dodge")+
        labs(title = "Poisson Regression: Histogram Counts by Group",
             y = "Frequency", 
             x = "Event counts") +
        scale_y_continuous(expand = expansion(mult = c(0, 1)))+
        theme_barbie()
      }
  })
  
  # Logistic Regression Performance Metrics
  output$poisson_metrics <- renderText({
    data <- poisson_data()

    if (is.null(data)) return("No data available. Adjust sample size or lambda values.")
    
    # Poisson regression model
    poisson_model <- glm(counts ~ group, family = poisson, data = data)
    
    # Extract model coefficient
    Incidence_rate_ratio <- exp(coef(poisson_model)["group1"])  # Rate Ratio (RR)
    
    # Average marginal effect
    AME <- summary(margins(poisson_model))
    
    # Output the results
    paste(
      "Poisson Regression Effect Sizes:",
      "",
      sprintf("1. Incidence_rate_ratio: %.2f", Incidence_rate_ratio),
      "   Interpretation: The rate of events in Group 1 is RR times that of Group 0.",
      "",
      sprintf("2. Average marginal effect: %.2f", AME$AME),
      "   Interpretation: The rate of events in Group 1 is RR times that of Group 0.",
      "",
      

      sep = "\n"
    )
  })
}

shinyApp(ui, server)
