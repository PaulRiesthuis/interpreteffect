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
  "Interpret effect sizes",
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
        numericInput("CER", "Control group event rate", value = 20, min = 1, max = 100, step = 1),
        radioButtons("plot_type", "Select Plot Type:",
                     choices = c("Density Plot" = "density", "Dotplot" = "dotplot"),
                     selected = "density")
      ),
      mainPanel(
        plotOutput("distribution_plot"),
        verbatimTextOutput("metrics_output")
      )
    )
  ),
  tabPanel(
    "Within subject design",
    titlePanel("Compare Two Groups with Distributions"),
    sidebarLayout(
      sidebarPanel(
        numericInput("mean1_w", "Group 1 (Timepoint 1) Mean:", value = 5),
        numericInput("sd1_w", "Group 1 (Timepoint 1) SD:", value = 2),
        numericInput("mean2_w", "Group 2 (Timepoint 2) Mean:", value = 6.6),
        numericInput("sd2_w", "Group 2 (Timepoint 2) SD:", value = 2),
        numericInput("r", "Correlation", value = .5, min = -.99, max= .99, step=.01),
        numericInput("N", "Sample Size:", value = 30),
        numericInput("SESOI", "Region of equivalent scores between T1 and T2", value = 0),
        radioButtons("plot_type_w", "Select Plot Type:",
                     choices = c("Density Plot" = "density_w", "Dotplot" = "dotplot_w"),
                     selected = "density_w")
      ),
      mainPanel(
        plotOutput("distribution_plot_w"),
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
    CER <- input$CER
    mean_dif <- mean2 - mean1
    pooled_sd <- sqrt((sd1^2 + sd2^2) / 2)
    cohend <- mean_dif / pooled_sd
    auc <- pnorm(cohend / sqrt(2))
    nnt <- 1 / (pnorm(cohend + qnorm(CER/100)) - CER/100)
    nnt_nnh_label <- ifelse(nnt > 0, "Number Needed to Treat (NNT)", "Number Needed to Harm (NNH)")


    list(mean_dif = mean_dif, "Cohen's d" = cohend, AUC = auc, NNT = abs(nnt), nnt_nnh_label = nnt_nnh_label)
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
           xlab = "Value", ylab = "Density")
      lines(x, y2, col = "red", lwd = 2)
      legend("topright", legend = c("Group 1", "Group 2"),
             col = c("blue", "red"), lwd = 2)
    } else if (input$plot_type == "dotplot") {
      set.seed(2794)
      correct2 <- sim_design(
        between = list(group = c("Group 1", "Group 2")), # Independent variables
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
        ) +
        ylim(y_min, y_max) + # Adjust Y-axis limits
        theme_classic() +
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
      "Between subject effect sizes:",
      "",
      sprintf("1. Raw mean difference: %.2f", m[["mean_dif"]]),
      sprintf("   Interpretation: Participants in group 2 have, on average, %.2f points more than participants in the group 1", m[["mean_dif"]]),
      "",
      sprintf("2. Cohen's d: %.2f", m[["Cohen's d"]]),
      sprintf("   Interpretation: Participants in group 2 have, on average, %.2f standard deviations more than participants in the group 1", m[["Cohen's d"]]),
      "",
      sprintf("3. Probability of superiority (AUC): %.2f", m[["AUC"]]),
      "   Interpretation: The probability of superiority (Area Under the Curve) reflects the probability that a randomly chosen individual from Group 2 has a higher value than a randomly chosen individual from Group 1. Values closer to 1 indicate strong discrimination.",
      "",
      sprintf("4. %s: %.2f", m[["nnt_nnh_label"]], m[["NNT"]]),
      interpretation <- if (m[["nnt_nnh_label"]] == "Number Needed to Treat (NNT)") {
        sprintf("   Interpretation: The %s of %.2f means that on average, %.2f individuals need to be treated to have one additional success in Group 2 compared to Group 1.", m[["nnt_nnh_label"]], m[["NNT"]], m[["NNT"]])
      } else {
        sprintf("   Interpretation: The %s of %.2f means that on average, %.2f individuals need to be exposed to the condition in Group 2 for one additional adverse outcome compared to Group 1.", m[["nnt_nnh_label"]], m[["NNT"]], m[["NNT"]])
      },
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
    mean_dif <- mean2-mean1
    pooled_sd_w <- sqrt(sd1^2 + sd2^2 - 2 * r *sd1 * sd2)
    cohend_w <- mean_dif/pooled_sd_w
    auc_w <- pnorm(cohend_w/sqrt(2))
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

    list(mean_dif = mean_dif, "Cohen's dz" = cohend_w, "t2 >= t1" = greater, "t2 <= t1" = smaller, "t2 = t1" = equal, AUC = auc_w)
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
           xlab = "Value", ylab = "Density")
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
        ) +
        ylim(y_min, y_max) + # Adjust Y-axis limits
        theme_classic() +
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
      "Within subject effect sizes:",
      "",
      sprintf("1. Raw mean difference: %.2f", m[["mean_dif"]]),
      sprintf("   Interpretation: Participants in group 2 have, on average, %.2f points more than participants in the group 1", m[["mean_dif"]]),
      "",
      sprintf("2. Cohen's dz: %.2f", m[["Cohen's dz"]]),
      sprintf("   Interpretation: Participants in group 2 have, on average, %.2f standard deviations more than participants in the group 1", m[["Cohen's dz"]]),
      "",
      sprintf("3. Proportion of scores Group 2 (Timpoint 2) > Group 1 (Timpoint 1): %.2f", m[["t2 >= t1"]]),
      "   Interpretation: The proportion of participants that have scores at timepoint 2 (group 2) greater than  timepoint 1 (group 1).",
      "",
      sprintf("4. Proportion of scores Group 2 (Timpoint 2) < Group 1 (Timpoint 1): %.2f", m[["t2 <= t1"]]),
      "   Interpretation: The proportion of participants that have scores at timepoint 2 (group 2) less than  timepoint 1 (group 1).",
      "",
      sprintf("5. Proportion of scores Group 2 (Timpoint 2) = Group 1 (Timpoint 1): %.2f", m[["t2 = t1"]]),
      "   Interpretation: The proportion of participants that have equivalent scores on timepoint 1 (group 1) and 2 (group 2) based on the region of equivalent scores",
      "",
      sprintf("6. Probability of superiority (AUC): %.2f", m[["AUC"]]),
      sprintf("   Interpretation: The probability of superiority (Area Under the Curve) of %.2f indicates the probability that a randomly chosen individual from t2 has a higher value than a randomly chosen individual from t1.", m[["AUC"]]),
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
    CLES <- asin(sqrt(r))/pi + .5
    NNT <- pi / (atan((abs(cor)/sqrt(1-(cor^2)))) + asin(cor))
    r2 <- summary(mod)$r.squared

    list("Unstandardized regression coefficient" = simplES, Correlation = cor, CLES = CLES, NNT=NNT, r2 = r2)
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
          x = "Predictor (X)",
          y = "Outcome (Y)"
        ) +
        theme_classic()
    } else if(input$line_type == "LOESS") {
      ggplot(correct2, aes(x = X, y = Y)) +
        geom_point(alpha = 0.6) +
        geom_smooth(method = "loess", se = TRUE, color = "red") +
        labs(
          x = "Predictor (X)",
          y = "Outcome (Y)"
        ) +
        theme_classic()}
  })

  output$metrics_output_c <- renderText({
    paste(
      "Correlation effect sizes:",
      "",
      sprintf("1. Unstandardized Regression Coefficient: %.2f", metrics_c()[["Unstandardized regression coefficient"]]),
      sprintf("   Interpretation: A one-unit increase in the predictor (X) is associated with an expected change of %.2f units in the outcome (Y), holding all other variables constant.", metrics_c()[["Unstandardized regression coefficient"]]),
      "",
      sprintf("2. Correlation (of simulated data; might not equal the given correlation): %.2f", metrics_c()[["Correlation"]]),
      sprintf("   Interpretation: The correlation of %.2f indicates the strength and direction of the linear relationship between X and Y, with values closer to 1 or -1 indicating a stronger relationship.", metrics_c()[["Correlation"]]),
      "",
      sprintf("3. R-squared: %.2f", metrics_c()[["r2"]]),
      sprintf("   Interpretation: R-squared of %.2f represents the proportion of variance in the dependent variable that is explained by the regression model.", metrics_c()[["r2"]]),
      "",
      sprintf("4. Dunlap CLES: %.2f", metrics_c()[["CLES"]]),
      sprintf("   Interpretation: Based on the Dunlap CLES, there is a %.2f%% chance that in a pair of randomly selected people, the participant with the higher score on one variable will also have a higher score on the other variable.", metrics_c()[["CLES"]] * 100),
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
      ggplot(data, aes(x = as.factor(group), y = y, color = as.factor(y))) +
        geom_jitter(alpha = 0.6, position = position_jitter(width = 0.1, height = 0.1)) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 1), labels = c("Fail", "Success")) +
        scale_x_discrete(labels = c("1" = "Group 1", "2" = "Group 2")) +
        labs(
          x = "Group",
          y = "Predicted Probability",
          color = "Observed Outcome"
        ) +
        theme_classic()
    } else if (input$graph_type_log == "bargraph") {
      ggplot(data, aes(x = as.factor(group), fill = as.factor(y))) +
        geom_bar(position = "fill") +  # Stacked bars, scaled to proportions
        scale_y_continuous(labels = scales::percent_format()) +  # Convert to percentages
        labs(
          x = "Group",
          y = "Proportion",
          fill = "Outcome"
        ) +
        theme_classic()}
  })

  # Logistic Regression Performance Metrics
  output$logistic_metrics <- renderText({
    data <- logistic_data()

    # Handle NULL data scenario
    if (is.null(data)) return("No data available. Adjust N or probabilities.")

    # Fit logistic regression model
    logit_model <- glm(y ~ group, family = binomial, data = data)

    # Compute Odds Ratio (OR)
    OR <- (input$probgroup1/(1-input$probgroup1))/(input$probgroup0/(1-input$probgroup0))
    OR_label <- ifelse(OR > 1, "higher", "lower")

    # Compute Average Marginal Effect (AME)
    AME <- summary(margins(logit_model))
    AME_label <- ifelse(AME$AME > 0, "increases", "decreases")

    # Compute Risk Ratio (RR)
    RR <- input$probgroup1 / input$probgroup0  # Risk Ratio formula
    RR_label <- ifelse(RR > 1, "higher", "lower")

    # Risk difference
    risk_difference <- input$probgroup1 - input$probgroup0

    # AUC calculation
    AUC <- roc(data$y, data$group, direction = "<")$auc

    # Compute NNT / NNH
    nnt_nnh <- ifelse(abs(risk_difference) > 0, 1 / abs(risk_difference), Inf)
    nnt_nnh_label <- ifelse(risk_difference > 0, "Number Needed to Treat (NNT)", "Number Needed to Harm (NNH)")

    # Interpretations
    paste(
      "Logistic Regression Effect Sizes:",
      "",
      sprintf("1. Odds Ratio (OR): %.2f", OR),
      sprintf("   Interpretation: The odds ratio of %.2f means that the odds of success are %.2f times %s in Group 2 compared to Group 1.", OR, OR, OR_label),
      "",
      sprintf("2. Risk Ratio (RR): %.2f", RR),
      sprintf("   Interpretation: The risk ratio of %.2f indicates that the probability of success is %.2f %s in Group 2 compared to Group 1.", RR, RR, RR_label),
      "",
      sprintf("3. Average Marginal Effect (AME): %.3f", AME$AME),
      sprintf("   Interpretation: On average, being in Group 2 %s the probability of success by %.2f percentage points compared to Group 1.", AME_label, abs(AME$AME) * 100),
      "",
      sprintf("5. Area Under Curve (AUC): %.2f", AUC),
      sprintf("   Interpretation: An AUC of %.2f reflects the probability that a randomly chosen individual from Group 2 has a higher value than a randomly chosen individual from Group 1. Values closer to 1 indicate stronger discrimination.", AUC),
      "",
      sprintf("6. %s: %.2f", nnt_nnh_label, nnt_nnh),
      interpretation <- if (nnt_nnh_label == "Number Needed to Treat (NNT)") {
        sprintf("   Interpretation: The %s of %.2f means that on average, %.2f individuals need to be treated to have one additional success in Group 2 compared to Group 1.", nnt_nnh_label, nnt_nnh, nnt_nnh)
      } else {
        sprintf("   Interpretation: The %s of %.2f means that on average, %.2f individuals need to be exposed to the condition in Group 2 for one additional adverse outcome compared to Group 1.", nnt_nnh_label, nnt_nnh, nnt_nnh)
      },
      sep = "\n"
    )
  })

  # Poisson Regression Data Generation
  poisson_data <- reactive({
    set.seed(2794)

    N1 <- input$N_pr1  # Sample size group 1
    N2 <- input$N_pr2  # Sample size group 1


    lambda_group1 <- input$countgroup1
    lambda_group2 <- input$countgroup2

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
        scale_x_discrete(labels = c("Group 1", "Group 2")) +
        labs(x = "Group",
             y = "Event Counts") +
        theme_classic() +
        theme(legend.position = "none")
    } else if (input$graph_type == "histogram"){
      ggplot(data, aes(counts, fill = group)) +
        geom_histogram(binwidth=.5, position="dodge")+
        labs(y = "Frequency",
             x = "Event counts",
             fill = NULL) +
        scale_y_continuous(expand = expansion(mult = c(0, 1)))+
        scale_fill_discrete(labels = c("0" = "Group 1", "1" = "Group 2")) +
        theme_classic()
    }
  })

  # Logistic Regression Performance Metrics
  output$poisson_metrics <- renderText({
    data <- poisson_data()

    if (is.null(data)) return("No data available. Adjust sample size or lambda values.")

    # Poisson regression model
    poisson_model <- glm(counts ~ group, family = poisson, data = data)

    # Extract model coefficient
    Incidence_rate_ratio <- (input$countgroup2/input$N_pr2)/(input$countgroup1/input$N_pr1)

    # Average marginal effect
    AME <- summary(margins(poisson_model))
    AME_label <- ifelse(AME$AME > 0, "increase", "decrease")

    # Output the results
    paste(
      "Poisson Regression Effect Sizes:",
      "",
      sprintf("1. Incidence Rate Ratio: %.2f", Incidence_rate_ratio),
      sprintf("   Interpretation: The rate of events in Group 2 is %.2f times that of Group 1.", Incidence_rate_ratio),
      "",
      sprintf("2. Average Marginal Effect: %.2f", AME$AME),
      sprintf("   Interpretation: On average, being in Group 2 leads to an %s of %.2f events compared to Group 1.", AME_label, AME$AME),
      "",
      sep = "\n"
    )
  })
}

shinyApp(ui, server)
