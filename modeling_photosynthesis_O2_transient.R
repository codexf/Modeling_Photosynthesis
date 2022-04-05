# R script to analyze dynamic photosynthesis during oxygen transient. 
# The first half of the script analyzes the low oxygen (2% O2) to high oxygen (40% O2) transition.
# The second half of the script analyzes the high oxygen (40% O2) to low oxygen (2% O2) transition.
# Input csv file should have a column for time and columns for net CO2 assimilation rate,
# with replicates in separate columns.

# Author: XF
# V 1.0

# Libraries required
library(tidyverse)
library(scales)

# plot theme
plot_theme <- {
  theme(
    panel.background = element_rect(
      fill = NA,
      color = "black",
      size = 1.2
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(size = 0.5, colour = "black"),
    axis.text.x = element_text(
      angle = 0,
      size = 12,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      angle = 0,
      size = 12,
      hjust = 1,
      color = "black"
    ),
    axis.title = element_text(size = 12),
    legend.key = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    axis.ticks.length.x = unit(.08, "inch"),
    axis.ticks.length.y = unit(.08, "inch")
  )
}

# define output path
out_path <- 'output_O2_transient'
dir.create(out_path)


###############################################################################
##    Analysis for low oxygen (2% O2) to high oxygen (40% O2) transition     ##
###############################################################################


# read data for low (2% O2) to high oxygen (40% O2) transition (2-40)
data_2_40 <-
  read_csv("data_O2_transient/low_to_high_O2_transition_2-40.csv") %>%
  filter(Time < 240) # select time up to 240s
rep_name <- c(paste0('rep', (1:(ncol(
  data_2_40
) - 1))))
names(data_2_40) <- c('Time', rep_name)

# function to organize data for low (2% O2) to high oxygen (40% O2) transition (2-40)
makeDf_2_40 <- function(data, rep_no) {
  # compute the mean A before O2 transition
  mean_pre <- data %>%
    filter(Time < 0) %>%
    summarise_at(vars(-Time), ~ mean(., na.rm = TRUE))
  
  # compute the time of A apex
  time_maxA <- data %>% gather(key = Rep, value = A,-Time) %>%
    group_by(Rep) %>%
    slice(which.max(A)) %>%
    select(Rep, Time) %>%
    ungroup() %>%
    column_to_rownames("Rep") %>%
    t %>%
    as.data.frame()
  
  df <- data %>% select('Time', rep_no) %>%
    filter(Time > as.numeric(time_maxA[rep_no])) %>%  #filter points after the apex
    filter((!!as.symbol(rep_no))  <= as.numeric(mean_pre[rep_no])) %>% # filter points <= the mean A before O2 transition
    select(rep_no) %>%
    rename(y = rep_no) %>%
    mutate(t = seq(0, 2 * nrow(.) - 1, 2))  #add time_new <- seq(0, 2* nrow(rep1)-1, 2)
  return (df)
}


# Function to fit exponential model and calculated area for low (2% O2) to high oxygen (40% O2) transition (2-40)
getArea_2_40 <- function(df) {
  # Select an approximate $\theta$, since theta must be lower than min(y), and greater than zero
  theta.0 <- min(df$y) * 0.5
  
  # Estimate the rest parameters using a linear model
  model.0 <- lm(log(y - theta.0) ~ t, data = df)
  alpha.0 <- exp(coef(model.0)[1])
  beta.0 <- coef(model.0)[2]
  
  # Starting parameters
  start <- list(alpha = alpha.0,
                beta = beta.0,
                theta = theta.0)
  
  # Fit the model (with estimated starting parameters)
  model <-
    nls(y ~ alpha * exp(beta * t) + theta ,
        data = df,
        start = start)
  
  # integral of function
  ## define the integrated function
  alpha <- coef(model)['alpha']
  beta <- coef(model)['beta']
  theta <- coef(model)['theta']
  integrand <- function(x) {
    alpha * exp(beta * x) + theta
  }
  # integrate the function from 0 to 150
  upper_int <- 150
  intval <- integrate(integrand, lower = 0, upper = upper_int)
  
  # area between
  area <- intval$value - upper_int * theta
  names(area) <- 'area'
  
  output <- c(area, alpha, beta, theta)
  return (output)
}

# function to plot shaded curves for low (2% O2) to high oxygen (40% O2) transition (2-40)
plotShaded_2_40 <- function(df) {
  area_coef <- getArea_2_40(df)
  area <- area_coef['area']
  alpha <- area_coef['alpha']
  beta <- area_coef['beta']
  theta <- area_coef['theta']
  expFunc <- function(x) {
    alpha * exp(beta * x) + theta
  }
  df$fit <-
    expFunc(df$t)
  df$theta <- theta
  
  # shading using geom_ribbon
  plot <- ggplot(df, aes(x = t, y = y)) +
    geom_line(aes(x = t, y = fit), size = 1, color = 'black') +
    geom_line(aes(x = 0, y = fit), size = 1, color = 'black') +
    geom_line(aes(x = t, y = theta), size = 1, color = 'black') +
    geom_ribbon(aes(ymin = theta, ymax = fit),
                fill = "gray",
                alpha = .8) +
    geom_point(size = 1) +
    theme_bw() +
    ggtitle(paste0(
      "area =",
      round(area, 1),
      "     y = ",
      round(alpha, 3),
      "* exp(",
      round(beta, 3),
      "*x) +" ,
      round(theta, 3)
    )) +
    scale_x_continuous(name = 'Time after switching oxygen (s)',
                       #expand = c(0, 0),
                       limits = c(0, 155)) +
    scale_y_continuous(
      name = expression("A (\u03BCmol CO"[2] ~ m ^ -2 ~ s ^ -1 ~ ')'),
      #expand = c(0, 0),
      #limits = c(2, 7),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    theme(legend.position = "none") +
    plot_theme
  
  return (plot)
}

# combine traces from all reps for low (2% O2) to high oxygen (40% O2) transition (2-40)
datalist = list()
for (i in 1:10) {
  tem <- makeDf_2_40(data_2_40, rep_name[i])
  tem['id'] <- rep_name[i]   # Append new row
  datalist[[i]] <- tem
}
all_traces_2_40 = do.call(rbind, datalist)


# plot all traces with fitted lines for low (2% O2) to high oxygen (40% O2) transition (2-40)
# sort id by A value
id_order <- all_traces_2_40 %>%
  filter(t == 100) %>%
  arrange(y) %>%
  select(id) %>%
  pull()

# add the order to a new column
all_traces_2_40 <- all_traces_2_40 %>%
  rowwise() %>%
  mutate(ordered = which(id_order == id))


plot_traces <- all_traces_2_40 %>%
  ggplot(aes(
    x = t,
    y = y,
    group = ordered,
    color = ordered
  )) +
  geom_point(size = .5) +
  # SSasymp: Self-Starting Nls Asymptotic Regression Model
  stat_smooth(
    method = "nls",
    formula = y ~ SSasymp(x, Asym, R0, lrc),
    se = FALSE,
    size = 0.5
  ) +
  theme_bw() +
  scale_x_continuous(name = 'Time after switching oxygen (s)',
                     #expand = c(0, 0),
                     limits = c(0, 205)) +
  scale_y_continuous(
    name = expression("A (\u03BCmol CO"[2] ~ m ^ -2 ~ s ^ -1 ~ ')'),
    #expand = c(0, 0),
    limits = c(2.5, 7),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  #scale_colour_gradient(low = "grey", high = "black") +
  theme(legend.position = "none") +
  plot_theme


ggsave(
  paste0(out_path, "/figure_combined_traces_2_40.jpeg"),
  plot_traces,
  width = 6,
  height = 5,
  dpi = 600,
  units = "in"
)

# combine integrated areas from all reps for low (2% O2) to high oxygen (40% O2) transition (2-40)
datalist = list()

for (i in 1:10) {
  tem <-
    makeDf_2_40(data_2_40, rep_name[i]) %>% getArea_2_40()  # Create new row
  tem['id'] <- rep_name[i] # Append new row
  datalist[[i]] <- tem
}

all_area_2_40 = do.call(rbind, datalist)
write.csv(all_area_2_40, paste0(out_path, "/fitting_parameters_2_40.csv"))


# plot each shaded curve for low (2% O2) to high oxygen (40% O2) transition (2-40)
for (i in 1:10) {
  plot <- makeDf_2_40(data_2_40, rep_name[i]) %>% plotShaded_2_40()
  ggsave(
    paste0(out_path, "/figure_2_40_", rep_name[i], ".jpeg"),
    plot,
    width = 6,
    height = 5,
    dpi = 600,
    units = "in"
  )
}

###############################################################################
##    Analysis for high oxygen (40% O2) to low oxygen (2% O2) transition     ##
###############################################################################

#---------------------------------------------------------------------
# read data for high (40% O2) to low oxygen (2% O2) transition (40-2)
data_40_2 <- read_csv("data_O2_transient/high_to_low_O2_transition_40-2.csv") %>%
  filter(Time < 240) # select time up to 240s
rep_name <- c(paste0('rep', (1:(ncol(
  data_40_2
) - 1))))
names(data_40_2) <- c('Time', rep_name)


# function to organize data for high (40% O2) to low oxygen (2% O2) transition (40-2)
makeDf_40_2 <- function(data, rep_no) {
  # compute the mean A before O2 transition
  mean_pre <- data %>%
    filter(Time < 0) %>%
    summarise_at(vars(-Time), ~ mean(., na.rm = TRUE))
  
  # compute the time of A apex
  time_minA <- data %>% gather(key = Rep, value = A, -Time) %>%
    group_by(Rep) %>%
    slice(which.min(A)) %>%
    select(Rep, Time) %>%
    ungroup() %>%
    column_to_rownames("Rep") %>%
    t %>%
    as.data.frame()
  
  df <- data %>% select('Time', rep_no) %>%
    filter(Time > as.numeric(time_minA[rep_no])) %>%  #filter points after the apex
    filter((!!as.symbol(rep_no))  >= as.numeric(mean_pre[rep_no])) %>% # filter points <= the mean A before O2 transition
    select(rep_no) %>%
    rename(y = rep_no) %>%
    mutate(t = seq(0, 2 * nrow(.) - 1, 2))  #add time_new <- seq(0, 2* nrow(rep1)-1, 2)
  return (df)
}

# Function to fit negative exponential model and calculated area for high (40% O2) to low oxygen (2% O2) transition (40-2)
getArea_40_2 <- function(df) {
  # Exponential Model Fitting

  # Prepare a good initial state
  theta.0 <- max(df$y) * 1.1
  
  # Estimate the rest parameters using a linear model
  model.0 <- lm(log(-y + theta.0) ~ t, data = df)
  alpha.0 <- -exp(coef(model.0)[1])
  beta.0 <- coef(model.0)[2]
  
  # starting parameters
  start <- list(alpha = alpha.0,
                beta = beta.0,
                theta = theta.0)
  
  # Fit the model
  model <-
    nls(y ~ alpha * exp(beta * t) + theta ,
        data = df,
        start = start)
  
  
  # integral of function
  ## define the integrated function
  alpha <- coef(model)['alpha']
  beta <- coef(model)['beta']
  theta <- coef(model)['theta']
  integrand <- function(x) {
    alpha * exp(beta * x) + theta
  }
  # integrate the function from 0 to 150
  upper_int <- 150
  intval <- integrate(integrand, lower = 0, upper = upper_int)
  
  # area between
  area <- upper_int * theta - intval$value
  names(area) <- 'area'
  
  output <- c(area, alpha, beta, theta)
  return (output)
}

# function to plot for high (40% O2) to low oxygen (2% O2) transition (40-2)
plotShaded_40_2 <- function(df) {
  area_coef <- getArea_40_2(df)
  area <- area_coef['area']
  alpha <- area_coef['alpha']
  beta <- area_coef['beta']
  theta <- area_coef['theta']
  expFunc <- function(x) {
    alpha * exp(beta * x) + theta
  }
  df$fit <-
    expFunc(df$t)   # df$fit <- predict(model, list(x = df$t))
  df$theta <- theta # y = 3.096 * exp(-0.056 * x) + 3.328
  
  # shading using geom_ribbon
  plot <- ggplot(df, aes(x = t, y = y)) +
    geom_line(aes(x = t, y = fit), size = 1, color = 'black') +
    geom_line(aes(x = 0, y = fit), size = 1, color = 'black') +
    geom_line(aes(x = t, y = theta), size = 1, color = 'black') +
    geom_ribbon(aes(ymin = theta, ymax = fit),
                fill = "gray",
                alpha = .8) +
    geom_point(size = 1) +
    theme_bw() +
    ggtitle(paste0(
      "area =",
      round(area, 1),
      "     y = ",
      round(alpha, 3),
      "* exp(",
      round(beta, 3),
      "*x) +" ,
      round(theta, 3)
    )) +
    scale_x_continuous(name = 'Time after switching oxygen (s)',
                       #expand = c(0, 0),
                       limits = c(0, 155)) +
    scale_y_continuous(
      name = expression("A (\u03BCmol CO"[2] ~ m ^ -2 ~ s ^ -1 ~ ')'),
      #expand = c(0, 0),
      #limits = c(2, 7),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    theme(legend.position = "none") +
    plot_theme
  return (plot)
}



# combine traces from all reps for high (40% O2) to low oxygen (2% O2) transition (40-2)
datalist = list()
for (i in 1:10) {
  tem <- makeDf_40_2(data_40_2, rep_name[i])
  tem['id'] <- rep_name[i]  # Append new row
  datalist[[i]] <- tem
}
all_traces_40_2 = do.call(rbind, datalist)

# plot all traces with fitted lines
# sort id by A value
id_order <- all_traces_40_2 %>%
  filter(t == 100) %>%
  arrange(y) %>%
  select(id) %>%
  pull()

# add the order to a new column
all_traces_40_2 <- all_traces_40_2 %>%
  rowwise() %>%
  mutate(ordered = which(id_order == id))


plot_traces <- all_traces_40_2 %>%
  ggplot(aes(
    x = t,
    y = y,
    group = ordered,
    color = ordered
  )) +
  geom_point(size = .5) +
  # SSasymp: Self-Starting Nls Asymptotic Regression Model
  stat_smooth(
    method = "nls",
    formula = y ~ SSasymp(x, Asym, R0, lrc),
    se = FALSE,
    size = 0.5
  ) +
  theme_bw() +
  scale_x_continuous(name = 'Time after switching oxygen (s)',
                     #expand = c(0, 0),
                     limits = c(0, 205)) +
  scale_y_continuous(
    name = expression("A (\u03BCmol CO"[2] ~ m ^ -2 ~ s ^ -1 ~ ')'),
    #expand = c(0, 0),
    limits = c(2.5, 7),
    breaks = scales::pretty_breaks(n = 6)
  ) +
  #scale_colour_gradient(low = "grey", high = "black") +
  theme(legend.position = "none") +
  plot_theme


ggsave(
  paste0(out_path, "/figure_combined_traces_40_2.jpeg"),
  plot_traces,
  width = 6,
  height = 5,
  dpi = 600,
  units = "in"
)

# combine results from all reps for high (40% O2) to low oxygen (2% O2) transition (40-2)
datalist = list()

for (i in 1:10) {
  tem <-
    makeDf_40_2(data_40_2, rep_name[i]) %>% getArea_40_2()  # Create new row
  tem['id'] <- rep_name[i]                 # Append new row
  datalist[[i]] <- tem
}

all_area_40_2 = do.call(rbind, datalist)
write.csv(all_area_40_2, paste0(out_path, "/table_fitting_parameters_40_2.csv"))

# plot each shaded curve for high (40% O2) to low oxygen (2% O2) transition (40-2)
for (i in 1:10) {
  plot <- makeDf_40_2(data_40_2, rep_name[i]) %>% plotShaded_40_2()
  ggsave(
    paste0(out_path, "/figure_40_2_", rep_name[i], ".jpeg"),
    plot,
    width = 6,
    height = 5,
    dpi = 300,
    units = "in"
  )
}

