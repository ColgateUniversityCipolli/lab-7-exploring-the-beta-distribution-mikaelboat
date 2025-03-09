############################
# LAB 07
############################
# load libraries
library(tidyverse)
library(e1071)
library(cumstats)
library(patchwork)
############################
# task one: describe population distribution

# functions to calculate mean, variance, skewness and kurtosis
compute <- function(alpha, beta){
  mean <- alpha / (alpha + beta)
  variance <- (alpha*beta)/(((alpha+beta)^2)*(alpha+beta+1))
  skewness <- (2*(beta-alpha)*sqrt(alpha+beta+1))/((alpha+beta+2)*sqrt(alpha*beta))
  kurtosis <- (6*((alpha-beta)^2*(alpha+beta+1)-alpha*beta*(alpha+beta+2)))/(alpha*beta*(alpha+beta+2)*(alpha+beta+3))
  
  rowtoAdd <- c(alpha, beta, mean, variance, skewness, kurtosis)
  names(rowtoAdd) <- c("alpha", "beta", "mean", "variance", "skewness", "kurtosis")
  return(rowtoAdd)
}

# distributions
main.df <- tibble(
  alpha = numeric(), beta = numeric(), mean = numeric(),
  variance = numeric(), skewness = numeric(), kurtosis = numeric()
)
# Beta (alpha = 2, beta = 5)
caseOne <- compute(2,5)
main.df <- bind_rows(main.df, caseOne)
# Beta (alpha = 5, beta = 5)
caseTwo <- compute(5, 5)
main.df <- bind_rows(main.df, caseTwo)
# Beta (alpha = 5, beta = 2)
caseThree <- compute(5, 2)
main.df <- bind_rows(main.df, caseThree)
# Beta (alpha = 0.5, beta = 0.5)
caseFour <- compute(0.5, 0.5)
main.df <- bind_rows(main.df, caseFour)

#########
# plots
# beta(2,5)
alpha <- 2
beta <- 5
beta2_5 <- tibble(x = seq(0, 1, length.out=1000))|>   # generate a grid of points
  mutate(beta.pdf = dbeta(x, alpha, beta),                      # compute the beta PDF
         norm.pdf = dnorm(x,                                    # Gaussian distribution with
                          mean = alpha/(alpha+beta),            # same mean and variance
                          sd = sqrt((alpha*beta)/((alpha+beta)^2*(alpha+beta+1)))))

figure1 <- ggplot(data= beta2_5) +                                              # specify data
  geom_line(aes(x=x, y=beta.pdf, color="Beta (2, 5)")) +                 # plot beta dist
  geom_line(aes(x=x, y=norm.pdf, color= "Normal (0.285, 0.0255)")) +  # plot gaussian dist
geom_hline(yintercept = 0)+                                            # plot x axis
  theme_bw()+                                                          # change theme
  xlab("x")+                                                           # label x axis
  ylab("Density")+                                                     # label y axis
  scale_color_manual("", values = c("black", "grey"))+                 # change colors
  theme(legend.position = "bottom") 


# beta(5,5)
alpha <- 5
beta <- 5
beta5_5 <- tibble(x = seq(0, 1, length.out=1000))|>   # generate a grid of points
  mutate(beta.pdf = dbeta(x, alpha, beta),                      # compute the beta PDF
         norm.pdf = dnorm(x,                                    # Gaussian distribution with
                          mean = alpha/(alpha+beta),            # same mean and variance
                          sd = sqrt((alpha*beta)/((alpha+beta)^2*(alpha+beta+1)))))

figure2 <- ggplot(data = beta5_5) +                                              # specify data
  geom_line(aes(x=x, y=beta.pdf, color="Beta (5, 5)")) +                 # plot beta dist
  geom_line(aes(x=x, y=norm.pdf, color= "Normal (0.5, 0.022)")) +  # plot gaussian dist
  geom_hline(yintercept = 0)+                                            # plot x axis
  theme_bw()+                                                          # change theme
  xlab("x")+                                                           # label x axis
  ylab("Density")+                                                     # label y axis
  scale_color_manual("", values = c("black", "grey"))+                 # change colors
  theme(legend.position = "bottom") 


# beta(5, 2)
alpha <- 5
beta <- 2
beta5_2 <- tibble(x = seq(0, 1, length.out=1000))|>   # generate a grid of points
  mutate(beta.pdf = dbeta(x, alpha, beta),                      # compute the beta PDF
         norm.pdf = dnorm(x,                                    # Gaussian distribution with
                          mean = alpha/(alpha+beta),            # same mean and variance
                          sd = sqrt((alpha*beta)/((alpha+beta)^2*(alpha+beta+1)))))

figure3 <- ggplot(data= beta5_2) +                                              # specify data
  geom_line(aes(x=x, y=beta.pdf, color="Beta (5, 2)")) +                 # plot beta dist
  geom_line(aes(x=x, y=norm.pdf, color= "Normal (0.714, 0.0255)")) +  # plot gaussian dist
  geom_hline(yintercept = 0)+                                            # plot x axis
  theme_bw()+                                                          # change theme
  xlab("x")+                                                           # label x axis
  ylab("Density")+                                                     # label y axis
  scale_color_manual("", values = c("black", "grey"))+                 # change colors
  theme(legend.position = "bottom") 


# beta (0.5, 0.5)
alpha <- 0.5
beta <- 0.5
beta0.5_0.5 <- tibble(x = seq(0, 1, length.out=1000))|>   # generate a grid of points
  mutate(beta.pdf = dbeta(x, alpha, beta),                      # compute the beta PDF
         norm.pdf = dnorm(x,                                    # Gaussian distribution with
                          mean = alpha/(alpha+beta),            # same mean and variance
                          sd = sqrt((alpha*beta)/((alpha+beta)^2*(alpha+beta+1)))))

figure4 <- ggplot(data= beta0.5_0.5) +                                              # specify data
  geom_line(aes(x=x, y=beta.pdf, color="Beta (0.5, 0.5)")) +                 # plot beta dist
  geom_line(aes(x=x, y=norm.pdf, color= "Normal (0.5, 0.125)")) +  # plot gaussian dist
  geom_hline(yintercept = 0)+                                            # plot x axis
  theme_bw()+                                                          # change theme
  xlab("x")+                                                           # label x axis
  ylab("Density")+                                                     # label y axis
  scale_color_manual("", values = c("black", "grey"))+                 # change colors
  theme(legend.position = "bottom") 


############################
# task two: compute moments
beta.summary <- function(alpha, beta, k, centered){
  if(centered){
    integrand <- function(x) (x)*dbeta(x,alpha,beta)
    mu_x = integrate(integrand, lower = 0, upper = 1)$value
    
    integrand <-  function(x) (x - mu_x)^k * dbeta(x,alpha,beta)
    kth.uncentered = integrate(integrand, lower = 0, upper = 1)
    
  }
  else{
    integrand <- function(x) (x^k)*dbeta(x,alpha,beta)
    kth.centered = integrate(integrand, lower = 0, upper = 1)
  }
}

# computing beta(2,5) population-level characteristics using beta.summary()
mean <- beta.summary(2,5,1,F)$value
variance <- beta.summary(2,5,2,T)$value 
skewness <- (beta.summary(2,5,3,T)$value) / ((beta.summary(2,5,2,T)$value)^(3/2))
kurtosis <- (beta.summary(2,5,4,T)$value) / ((beta.summary(2,5,2,T)$value)^(2)) - 3


############################
# task three
numeric.summary <- tibble(
  alpha = numeric(), beta = numeric(), mean = numeric(),
  variance = numeric(), skewness = numeric(), kurtosis = numeric()
)

set.seed(7272) # Set seed so we all get the same results.
sample.size <- 500 # Specify sample details


# beta(2,5)
alpha <- 2
beta <- 5
caseOneSample <- rbeta(n = sample.size,  # sample size
                     shape1 = alpha,   # alpha parameter
                     shape2 = beta)    # beta parameter

caseOne.summary <- data.frame(caseOneSample) |>
  summarize(alpha = alpha,
            beta = beta,
            mean = mean(caseOneSample),
            variance = var(caseOneSample),
            skewness = skewness(caseOneSample),
            kurtosis = kurtosis(caseOneSample))

numeric.summary <- bind_rows(numeric.summary, caseOne.summary)



p1 <- ggplot(data.frame(caseOneSample), aes(x = caseOneSample)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 15, 
                 fill = "lightblue",
                 color = "black") +
  stat_density(aes(color = "Estimated Density"),
               linewidth = 0.65,
               geom = "line") + 
  stat_function(fun = dbeta,
                args = list(shape1 = alpha, shape2 = beta),
                aes(color = "Beta (2,5)"),
                linewidth = 1) + 
  scale_x_continuous(limits=c(0, 1)) +
  labs(x = "Value",
       y = "Frequency",
       color = "Legend") +
  scale_color_manual(values = c("Estimated Density" = "black", 
                               "Beta (2,5)" = "grey")) +
  theme_minimal()


# beta(5,5)
alpha <- 5
beta <- 5
caseTwoSample <- rbeta(n = sample.size,  # sample size
                       shape1 = alpha,   # alpha parameter
                       shape2 = beta)    # beta parameter


caseTwo.summary <- data.frame(caseTwoSample) |>
  summarize(alpha = alpha,
            beta = beta,
            mean = mean(caseTwoSample),
            variance = var(caseTwoSample),
            skewness = skewness(caseTwoSample),
            kurtosis = kurtosis(caseTwoSample))

numeric.summary <- bind_rows(numeric.summary, caseTwo.summary)


p2 <- ggplot(data.frame(caseTwoSample), aes(x = caseTwoSample)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 15, 
                 fill = "lightblue",
                 color = "black") +
  stat_density(aes(color = "Estimated Density"), geom = "line",
               linewidth = 0.65) + 
  stat_function(fun = dbeta,
                args = list(shape1 = alpha, shape2 = beta),
                aes(color = "Beta (5,5)"),
                linewidth = 1) + 
  scale_x_continuous(limits=c(0, 1)) +
  labs(x = "Value",
       y = "Frequency",
       color = "Legend") +
  scale_color_manual(values = c("Estimated Density" = "black", 
                                "Beta (5,5)" = "grey")) +
  theme_minimal()


# beta(5, 2)
alpha <- 5
beta <- 2
caseThreeSample <- rbeta(n = sample.size,  # sample size
                       shape1 = alpha,   # alpha parameter
                       shape2 = beta)    # beta parameter

caseThree.summary <- data.frame(caseThreeSample) |>
  summarize(alpha = alpha,
            beta = beta,
            mean = mean(caseThreeSample),
            variance = var(caseThreeSample),
            skewness = skewness(caseThreeSample),
            kurtosis = kurtosis(caseThreeSample))

numeric.summary <- bind_rows(numeric.summary, caseThree.summary)


p3 <- ggplot(data.frame(caseThreeSample), aes(x = caseThreeSample)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 15, 
                 fill = "lightblue",
                 color = "black") +
  stat_density(aes(color = "Estimated Density"), geom = "line",
               linewidth = 0.65) + 
  stat_function(fun = dbeta,
                args = list(shape1 = alpha, shape2 = beta),
                aes(color = "Beta (5,2)"),
                linewidth = 1) + 
  scale_x_continuous(limits=c(0, 1)) +
  labs(x = "Value",
       y = "Frequency",
       color = "Legend") +
  scale_color_manual(values = c("Estimated Density" = "black", 
                                "Beta (5,2)" = "grey")) +
  theme_minimal()


# beta(0.5,0.5)
alpha <- 0.5
beta <- 0.5
caseFourSample <- rbeta(n = sample.size,  # sample size
                       shape1 = alpha,   # alpha parameter
                       shape2 = beta)    # beta parameter

caseFour.summary <- data.frame(caseFourSample) |>
  summarize(alpha = alpha,
            beta = beta,
            mean = mean(caseFourSample),
            variance = var(caseFourSample),
            skewness = skewness(caseFourSample),
            kurtosis = kurtosis(caseFourSample))

numeric.summary <- bind_rows(numeric.summary, caseFour.summary)


p4 <- ggplot(data.frame(caseFourSample), aes(x = caseFourSample)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 15, 
                 fill = "lightblue",
                 color = "black") +
  stat_density(aes(color = "Estimated Density"), geom = "line",
               linewidth = 0.65) + 
  stat_function(fun = dbeta,
                args = list(shape1 = alpha, shape2 = beta),
                aes(color = "Beta (0.5,0.5)"),
                #color = "grey", 
                linewidth = 1) + 
  scale_x_continuous(limits=c(0, 1)) +
  labs(x = "Value",
       y = "Frequency",
       color = "Legend") +
  scale_color_manual(values = c("Estimated Density" = "black", 
                                "Beta (0.5,0.5)" = "grey")) +
  theme_minimal()



##########################
# task four
cum.numerical <- data.frame(
  mean = cummean(caseOneSample),
  variance = cumvar(caseOneSample),
  skewness = cumskew(caseOneSample),
  kurtosis = cumkurt(caseOneSample) - 3 # excess cumulative kurtosis
)

# mean plot
mean_p <- ggplot(cum.numerical, aes(x = 1:nrow(cum.numerical), y = mean)) +
  geom_line() +
  geom_hline(yintercept = main.df[[1, "mean"]]) + # get true value from our tibble which contains it
  labs(x = "x",
        title = "Cumulative Mean")

# variance plot
variance_p <- ggplot(cum.numerical, aes(x = 1:nrow(cum.numerical), y = variance)) +
  geom_line() +
  geom_hline(yintercept = main.df[[1, "variance"]]) +
  labs(x = "x",
        title = "Cumulative Variance")

# skewness plot
skewness_p <- ggplot(cum.numerical, aes(x = 1:nrow(cum.numerical), y = skewness)) +
  geom_line() +
  geom_hline(yintercept = main.df[[1, "skewness"]]) +
  labs(x = "x",
        title = "Cumulative Skewness")

# kurtosis plot
kurtosis_p <- ggplot(cum.numerical, aes(x = 1:nrow(cum.numerical), y = kurtosis)) +
  geom_line() +
  geom_hline(yintercept = main.df[[1, "kurtosis"]]) +
  labs(x = "x",
        title = "Cumulative Kurtosis")

# 2x2 grid
combined_plot <- (mean_p + variance_p) / (skewness_p + kurtosis_p)

# writing a for() loop to simulate new data

for (i in 2:50){
  
  set.seed(7272 + i)
  new.data <- rbeta(n = 500, shape1 = 2, shape2 = 5)
  
  cum.data <- data.frame(
    mean = cummean(new.data),
    variance = cumvar(new.data),
    skewness = cumskew(new.data),
    kurtosis = cumkurt(new.data) - 3 # excess cumulative kurtosis
  )
  
  mean_p = mean_p +
    geom_line(data = cum.data, aes(x = 1:nrow(cum.data), y = mean), color = i)
  variance_p = variance_p +
    geom_line(data = cum.data, aes(x = 1:nrow(cum.data), y = variance), color = i)
  skewness_p = skewness_p +
    geom_line(data = cum.data, aes(x = 1:nrow(cum.data), y = skewness), color = i)
  kurtosis_p = kurtosis_p +
    geom_line(data = cum.data, aes(x = 1:nrow(cum.data), y = kurtosis), color = i)

}

combined_plot2 <- (mean_p + variance_p) / (skewness_p + kurtosis_p)
