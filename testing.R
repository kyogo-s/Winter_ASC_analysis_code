# Kyogo Sakai September 2023 
# Analysis files for 2023 Winter Session ASC

#Importing
library(tidyverse)

jul11_cell1 <- read.csv("17_07_2023_cell1_amplitudes.csv")
jul11_cell1 <- gather(jul11_cell1, type, value)
jul11_cell1 <- na.omit(jul11_cell1)

jul11_cell1$location <- substring(jul11_cell1$type, 4, 4)
jul11_cell1$when <- as.factor(substring(jul11_cell1$type, 6))
jul11_cell1$abs <- abs(jul11_cell1$value)

# testing

normalise <- function(n){
    return((cumsum(n)/sum(cumsum(n)))/max(cumsum(n)/sum(cumsum(n))))
}

plot_cum_hist <- function(source, number, label){
    breaks = seq(0, 70, l = 49)
    # Plot the three histograms
    h1 <- hist(source$abs[source$location == number & source$when == "before"], breaks = breaks, plot = FALSE)
    h2 <- hist(source$abs[source$location == number & source$when == "during"], breaks = breaks, plot = FALSE)
    h3 <- hist(source$abs[source$location == number & source$when == "after"], breaks = breaks, plot = FALSE)

    # Plot Subtracted histogram
    hdiff <- h2$counts - (h1$counts + h3$counts)/2

    # Make Plots

    df <- data.frame(before = normalise(h1$counts), during = normalise(h2$counts), after = normalise(h3$counts))
    df$index <- breaks[1:48]
    df <- gather(df, type, value, -index)

    fig <- ggplot(df, aes(x = index, y = value, color = type)) + geom_point() + geom_line(size = 2) + 
        labs(title = paste(label), x = "mEPSC Amplitude (pA)", y = "Cumulative Density")
    plot(fig)

    df <- data.frame(diff = hdiff, index = breaks[1:48])
    fig <- ggplot(df, aes(x = index, y = diff)) + geom_point() + geom_line(size = 2) + 
        labs(title = paste(label), x = "mEPSC Amplitude (pA)", y = "during - AVG(before, after)") + scale_y_continuous(limits = c(-10, 120))
    plot(fig)

    # KS Test
    ks1_3 <- ks.test(h1$counts, h3$counts, alternative = "two.sided")
    ks1_2 <- ks.test(h1$counts, h2$counts, alternative = "two.sided")
    ks2_3 <- ks.test(h2$counts, h3$counts, alternative = "two.sided")
    #print(ks1_3)
    #print(ks1_2)
    #print(ks2_3)
}

# Execute commands
for (i in seq(1, 6)) {
    plot_cum_hist(jul11_cell1, i, paste("July 11th", as.character(i)))
}
