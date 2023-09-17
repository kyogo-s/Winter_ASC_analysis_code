# Kyogo Sakai

# Importing
library(tidyverse)

# Sorting Files
overview <- read.csv("../Overview4.csv")

# functions
normalise <- function(n){
    return((cumsum(n)/sum(cumsum(n)))/max(cumsum(n)/sum(cumsum(n))))
}

get_raw_data <- function(source, number, label){
    h1 <- source$abs[source$location == number & source$when == "before"]
    h2 <- source$abs[source$location == number & source$when == "during"]
    h3 <- source$abs[source$location == number & source$when == "after"]
    return(list(h1, h2, h3))
}

get_hist_data <- function(raw_data, j, breaks){
    breaks <- seq(0, max(unlist(raw_data)), l = breaks)
    h1 <- hist(raw_data[[1]], breaks = breaks, plot = FALSE)
    h2 <- hist(raw_data[[2]], breaks = breaks, plot = FALSE)
    h3 <- hist(raw_data[[3]], breaks = breaks, plot = FALSE)
    return(list(h1, h2, h3))
}

plot_cum_hist <- function(hist_data, number, label, plot = TRUE){
    h1 <- hist_data[[1]]
    h2 <- hist_data[[2]]
    h3 <- hist_data[[3]]

    hdiff <- h2$counts - (h1$counts + h3$counts)/2
    df <- data.frame(index = seq(1, length(h2$counts)),before = normalise(h1$counts), during = normalise(h2$counts), after = normalise(h3$counts))
    df <- gather(df, type, value, -index)

    dfalt <- data.frame(index = seq(1, length(h2$counts)), before = h1$counts, during = h2$counts, after = h3$counts)
    dfalt <- gather(dfalt, type, value, -index)
    fig1 <- ggplot(df, aes(x = index, y = value, color = type)) + geom_point() + geom_line(size = 2) + 
        labs(title = paste(label), x = "mEPSC Amplitude (pA)", y = "Normalised Cumulative Density")

    fig2 <- ggplot(dfalt, aes(x = index, y = value, color = type)) + geom_point() + geom_line(size = 2) + 
        labs(title = paste(label), x = "mEPSC Amplitude (pA)", y = "Density")
    if (plot){
        plot(fig2)
    }
}

KS_test <- function(raw_data, number, label){
    h1 <- raw_data[[1]]
    h2 <- raw_data[[2]]
    h3 <- raw_data[[3]]

    ks1_3 <- ks.test(h1, h3, alternative = "two.sided")
    ks1_2 <- ks.test(h1, h2, alternative = "two.sided")
    ks2_3 <- ks.test(h2, h3, alternative = "two.sided")
    ks2_avg <- ks.test(h2, (h1 + h3)/2, alternative = "two.sided")
    stats <- data.frame(ks1_3 = ks1_3$p.value, ks1_2 = ks1_2$p.value, ks2_3 = ks2_3$p.value, ks2_avg = ks2_avg$p.value, location = number)
    return(stats)
}

subtraction_analysis <- function(raw_data, hist_data, j, eps, breaks, label){
    breaks <- seq(0, max(unlist(raw_data)), l = breaks)

    h1 <- hist_data[[1]]$count/eps[1]
    h2 <- hist_data[[2]]$count/eps[2]
    h3 <- hist_data[[3]]$count/eps[3]

    h_avg <- (h1 + h3)/2
    hdiff <- h2 - h_avg

    df <- data.frame(index = breaks[1:48], hdiff = hdiff)
    fig1 <- ggplot(df, aes(x = index, y = hdiff)) + geom_point() + geom_line(size = 2) + 
        labs(title = paste(label), x = "Amplitude Histogram Difference, Normalised by Episodes", y = "during - AVG(before, after)")
    plot(fig1)

    #mean <- data.frame()

    return(list(df, mean(hdiff)))
}

# Execute commands

all_mean_amplitude <- data.frame()

for (i in seq(1, length(overview$name))){
    # Dealing with particular cell:
    filename <- paste('../', overview$name[i], sep = '')
    celltype <- overview$cell_type[i]
    num_eps <- strsplit(overview$episodes[i], split = '%')
    distances <- strsplit(overview$distance[i], split = ' ')
    no_locs <- overview$no_locs[i]

    df <- read.csv(paste(filename, '.csv', sep = ''))
    df <- gather(df, type, value)
    df <- na.omit(df)
    df$location <- as.factor(substring(df$type, 4, 4))
    df$when <- as.factor(substring(df$type, 6))
    df$abs <- abs(df$value)

    hdiff_df <- data.frame()
    KS_results <- data.frame()


    for (j in seq(1, length(levels(df$location)))) {
        # Dealing with particular location of cell:
        cell_label <- paste(filename, "Loc:", as.character(j), "cellType:", celltype)

        this_eps <- as.numeric(strsplit(as.character(num_eps[[1]]), " ")[[j]])

        this_dist <- as.numeric(distances[[1]])[j]

        raw_data <- get_raw_data(df, j, cell_label)
        hist_data <- get_hist_data(raw_data, j, 49)
        plot_cum_hist(hist_data, j, cell_label)
        KS_res <- KS_test(raw_data, j, cell_label)
        hdiff <- subtraction_analysis(raw_data, hist_data, j, this_eps, 49, cell_label)

        dist_mean <- data.frame(dist = this_dist, mean = hdiff[[2]], location = j, celltype = celltype, file = filename)

        hdiff_df <- rbind(hdiff_df, hdiff[[1]])
        all_mean_amplitude <- rbind(all_mean_amplitude, dist_mean)
    }
    print(paste("finished", filename))
}

a <- all_mean_amplitude %>%
    group_by(file) %>%
    summarise()
 



ggplot(all_mean_amplitude, aes(x = dist, y = mean, color = celltype)) + 
    geom_point() 

test <- all_mean_amplitude %>% 
    group_by(celltype, file) %>%
    arrange(dist) %>%
    mutate(norm_mean = mean/mean[1]) 


fig <- all_mean_amplitude %>% 
    group_by(celltype, file) %>%
    arrange(dist) %>%
    mutate(norm_mean = mean/mean[1]) %>%
    ggplot(aes(x = dist, y = norm_mean, color = celltype)) +
        geom_point(size = 2) + labs(title = "Mean amplitude by distance, normalised to first cell", x = "Distance (uM)", y = "Normalised mean mEPSC Amplitude (a.u.)")

plot(fig)