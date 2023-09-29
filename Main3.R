# Kyogo Sakai

# Importing
library(tidyverse)
library(ggsignif)
theme_set(theme_classic())

# Sorting Files
overview <- read.csv("../dataset2/Overview4.csv")

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

plot_cum_hist <- function(hist_data, number, label, eps, plot = TRUE){
    h1 <- hist_data[[1]]$count
    h2 <- hist_data[[2]]$count
    h3 <- hist_data[[3]]$count

    h1 <- hist_data[[1]]$count/eps[1]
    h2 <- hist_data[[2]]$count/eps[2]
    h3 <- hist_data[[3]]$count/eps[3]

    hdiff <- h2 - (h1 + h3)/2
    df <- data.frame(index = seq(1, length(h2)),before = normalise(h1), during = normalise(h2), after = normalise(h3))
    df <- gather(df, type, value, -index)

    dfalt <- data.frame(index = seq(1, length(h2)), before = h1, during = h2, after = h3)
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
    stats <- data.frame(ks1_3 = ks1_3$p.value, ks1_2 = ks1_2$p.value, ks2_3 = ks2_3$p.value, ks2_avg = ks2_avg$p.value, location = number, label = label)
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

    mean <- mean_from_hist(breaks[1:48], hdiff)

    return(list(df, mean))
}

mean_from_hist <- function(bins, freqs){
    #return(sum(bins * freqs)/sum(freqs))

    # # Alternative summary statistic - if using mode instead of mean
    # df <- data.frame(bins = bins, freqs = freqs)
    # df <- arrange(df, freqs)
    # return(df$bins[length(df$bins)])

    # alternative mean if capping a max
    bins <- bins[1:30]
    freqs <- freqs[1:30]

    return(sum(bins * freqs)/sum(freqs))

}

# Execute commands
KS_results_all <- data.frame()

all_mean_amplitude <- data.frame()

for (i in seq(1, length(overview$name))){
    # Dealing with particular cell:
    filename <- paste('../dataset2/', overview$name[i], sep = '')
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
        plot_cum_hist(hist_data, j, cell_label, this_eps)
        KS_res <- KS_test(raw_data, j, cell_label)
        KS_results <- rbind(KS_results, KS_res)
        hdiff <- subtraction_analysis(raw_data, hist_data, j, this_eps, 49, cell_label)

        dist_mean <- data.frame(dist = this_dist, mean = hdiff[[2]], location = j, celltype = celltype, file = filename)

        hdiff_df <- rbind(hdiff_df, hdiff[[1]])
        all_mean_amplitude <- rbind(all_mean_amplitude, dist_mean)
    }
    print(paste("finished", filename))
    print(KS_results)
    KS_results_all <- rbind(KS_results_all, KS_results)
}

    write.csv(KS_results_all, "KS_results_all.csv")


a <- all_mean_amplitude %>%
    group_by(file) %>%
    arrange(dist) %>%
    summarise(furthest = mean[length(mean)], celltype = celltype[length(celltype)], dist = dist[length(dist)])

# Fig 1a
ggplot(a, aes(x = as.factor(celltype), y = furthest, fill = as.factor(celltype)), size = 2) + 
    geom_bar(stat = "summary", fun.y = "mean", width = 0.5) + 
    geom_jitter(color = '#00000088') + 
    geom_errorbar(stat = 'summary', width = 0.2) + 
    labs(x = "Cell Type", y = "Mean Amplitude (pA)") + 
    scale_x_discrete(labels = c("SL", "SP")) + 
    scale_fill_discrete(name = "Cell Type", labels = c("SL", "SP")) +
    geom_signif(
        annotation = formatC("*"),
        y_position = 16, xmin = 1, xmax = 2,
        tip_length = 0
    )

ggsave("SP_SL_bar.png", width = 4, height = 4, dpi = 300)

a %>%
    group_by(celltype) %>%
    summarise(mean = mean(furthest), sd = sd(furthest), sem = sd(furthest)/sqrt(length(furthest))) %>%
    ggplot(aes(x = celltype, y = mean, fill = celltype)) + 
        geom_bar(stat = "identity") + 
        geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) + 
        labs(title = "Mean amplitude of furthest location", x = "Cell Type", y = "Mean Amplitude (pA)")

ggplot(a, aes(x = as.factor(celltype), y = furthest, fill = as.factor(celltype))) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(color = 'red') +
    labs(title = "Mean amplitude of furthest location", x = "Cell Type", y = "Mean Amplitude (pA)")

ggplot(a, aes(x = as.factor(celltype), y = furthest, fill = as.factor(celltype))) + 
    geom_bar(stat = "identity") + 
    labs(title = "Mean amplitude of furthest location", x = "Cell Type", y = "Mean Amplitude (pA)")   

ggplot(all_mean_amplitude, aes(x = dist, y = mean, color = celltype)) + 
    geom_point() 

ggsave("dist_amp_scatter.png", width = 4, height = 4, dpi = 300)

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

res <- t.test(a$furthest[a$celltype == "SL"], a$furthest[a$celltype == "SP"])
res

for (file in unique(all_mean_amplitude$file)){
    data <- all_mean_amplitude[all_mean_amplitude$file == file, ] 
    print(file)
    plot <- ggplot(data, aes(x = dist, y = mean)) + geom_point(color = 'red', size = 5) + 
        labs(x = "Distance (um)", y = "Mean mEPSC Amplitude (pA)", title = file)
    plot(plot)
}