
trim_trace <- function(df){
    df <- df[df$time > 1.3 & df$time < 2.3, ]
    return(df)
}

# 11 Jul cell 1 Loc 4 trace number 2 and 21
puff <- read.csv("../sample_trace_test_w_puff.csv")
nopuff <- read.csv("../sample_trace_test_no_puff.csv")

ggplot(trim_trace(puff), aes(x = time, y = current)) + geom_line() + 
    scale_y_continuous(limits = c(0e-12, 60e-12))

ggplot(trim_trace(nopuff), aes(x = time, y = current)) + geom_line() + 
    scale_y_continuous(limits = c(0e-12, 60e-12))
