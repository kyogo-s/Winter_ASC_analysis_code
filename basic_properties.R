library(tidyverse)
library(ggsignif)

df <- read.csv("../electrophysiological.csv", stringsAsFactors = TRUE)

df <- gather(df, measure, value, input.resistance:resting)
df <- na.omit(df)

ggplot(subset(df, df$measure == "input.resistance"), aes(x = celltype, y = value, fill = celltype)) + 
    geom_bar(stat = "summary", fun.y = "mean", width = 0.5) + 
    geom_jitter(color = "#ff000088") + 
    geom_errorbar(stat = "summary", width = 0.2) + 
    labs(x = "Cell Type", y = "Input Resistance (MOhm)") + 
    scale_fill_discrete(name = "Cell Type", labels = c("SL", "SP")) + 
    geom_signif(
        annotation = formatC("****"),
        y_position = 600, xmin = 1, xmax = 2,
        tip_length = 0
    )

ggsave("SL_SP_input_resistance.png", width = 4, height = 4, dpi = 300)

ggplot(subset(df, df$measure == "resting"), aes(x = celltype, y = value, fill = celltype)) + 
    geom_bar(stat = "summary", fun.y = "mean", width = 0.5) + 
    geom_jitter(color = "#ff000088") + 
    geom_errorbar(stat = "summary", width = 0.2) + 
    labs(x = "Cell Type", y = "Resting Potential (mV)") + 
    scale_fill_discrete(name = "Cell Type", labels = c("SL", "SP")) + 
    scale_y_continuous(trans = "reverse") +
    geom_signif(
        annotation = formatC("****"),
        y_position = 90, xmin = 1, xmax = 2,
        tip_length = 0
    )

ggsave("SL_SP_resting.png", width = 4, height = 4, dpi = 300)

t.test(df$value[df$measure == "input.resistance" & df$celltype == 'SP'], df$value[df$measure == "input.resistance" & df$celltype == 'SL'])
t.test(df$value[df$measure == "resting" & df$celltype == 'SP'], df$value[df$measure == "resting" & df$celltype == 'SL'])

