library(ggplot2)
library(RColorBrewer)

# df: a MetaData matrix with scores/gene expressions
# SelectedValue: value to plot

refcolor <- brewer.pal(5, 'Greys')[-1]
names(refcolor) <-
  c("FCA_Sampled",
    "HCA_CB_Sampled",
    "JUV_Sampled",
    "HCA_BM_Sampled")
samplecolor = rev(
  c(
    rev(refcolor),
    D117 = "#0058b4",
    D129 = "#2188c9",
    D217 = "#fbbb25",
    I217 = "#fca349",
    D213 = "#ff6b36",
    D124 = "#e34e2e",
    D123 = "#c33126",
    D360 = "#a41220"
  )
)
df$Donor <- factor(df$Donor, levels = names(samplecolor))

ggplot(df,
       aes_string(
         y = SelectedValue,
         x = 'Donor',
         color = 'Donor',
         fill = 'Donor'
       )) +
  geom_violin() +
  theme_classic() +
  scale_color_manual(values = samplecolor) +
  scale_fill_manual(values = samplecolor) +
  theme(axis.text.x = element_text(colour = 'black')) +
  theme(axis.text.y = element_text(colour = 'black')) +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0.5
  ))
