# 07_figure2.R -- redraw Figure 2 with the weighted-MI standardized risks
suppressMessages({library(ggplot2); library(dplyr); library(patchwork); library(scales)})
GRP  <- c("never","remitted","incident","persistent")
LAB  <- c(never="Never low", remitted="Remitted", incident="Incident", persistent="Persistent")
COL  <- c(never="#2f5f9e", remitted="#8fb2d9", incident="#e8a33d", persistent="#a53b2a")
SHAPE<- c(never=16, remitted=15, incident=17, persistent=18)
NS   <- c(never="3,833", remitted="325", incident="376", persistent="205")

## ---- panel A: observed outcomes ----
a <- data.frame(
  grp = factor(rep(c("Never low","Remitted","Incident","Persistent"), each=3),
               levels=c("Persistent","Incident","Remitted","Never low")),
  cat = factor(rep(c("Alive, ADL-free","ADL disability","Death"), 4),
               levels=c("Alive, ADL-free","ADL disability","Death")),
  val = c(87.8,7.7,4.4, 79.1,11.7,9.2, 73.4,15.2,11.4, 62.0,17.1,21.0))
labs <- data.frame(
  grp = factor(c("Never low","Remitted","Incident","Persistent"),
               levels=c("Persistent","Incident","Remitted","Never low")),
  n   = c("3,833","325","376","205"))
a$ylab <- factor(paste0(a$grp, "\n(n = ", labs$n[match(a$grp, labs$grp)], ")"),
                 levels=paste0(levels(a$grp), "\n(n = ", labs$n[match(levels(a$grp), labs$grp)], ")"))
pa <- ggplot(a, aes(x=val, y=ylab, fill=cat)) +
  geom_col(width=0.62, colour="white", linewidth=0.3, position=position_stack(reverse=TRUE)) +
  geom_text(aes(label=sprintf("%.1f", val), colour=cat), size=3.7, show.legend=FALSE,
            position=position_stack(reverse=TRUE, vjust=0.5), data=subset(a, cat!="Death")) +
  geom_text(data=subset(a, cat=="Death"), aes(x=102, y=ylab, label=sprintf("%.1f", val)),
            size=3.7, colour="#a53b2a", hjust=0, inherit.aes=FALSE) +
  scale_fill_manual(values=c("Alive, ADL-free"="#a8c4e0","ADL disability"="#e8a33d","Death"="#a53b2a"), name=NULL) +
  scale_colour_manual(values=c("Alive, ADL-free"="#1d3f66","ADL disability"="#7a4b00","Death"="#a53b2a"), guide="none") +
  scale_x_continuous(limits=c(0,110), breaks=seq(0,100,20), expand=c(0,0)) +
  labs(title="A. Observed five-year outcomes", x="Five-year outcome (%)", y=NULL) +
  theme_minimal(base_size=10, base_family="sans") +
  theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(colour="grey90", linewidth=0.35),
        plot.title=element_text(face="bold", size=11, margin=margin(b=6)),
        axis.text.y=element_text(size=9, colour="grey15"),
        axis.text.x=element_text(size=9, colour="grey25"),
        axis.title.x=element_text(size=9.5, margin=margin(t=4)),
        legend.position="bottom", legend.justification="center",
        legend.text=element_text(size=9), legend.key.size=unit(0.85,"lines"),
        plot.margin=margin(4,10,0,4))

## ---- panels B and C ----
ab <- read.csv("out/S2_weighted_MI_absolute_risk.csv")
ab <- ab[ab$cutoff==19, ]
ab$grp <- factor(ab$group, levels=GRP)
mk <- function(yv, ttl, ylab, ymax, brk){
  ggplot(ab, aes(age, .data[[yv]], colour=grp, shape=grp)) +
    geom_line(linewidth=0.75) + geom_point(size=2.1) +
    scale_colour_manual(values=COL, labels=LAB, name=NULL) +
    scale_shape_manual(values=SHAPE, labels=LAB, name=NULL) +
    scale_y_continuous(limits=c(0,ymax), breaks=brk, expand=c(0,0)) +
    scale_x_continuous(breaks=c(60,65,70,75,80), limits=c(57,83), expand=c(0,0)) +
    labs(title=ttl, x="Age at the 2013 landmark (years)", y=ylab) +
    theme_minimal(base_size=10, base_family="sans") +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_line(colour="grey91", linewidth=0.3),
          panel.grid.major.x=element_blank(),
          plot.title=element_text(face="bold", size=11, margin=margin(b=6)),
          axis.text=element_text(size=9, colour="grey25"),
          axis.title=element_text(size=9.5),
          axis.title.x=element_text(margin=margin(t=6)),
          axis.title.y=element_text(margin=margin(r=6)),
          legend.position=c(0.02,0.99), legend.justification=c(0,1),
          legend.background=element_blank(), legend.text=element_text(size=8.6),
          legend.key.height=unit(0.8,"lines"), legend.key.width=unit(1.1,"lines"),
          plot.margin=margin(4,10,4,4))
}
pb <- mk("disabled", "B. Standardized risk of ADL disability", "Five-year ADL disability risk (%)", 24, seq(0,20,5))
pc <- mk("dead", "C. Standardized risk of death", "Five-year mortality risk (%)", 32, seq(0,30,5))
pc <- pc + theme(legend.position="none")

fig <- pa / (pb | pc) + plot_layout(heights=c(0.92,1))
ggsave("figures/Figure2.png", fig, width=9.1, height=9.55, dpi=400, bg="white")
cat("saved figs/Figure2_rev_v2.png\n")
