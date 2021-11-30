##### Code used to generate the main text figures along with the data used
#The data has been extracted from output files from GREML analyses and compiled into summary files

library(tidyverse)
library(data.table)
library(ggpubr)


dir <- "${workdir}"
setwd(dir)


######  Figure 1  ######

recap_arrays <- fread("01_Arrays.hsq", h=T)
recap_arrays$LD <- factor(recap_arrays$LD, levels(factor(recap_arrays$LD))[c(2,1,3)])


ann_text <- data.frame(Variance = 0.48,Bin = "0.01 - 0.1",lab = "Text",
                       LD = factor("High LD", levels = levels(factor(recap_arrays$LD))),
                       Array="GSA", trait="Height")


plot_array_a <- recap_arrays %>% filter(Trait=="Height") %>% ggplot(aes(y=Variance, x=Bin, group=Array, fill=Array)) +
  geom_bar(position ="dodge", stat="identity") +
  facet_grid(~LD,  space = "free_x",  scale="free_x") +
  labs(title="Height", x ="MAF", y = expression(Variance %+-% SE)) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_discrete(name = "Array") +
  geom_errorbar(aes(ymin=Variance-SE, ymax=Variance+SE), color="black", width=.2, position=position_dodge(1)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.background = element_rect(fill="white", size=0.5, linetype="solid", color="black"),
        legend.position = c(0.2, 0.7),
        plot.title = element_text(hjust = 0.5),
        strip.background =element_rect(fill="white")) +
  geom_text(data = ann_text, label = expression(paste(italic("h"["G + Imp"]^"2" %~~% 0.50 - 0.56))))


ann_text <- data.frame(Variance = 0.19,Bin = "0.01 - 0.1",lab = "Text",
                       LD = factor("High LD", levels = levels(factor(recap_arrays$LD))),
                       Array="GSA", trait="BMI")


plot_array_b <- recap_arrays %>% filter(Trait=="BMI") %>% ggplot(aes(y=Variance, x=Bin, group=Array, fill=Array)) +
  geom_bar(position ="dodge", stat="identity") +
  facet_grid(~LD,  space = "free_x",  scale="free_x") +
  labs(title="BMI", x ="MAF", y = expression(Variance %+-% SE)) + 
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_discrete(name = "Array") +
  geom_errorbar(aes(ymin=Variance-SE, ymax=Variance+SE), color="black", width=.2, position=position_dodge(1)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white"),
        legend.background = element_rect(fill="white", size=0.5, linetype="solid", color="black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  geom_text(data = ann_text, label = expression(paste(italic("h"["G + Imp"]^"2" %~~% 0.16 - 0.21))))

plot_01_arrays <- ggarrange(plot_array_a, plot_array_b, labels=c("A", "B"), common.legend = F)



######  Figure 2  ######

ldms_split <- fread("./02_WGS.hsq", h=T)

ldms_split$PC <- factor(ldms_split$PC, levels(factor(ldms_split$PC))[c(3,1,2)])
ldms_split$Trait <- factor(ldms_split$Trait, levels(factor(ldms_split$Trait))[c(2,1)])
ldms_split$Split <- factor(ldms_split$Split, levels(factor(ldms_split$Split))[c(2,1)])

PCcolors_2pc <- setNames(c("#F8766D", "#00BA38"), levels(ldms_split$PC)[c(1:2)])

ldms_3bins_plot_a <- ldms_split %>% filter(Trait  == "Height" & Split == "Tertiles") %>% ggplot(aes(y=Variance, x=LD, group=PC, fill=PC)) +
  geom_bar(position ="dodge", stat="identity") +
  facet_grid(~Bin,  space = "free_x",  scale="free_x") +
  labs(title="Height", x ="LD grouping", y = expression(Variance %+-% SE)) + 
  expand_limits(y=c(NA, 0.8)) +  scale_y_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_fill_manual(values = PCcolors_2pc, drop=T) +
  geom_errorbar(aes(ymin=Variance-SE, ymax=Variance+SE), color="black", width=.2, position=position_dodge(1)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white"),
        legend.background = element_rect(fill="white", size=0.5, linetype="solid", color="black"),
        legend.position = c(0.3, 0.6),
        plot.title = element_text(hjust = 0.5)) 


ldms_3bins_plot_b <- ldms_split %>% filter(Trait  == "BMI" & Split == "Tertiles")  %>% ggplot(aes(y=Variance, x=LD, group=PC, fill=PC)) +
  geom_bar(position ="dodge", stat="identity") +
  facet_grid(~Bin,  space = "free_x",  scale="free_x") +
  labs(title="BMI", x ="LD grouping", y = expression(Variance %+-% SE)) + 
  expand_limits(y=c(NA, 0.5)) + scale_y_continuous(breaks = seq(0, 0.5, 0.1)) +
  scale_fill_manual(values = PCcolors_2pc) +
  geom_errorbar(aes(ymin=Variance-SE, ymax=Variance+SE), color="black", width=.2, position=position_dodge(1)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))

PCcolors <- setNames(c("#F8766D", "#00BA38", "#619CFF"), levels(ldms_split$PC))


ldms_4bins_plot_a <- ldms_split %>% filter(Trait  == "Height" & Split == "Quartiles") %>% ggplot(aes(y=Variance, x=LD, group=PC, fill=PC)) +
  geom_bar(position ="dodge", stat="identity") +
  facet_grid(~Bin,  space = "free_x",  scale="free_x") +
  labs(title="Height", x ="LD grouping", y = expression(Variance %+-% SE)) + 
  expand_limits(y=c(NA, 0.8)) +  scale_y_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_fill_manual(values = PCcolors) +
  geom_errorbar(aes(ymin=Variance-SE, ymax=Variance+SE), color="black", width=.2, position=position_dodge(1)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white"),
        legend.background = element_rect(fill="white", size=0.5, linetype="solid", color="black"),
        legend.position = c(0.3, 0.6),
        plot.title = element_text(hjust = 0.5))

ldms_4bins_plot_b <- ldms_split %>% filter(Trait  == "BMI" & Split == "Quartiles")  %>% ggplot(aes(y=Variance, x=LD, group=PC, fill=PC)) +
  geom_bar(position ="dodge", stat="identity") +
  facet_grid(~Bin,  space = "free_x",  scale="free_x") +
  labs(title="BMI", x ="LD grouping", y = expression(Variance %+-% SE)) + 
  expand_limits(y=c(NA, 0.5)) + scale_y_continuous(breaks = seq(0, 0.5, 0.1)) +
  scale_fill_manual(values = PCcolors) +
  geom_errorbar(aes(ymin=Variance-SE, ymax=Variance+SE), color="black", width=.2, position=position_dodge(1)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) 

ga <-  ggarrange(ldms_3bins_plot_a, ldms_3bins_plot_b,  labels=c("A", "B"), common.legend=F)
gb <- ggarrange(ldms_4bins_plot_a, ldms_4bins_plot_b, labels=c("C", "D"), common.legend = F)
plot_02_WGS <- ggarrange(ga, gb, common.legend = F, nrow = 2)


######  Figure 3  ######


ldms_enriched <- fread("./03_Enrichment.hsq", h=T)
ldms_enriched$LD <- factor(ldms_enriched$LD, levels(factor(ldms_enriched$LD))[c(4,3,2,1,5)])
ldms_enriched$PC <- factor(ldms_enriched$PC, levels(factor(ldms_enriched$PC))[c(2,1)])

tploth48 <-  ldms_enriched %>% filter(!LD == "Total") %>% filter(Trait == "Height") %>% filter(PC == "48 PCs") %>%  arrange(VarperSNP) %>% select(Bin)
tplotb48 <-  ldms_enriched %>% filter(!LD == "Total") %>% filter(Trait == "BMI") %>% filter(PC == "48 PCs") %>%  arrange(VarperSNP) %>% select(Bin)


ldms_enriched_plot_a <- ldms_enriched %>% filter(!LD == "Total") %>% filter(Trait == "Height") %>% filter(PC == "48 PCs") %>%
  arrange(VarperSNP) %>%
  unite("BinLD", c("Bin", "LD"), remove =FALSE) %>%
  mutate(BinLD=factor(BinLD, levels=BinLD)) %>%
  ggplot(aes(y=VarperSNP, x=BinLD, fill=LD)) +
  geom_bar(position ="dodge", stat="identity") +
  labs(title="Height 48PCs", x ="MAF", y = expression(Variance/variant %+-% SE)) +
  guides(fill = guide_legend(title = "Bin"), color = guide_legend(title = "Bin")) +
  geom_errorbar(aes(ymin=VarperSNP - SEperSNP, ymax=VarperSNP + SEperSNP), color="black", width=.2, position=position_dodge(1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        legend.background = element_rect(fill="white", size=0.5, linetype="solid", color="black"),
        legend.position = c(0.6, 0.4)) +
  scale_x_discrete(labels = tploth48$Bin) +
  coord_flip()

ldms_enriched_plot_b <- ldms_enriched %>% filter(!LD == "Total") %>% filter(Trait == "BMI") %>% filter(PC == "48 PCs") %>%
  arrange(VarperSNP) %>%
  unite("BinLD", c("Bin", "LD"), remove =FALSE) %>%
  mutate(BinLD=factor(BinLD, levels=BinLD)) %>%
  ggplot(aes(y=VarperSNP, x=BinLD, fill=LD)) +
  geom_bar(position ="dodge", stat="identity") +
  labs(title="BMI 48PCs", x ="MAF", y = expression(Variance/variant %+-% SE)) +
  guides(fill = guide_legend(title = "Bin"), color = guide_legend(title = "Bin")) +
  geom_errorbar(aes(ymin=VarperSNP - SEperSNP, ymax=VarperSNP + SEperSNP), color="black", width=.2, position=position_dodge(1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_x_discrete(labels = tplotb48$Bin) +
  coord_flip()

plot_03_enrichment <- ggarrange(ldms_enriched_plot_a, ldms_enriched_plot_b, common.legend = F, nrow = 1, ncol = 2)


######  Save figures  ######

ggsave("./Fig01_Arrays.eps", plot_01_arrays, width = 180, height = 100, dpi = 300, units = "mm",scale=1.3, device="eps")
ggsave("./Fig02_WGS.eps", plot_02_WGS, width = 180, height = 160, dpi = 300, units = "mm",scale=1.3, device="eps")
ggsave("./Fig03_Enrichment.eps", plot_03_enrichment, width = 180, height = 100, dpi = 300, units = "mm",scale=1.3, device="eps")
