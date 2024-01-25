###################################################################
# R codes for comparing methods for selecting the number of ICs
# 04: make PRMSE figures
###################################################################
rm(list = ls())

library(ggplot2)
library(ggpubr)

load("03_PRMSE_misQ.RData")

#################################################################################################
# For S estimates
S_PMSE = data.frame(PMSE=c(S_PMSE1$sparse2,S_PMSE1$sparse3,S_PMSE1$sparse4,S_PMSE1$fast2,S_PMSE1$fast3,S_PMSE1$fast4,
                           S_PMSE2$sparse2,S_PMSE2$sparse3,S_PMSE2$sparse4,S_PMSE2$fast2,S_PMSE2$fast3,S_PMSE2$fast4,
                           S_PMSE3$sparse2,S_PMSE3$sparse3,S_PMSE3$sparse4,S_PMSE3$fast2,S_PMSE3$fast3,S_PMSE3$fast4),
                    Number=c(rep("Q=2",length(S_PMSE1$sparse2)),rep("Q=3",length(S_PMSE1$sparse3)),rep("Q=4",length(S_PMSE1$sparse4)),
                             rep("Q=2",length(S_PMSE1$fast2)),rep("Q=3",length(S_PMSE1$fast3)),rep("Q=4",length(S_PMSE1$fast4)),
                             rep("Q=2",length(S_PMSE2$sparse2)),rep("Q=3",length(S_PMSE2$sparse3)),rep("Q=4",length(S_PMSE2$sparse4)),
                             rep("Q=2",length(S_PMSE2$fast2)),rep("Q=3",length(S_PMSE2$fast3)),rep("Q=4",length(S_PMSE2$fast4)),
                             rep("Q=2",length(S_PMSE3$sparse2)),rep("Q=3",length(S_PMSE3$sparse3)),rep("Q=4",length(S_PMSE3$sparse4)),
                             rep("Q=2",length(S_PMSE3$fast2)),rep("Q=3",length(S_PMSE3$fast3)),rep("Q=4",length(S_PMSE3$fast4))),
                    SNR=c(rep("Low",6*length(S_PMSE1$sparse2)),
                          rep("Medium",6*length(S_PMSE1$sparse2)),
                          rep("High",6*length(S_PMSE1$sparse2))),
                    Method=c(rep("Sparse ICA",3*length(S_PMSE1$sparse2)),rep("Fast ICA",3*length(S_PMSE1$fast2)),
                             rep("Sparse ICA",3*length(S_PMSE2$sparse2)),rep("Fast ICA",3*length(S_PMSE2$fast2)),
                             rep("Sparse ICA",3*length(S_PMSE3$sparse2)),rep("Fast ICA",3*length(S_PMSE3$fast2)))
)
S_PMSE$SNR = factor(S_PMSE$SNR, levels = c("Low", "Medium", "High"))

p1 = ggplot(S_PMSE, aes(x = Number, y = PMSE, fill = Method)) + 
  geom_boxplot() +
  facet_wrap(~SNR) +
  labs(x = "Number of Components", y = "PRMSE")+
  theme_bw()


# For M estimates
M_PMSE = data.frame(PMSE=c(M_PMSE1$sparse2,M_PMSE1$sparse3,M_PMSE1$sparse4,M_PMSE1$fast2,M_PMSE1$fast3,M_PMSE1$fast4,
                           M_PMSE2$sparse2,M_PMSE2$sparse3,M_PMSE2$sparse4,M_PMSE2$fast2,M_PMSE2$fast3,M_PMSE2$fast4,
                           M_PMSE3$sparse2,M_PMSE3$sparse3,M_PMSE3$sparse4,M_PMSE3$fast2,M_PMSE3$fast3,M_PMSE3$fast4),
                    Number=c(rep("Q=2",length(M_PMSE1$sparse2)),rep("Q=3",length(M_PMSE1$sparse3)),rep("Q=4",length(M_PMSE1$sparse4)),
                             rep("Q=2",length(M_PMSE1$fast2)),rep("Q=3",length(M_PMSE1$fast3)),rep("Q=4",length(M_PMSE1$fast4)),
                             rep("Q=2",length(M_PMSE2$sparse2)),rep("Q=3",length(M_PMSE2$sparse3)),rep("Q=4",length(M_PMSE2$sparse4)),
                             rep("Q=2",length(M_PMSE2$fast2)),rep("Q=3",length(M_PMSE2$fast3)),rep("Q=4",length(M_PMSE2$fast4)),
                             rep("Q=2",length(M_PMSE3$sparse2)),rep("Q=3",length(M_PMSE3$sparse3)),rep("Q=4",length(M_PMSE3$sparse4)),
                             rep("Q=2",length(M_PMSE3$fast2)),rep("Q=3",length(M_PMSE3$fast3)),rep("Q=4",length(M_PMSE3$fast4))),
                    SNR=c(rep("Low",6*length(M_PMSE1$sparse2)),
                          rep("Medium",6*length(M_PMSE1$sparse2)),
                          rep("High",6*length(M_PMSE1$sparse2))),
                    Method=c(rep("Sparse ICA",3*length(M_PMSE1$sparse2)),rep("Fast ICA",3*length(M_PMSE1$fast2)),
                             rep("Sparse ICA",3*length(M_PMSE2$sparse2)),rep("Fast ICA",3*length(M_PMSE2$fast2)),
                             rep("Sparse ICA",3*length(M_PMSE3$sparse2)),rep("Fast ICA",3*length(M_PMSE3$fast2)))
)
M_PMSE$SNR = factor(M_PMSE$SNR, levels = c("Low", "Medium", "High"))

p2 = ggplot(M_PMSE, aes(x = Number, y = PMSE, fill = Method)) + 
  geom_boxplot() +
  facet_wrap(~SNR) +
  labs(x = "Number of Components", y = "PRMSE")+
  theme_bw()


ggarrange(p1, p2,
          labels = c("A", "B"),
          ncol = 2, nrow = 1,
          common.legend = T,legend = "bottom")
ggsave("Supp_Figure11.pdf",width = 13, height = 7)
