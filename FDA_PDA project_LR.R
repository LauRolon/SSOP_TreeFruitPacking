# FDA/PDA Project #
#By Laura Rolon
#Last updated: 9/25/2023 MLR

#Attach libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(psych)
library(svglite)
library(agricolae)

#Set working directory
setwd("G:/My Drive/Penn State/Research/File for R/FDA_PDA/TwoYears") #Change as needed

#### FDA project ####
#Upload data
micro_data<-read_excel("Results_TwoYears.xlsx", sheet="Summary", col_names = TRUE)
micro_data$logList<-as.numeric(micro_data$logList)


#### FDA: Aerobic plate counts ####
#ANOVA - Four way ANOVA with interactions
anova_APC<-aov(logAPC ~ Treatment:Time:Year, data=micro_data)
summary(anova_APC) # significant 3 way interaction term p = 2e-16

# #tukey test
# tukey_APC<-HSD.test(anova_APC, trt="F_T_t_Y") 
# 
# #Make dataframe with Tukey groups
# tukey_APC_groups<-tukey_APC$groups
# tukey_APC_groups$Sample<-as.character(rownames(tukey_APC_groups))
# tukey_APC_groups_order<-tukey_APC_groups[with(tukey_APC_groups, order(Sample)),]
# 
# # Summary plots 
# #Calculate mean and sd by treatment
# 
# APC_stat<-describeBy(micro_data$logAPC, list(micro_data$F_T_t_Y), mat = TRUE) #By facility and treatment
# APC_stat$Facility<-c(rep("F1", 16), rep("F2",16), rep("F3",16))
# APC_stat$Treatment<-rep(c(rep("T1", 4), rep("T2",4), rep("T3",4), rep("T4", 4)),3)
# APC_stat$Time<-rep(c(rep("After",2), rep("Before",2)),12)
# APC_stat$Year<-rep(c("Y1","Y2"),24)
# APC_stat$Time <- factor(APC_stat$Time, levels = rev(levels(factor(APC_stat$Time))))
# APC_stat$Tukey<-tukey_APC_groups_order$groups
# 
# 
# APC_stat_plot<-ggplot(APC_stat, aes(x=Time, y=mean, fill=Time))+
#   geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
#   geom_text(aes(label=Tukey), position=position_dodge(width=0.9), vjust=-0.25)+
#   facet_grid(Facility~Year+Treatment)+ theme(legend.position = 'bottom')+
#   theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
#   ylab("log CFU/swab") + xlab("")+
#   theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
#   theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
#   theme(axis.title = element_text(size=15,color='black')) +
#   theme(panel.background = element_rect(fill='white', color = NA),
#         plot.background = element_rect(fill = 'white',color = NA), 
#         panel.border = element_rect(color='black',fill = NA,size=1))+
#   theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
#   theme(strip.background= element_blank(), strip.text = element_text(size=15),
#         panel.border = element_rect(color="black", fill=NA))+
#   scale_y_continuous(breaks=c(0,2,4,6,8,10,12), limits = c(0,12))+
#   ggtitle("Aerobic plate count - Summary")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
# APC_stat_plot
# ggsave("SummaryAPC.svg", plot=APC_stat_plot, device="svg", width=8, height=11, units="in", dpi=600)

#Seems like variance is too big between the two years. Split data by year and re analyze

micro_data_Y1<-subset(micro_data, Year=="Y1")
micro_data_Y2<-subset(micro_data, Year=="Y2")

#ANOVA
anova_APC_Y1<-aov(logAPC ~ F_T_t_Y, data=micro_data_Y1)
summary(anova_APC_Y1)

anova_APC_Y2<-aov(logAPC ~ F_T_t_Y, data=micro_data_Y2)
summary(anova_APC_Y2)

#tukey test
tukey_APC_Y1<-HSD.test(anova_APC_Y1, trt="F_T_t_Y") 
tukey_APC_Y2<-HSD.test(anova_APC_Y2, trt="F_T_t_Y") 

#Make dataframe with Tukey groups
tukey_APC_groups_Y1<-tukey_APC_Y1$groups
tukey_APC_groups_Y1$Sample<-as.character(rownames(tukey_APC_groups_Y1))
tukey_APC_groups_order_Y1<-tukey_APC_groups_Y1[with(tukey_APC_groups_Y1, order(Sample)),]

tukey_APC_groups_Y2<-tukey_APC_Y2$groups
tukey_APC_groups_Y2$Sample<-as.character(rownames(tukey_APC_groups_Y2))
tukey_APC_groups_order_Y2<-tukey_APC_groups_Y2[with(tukey_APC_groups_Y2, order(Sample)),]

# Summary plots 
#Calculate mean and sd by treatment
APC_stat_Y1<-describeBy(micro_data_Y1$logAPC, list(micro_data_Y1$F_T_t_Y), mat = TRUE) #By facility and treatment
APC_stat_Y1$Facility<-c(rep("F1", 8), rep("F2",8), rep("F3",8))
APC_stat_Y1$Treatment<-rep(c(rep("T1", 2), rep("T2",2), rep("T3",2), rep("T4", 2)),3)
APC_stat_Y1$Time<-rep(c("After","Before"),12)
APC_stat_Y1$Time <- factor(APC_stat_Y1$Time, levels = rev(levels(factor(APC_stat_Y1$Time))))
APC_stat_Y1$Tukey<-tukey_APC_groups_order_Y1$groups

APC_stat_Y2<-describeBy(micro_data_Y2$logAPC, list(micro_data_Y2$F_T_t_Y), mat = TRUE) #By facility and treatment
APC_stat_Y2$Facility<-c(rep("F1", 8), rep("F2",8), rep("F3",8))
APC_stat_Y2$Treatment<-rep(c(rep("T1", 2), rep("T2",2), rep("T3",2), rep("T4", 2)),3)
APC_stat_Y2$Time<-rep(c("After","Before"),12)
APC_stat_Y2$Time <- factor(APC_stat_Y2$Time, levels = rev(levels(factor(APC_stat_Y2$Time))))
APC_stat_Y2$Tukey<-tukey_APC_groups_order_Y2$groups

#Plot
APC_stat_plot_Y1<-ggplot(APC_stat_Y1, aes(x=Time, y=mean, fill=Time))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=Tukey), position=position_dodge(width=0.9), vjust=-0.25)+
  facet_grid(Facility~Treatment)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log CFU/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12), limits = c(0,12))+
  ggtitle("Aerobic plate count - Summary Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
APC_stat_plot_Y1
ggsave("SummaryAPC_Y1.svg", plot=APC_stat_plot_Y1, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryAPC_Y1.png", plot=APC_stat_plot_Y1, device="png", width=8, height=11, units="in", dpi=600)

APC_stat_plot_Y2<-ggplot(APC_stat_Y2, aes(x=Time, y=mean, fill=Time))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=Tukey), position=position_dodge(width=0.9), vjust=-0.25)+
  facet_grid(Facility~Treatment)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log CFU/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12), limits = c(0,12))+
  ggtitle("Aerobic plate count - Summary Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
APC_stat_plot_Y2
ggsave("SummaryAPC_Y2.svg", plot=APC_stat_plot_Y2, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryAPC_Y2.png", plot=APC_stat_plot_Y2, device="png", width=8, height=11, units="in", dpi=600)



# If the question is which treatment was most effective at reducing APC overall, subset by treatment and look at each independently
micro_data_Y1T1<-subset(micro_data, Treatment=="T1" & Year =="Y1")
micro_data_Y1T2<-subset(micro_data, Treatment=="T2" & Year =="Y1")
micro_data_Y1T3<-subset(micro_data, Treatment=="T3" & Year =="Y1")
micro_data_Y1T4<-subset(micro_data, Treatment=="T4" & Year =="Y1")

micro_data_Y2T1<-subset(micro_data, Treatment=="T1" & Year =="Y2")
micro_data_Y2T2<-subset(micro_data, Treatment=="T2" & Year =="Y2")
micro_data_Y2T3<-subset(micro_data, Treatment=="T3" & Year =="Y2")
micro_data_Y2T4<-subset(micro_data, Treatment=="T4" & Year =="Y2")

#Paired t-test for treatment effect regardless of facility
#Y1
t.test(micro_data_Y1T1$logAPC ~ micro_data_Y1T1$Time, paired=TRUE) # p=0.01349
t.test(micro_data_Y1T2$logAPC ~ micro_data_Y1T2$Time, paired=TRUE) # p=0.0008933
t.test(micro_data_Y1T3$logAPC ~ micro_data_Y1T3$Time, paired=TRUE) # p=4.902e-07
t.test(micro_data_Y1T4$logAPC ~ micro_data_Y1T4$Time, paired=TRUE) # p=0.0008139

#Y2
t.test(micro_data_Y2T1$logAPC ~ micro_data_Y2T1$Time, paired=TRUE) # p=0.071
t.test(micro_data_Y2T2$logAPC ~ micro_data_Y2T2$Time, paired=TRUE) # p=0.002476
t.test(micro_data_Y2T3$logAPC ~ micro_data_Y2T3$Time, paired=TRUE) # p=7.141e-06
t.test(micro_data_Y2T4$logAPC ~ micro_data_Y2T4$Time, paired=TRUE) # p=1.853e-06


# Summary plots 
#Calculate mean and sd by treatment
APC_stat<-describeBy(micro_data$logAPC, list(micro_data$Year,micro_data$Trt_time ), mat = TRUE) #By facility and treatment
APC_stat$Treatment<-c(rep("T1", 4), rep("T2",4), rep("T3",4), rep("T4", 4))
APC_stat$Time<-rep(c(rep("After",2),rep("Before",2)),4)
APC_stat$Year<-rep(c("Y1","Y2"),8)
APC_stat$Time <- factor(APC_stat$Time, levels = rev(levels(factor(APC_stat$Time))))


#Plot
APC_stat_plot<-ggplot(APC_stat, aes(x=Time, y=mean, fill=Time))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  facet_grid(Year~Treatment)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log CFU/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  #theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12), limits = c(0,12))+
  ggtitle("Aerobic plate count - Summary")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
APC_stat_plot
ggsave("SummaryAPC.svg", plot=APC_stat_plot, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryAPC.png", plot=APC_stat_plot, device="png", width=8, height=11, units="in", dpi=600)




#Subset data by treatment and year and evaluate overall effect of treatments on each year
micro_data_Y1<-subset(micro_data, Year=="Y1")
micro_data_Y2<-subset(micro_data, Year=="Y2")

#ANOVA
anova_APC_Y1<-aov(logAPC ~ F_T_t_Y, data=micro_data_Y1)
summary(anova_APC_Y1)

anova_APC_Y2<-aov(logAPC ~ F_T_t_Y, data=micro_data_Y2)
summary(anova_APC_Y2)

#tukey test
tukey_APC_Y1<-HSD.test(anova_APC_Y1, trt="F_T_t_Y") 
tukey_APC_Y2<-HSD.test(anova_APC_Y2, trt="F_T_t_Y") 

#Make dataframe with Tukey groups
tukey_APC_groups_Y1<-tukey_APC_Y1$groups
tukey_APC_groups_Y1$Sample<-as.character(rownames(tukey_APC_groups_Y1))
tukey_APC_groups_order_Y1<-tukey_APC_groups_Y1[with(tukey_APC_groups_Y1, order(Sample)),]

tukey_APC_groups_Y2<-tukey_APC_Y2$groups
tukey_APC_groups_Y2$Sample<-as.character(rownames(tukey_APC_groups_Y2))
tukey_APC_groups_order_Y2<-tukey_APC_groups_Y2[with(tukey_APC_groups_Y2, order(Sample)),]

# Summary plots 
#Calculate mean and sd by treatment
APC_stat_Y1<-describeBy(micro_data_Y1$logAPC, list(micro_data_Y1$F_T_t_Y), mat = TRUE) #By facility and treatment
APC_stat_Y1$Facility<-c(rep("F1", 8), rep("F2",8), rep("F3",8))
APC_stat_Y1$Treatment<-rep(c(rep("T1", 2), rep("T2",2), rep("T3",2), rep("T4", 2)),3)
APC_stat_Y1$Time<-rep(c("After","Before"),12)
APC_stat_Y1$Time <- factor(APC_stat_Y1$Time, levels = rev(levels(factor(APC_stat_Y1$Time))))
APC_stat_Y1$Tukey<-tukey_APC_groups_order_Y1$groups

APC_stat_Y2<-describeBy(micro_data_Y2$logAPC, list(micro_data_Y2$F_T_t_Y), mat = TRUE) #By facility and treatment
APC_stat_Y2$Facility<-c(rep("F1", 8), rep("F2",8), rep("F3",8))
APC_stat_Y2$Treatment<-rep(c(rep("T1", 2), rep("T2",2), rep("T3",2), rep("T4", 2)),3)
APC_stat_Y2$Time<-rep(c("After","Before"),12)
APC_stat_Y2$Time <- factor(APC_stat_Y2$Time, levels = rev(levels(factor(APC_stat_Y2$Time))))
APC_stat_Y2$Tukey<-tukey_APC_groups_order_Y2$groups

#Plot
APC_stat_plot_Y1<-ggplot(APC_stat_Y1, aes(x=Time, y=mean, fill=Time))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=Tukey), position=position_dodge(width=0.9), vjust=-0.25)+
  facet_grid(Facility~Treatment)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log CFU/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12), limits = c(0,12))+
  ggtitle("Aerobic plate count - Summary Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
APC_stat_plot_Y1
ggsave("SummaryAPC_Y1.svg", plot=APC_stat_plot_Y1, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryAPC_Y1.png", plot=APC_stat_plot_Y1, device="png", width=8, height=11, units="in", dpi=600)

APC_stat_plot_Y2<-ggplot(APC_stat_Y2, aes(x=Time, y=mean, fill=Time))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=Tukey), position=position_dodge(width=0.9), vjust=-0.25)+
  facet_grid(Facility~Treatment)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log CFU/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12), limits = c(0,12))+
  ggtitle("Aerobic plate count - Summary Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
APC_stat_plot_Y2
ggsave("SummaryAPC_Y2.svg", plot=APC_stat_plot_Y2, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryAPC_Y2.png", plot=APC_stat_plot_Y2, device="png", width=8, height=11, units="in", dpi=600)


#Paired t-test for each Treatment, facility and year
#Make Time column as a factor
micro_data_Y1$Time<-as.factor(micro_data_Y1$Time)
micro_data_Y2$Time<-as.factor(micro_data_Y2$Time)

#Subset data by facility and treatment
F1T1Y1<-subset(micro_data_Y1, Facility=="F1" & Treatment=="T1")
F1T2Y1<-subset(micro_data_Y1, Facility=="F1" & Treatment=="T2")
F1T3Y1<-subset(micro_data_Y1, Facility=="F1" & Treatment=="T3")
F1T4Y1<-subset(micro_data_Y1, Facility=="F1" & Treatment=="T4")

F2T1Y1<-subset(micro_data_Y1, Facility=="F2" & Treatment=="T1")
F2T2Y1<-subset(micro_data_Y1, Facility=="F2" & Treatment=="T2")
F2T3Y1<-subset(micro_data_Y1, Facility=="F2" & Treatment=="T3")
F2T4Y1<-subset(micro_data_Y1, Facility=="F2" & Treatment=="T4")

F3T1Y1<-subset(micro_data_Y1, Facility=="F3" & Treatment=="T1")
F3T2Y1<-subset(micro_data_Y1, Facility=="F3" & Treatment=="T2")
F3T3Y1<-subset(micro_data_Y1, Facility=="F3" & Treatment=="T3")
F3T4Y1<-subset(micro_data_Y1, Facility=="F3" & Treatment=="T4")

F1T1Y2<-subset(micro_data_Y2, Facility=="F1" & Treatment=="T1")
F1T2Y2<-subset(micro_data_Y2, Facility=="F1" & Treatment=="T2")
F1T3Y2<-subset(micro_data_Y2, Facility=="F1" & Treatment=="T3")
F1T4Y2<-subset(micro_data_Y2, Facility=="F1" & Treatment=="T4")

F2T1Y2<-subset(micro_data_Y2, Facility=="F2" & Treatment=="T1")
F2T2Y2<-subset(micro_data_Y2, Facility=="F2" & Treatment=="T2")
F2T3Y2<-subset(micro_data_Y2, Facility=="F2" & Treatment=="T3")
F2T4Y2<-subset(micro_data_Y2, Facility=="F2" & Treatment=="T4")

F3T1Y2<-subset(micro_data_Y2, Facility=="F3" & Treatment=="T1")
F3T2Y2<-subset(micro_data_Y2, Facility=="F3" & Treatment=="T2")
F3T3Y2<-subset(micro_data_Y2, Facility=="F3" & Treatment=="T3")
F3T4Y2<-subset(micro_data_Y2, Facility=="F3" & Treatment=="T4")

#Paired t-test by facility and treatment and year
t.test(F1T1Y1$logAPC ~ F1T1Y1$Time, paired=TRUE) #p = 0.4733
t.test(F1T2Y1$logAPC ~ F1T2Y1$Time, paired=TRUE) #p-value = 0.06483
t.test(F1T3Y1$logAPC ~ F1T3Y1$Time, paired=TRUE) #p-value = 0.0005237
t.test(F1T4Y1$logAPC ~ F1T4Y1$Time, paired=TRUE) # p-value = 0.1289

t.test(F2T1Y1$logAPC ~ F2T1Y1$Time, paired=TRUE) #p-value = 0.7086
t.test(F2T2Y1$logAPC ~ F2T2Y1$Time, paired=TRUE) #p-value = 0.0163
t.test(F2T3Y1$logAPC ~ F2T3Y1$Time, paired=TRUE) #p-value = 0.002988
t.test(F2T4Y1$logAPC ~ F2T4Y1$Time, paired=TRUE) # p-value = 0.1436

t.test(F3T1Y1$logAPC ~ F3T1Y1$Time, paired=TRUE) #p-value = 0.001395
t.test(F3T2Y1$logAPC ~ F3T2Y1$Time, paired=TRUE) #p-value = 0.0636
t.test(F3T3Y1$logAPC ~ F3T3Y1$Time, paired=TRUE) #p-value = 0.005812
t.test(F3T4Y1$logAPC ~ F3T4Y1$Time, paired=TRUE) # p-value = 0.0001817

t.test(F1T1Y2$logAPC ~ F1T1Y2$Time, paired=TRUE) #p-value = 0.6554
t.test(F1T2Y2$logAPC ~ F1T2Y2$Time, paired=TRUE) #p-value = 0.00425
t.test(F1T3Y2$logAPC ~ F1T3Y2$Time, paired=TRUE) #p-value = 0.003829
t.test(F1T4Y2$logAPC ~ F1T4Y2$Time, paired=TRUE) #p-value = 0.005193

t.test(F2T1Y2$logAPC ~ F2T1Y2$Time, paired=TRUE) #p-value = 0.01722
t.test(F2T2Y2$logAPC ~ F2T2Y2$Time, paired=TRUE) #p-value = 0.1243
t.test(F2T3Y2$logAPC ~ F2T3Y2$Time, paired=TRUE) #p-value = 0.001366
t.test(F2T4Y2$logAPC ~ F2T4Y2$Time, paired=TRUE) #p-value = 0.001127

t.test(F3T1Y2$logAPC ~ F3T1Y2$Time, paired=TRUE) #p-value = 0.8642
t.test(F3T2Y2$logAPC ~ F3T2Y2$Time, paired=TRUE) #p-value = 0.5528
t.test(F3T3Y2$logAPC ~ F3T3Y2$Time, paired=TRUE) #p-value = 0.05449
t.test(F3T4Y2$logAPC ~ F3T4Y2$Time, paired=TRUE) #p-value = 0.02876


#### FDA: Listeria spp. MPN ####
#Working with data separated by year
#In Y1 data, subset to remove F2 that doesn't have data for Listeria spp. quantification
ListMPN_Y1<-subset(micro_data_Y1, Facility != "F2")
ListMPN_Y1$logList<-as.numeric(ListMPN_Y1$logList)
micro_data_Y2$logList<-as.numeric(micro_data_Y2$logList)

#ANOVA
anova_List_Y1<-aov(logList ~ F_T_t_Y, data=ListMPN_Y1)
summary(anova_List_Y1) #p=0.0485

anova_List_Y2<-aov(logList ~ F_T_t_Y, data=micro_data_Y2)
summary(anova_List_Y2) #p=9.38e-07 ***


#tukey test
tukey_List_Y1<-HSD.test(anova_List_Y1, trt="F_T_t_Y") 
tukey_List_Y2<-HSD.test(anova_List_Y2, trt="F_T_t_Y") 

#Make dataframe with Tukey groups
tukey_List_groups_Y1<-tukey_List_Y1$groups
tukey_List_groups_Y1$Sample<-as.character(rownames(tukey_List_groups_Y1))
tukey_List_groups_order_Y1<-tukey_List_groups_Y1[with(tukey_List_groups_Y1, order(Sample)),]

tukey_List_groups_Y2<-tukey_List_Y2$groups
tukey_List_groups_Y2$Sample<-as.character(rownames(tukey_List_groups_Y2))
tukey_List_groups_order_Y2<-tukey_List_groups_Y2[with(tukey_List_groups_Y2, order(Sample)),]

# Summary plots 
#Calculate mean and sd by treatment
List_stat_Y1<-describeBy(ListMPN_Y1$logList, list(ListMPN_Y1$F_T_t_Y), mat = TRUE) #By facility and treatment
List_stat_Y1$Facility<-c(rep("F1", 8),rep("F3",8))
List_stat_Y1$Treatment<-rep(c(rep("T1", 2), rep("T2",2), rep("T3",2), rep("T4", 2)),2)
List_stat_Y1$Time<-rep(c("After","Before"),8)
List_stat_Y1$Time <- factor(List_stat_Y1$Time, levels = rev(levels(factor(List_stat_Y1$Time))))
List_stat_Y1$Tukey<-tukey_List_groups_order_Y1$groups

List_stat_Y2<-describeBy(micro_data_Y2$logList, list(micro_data_Y2$F_T_t_Y), mat = TRUE) #By facility and treatment
List_stat_Y2$Facility<-c(rep("F1", 8), rep("F2",8), rep("F3",8))
List_stat_Y2$Treatment<-rep(c(rep("T1", 2), rep("T2",2), rep("T3",2), rep("T4", 2)),3)
List_stat_Y2$Time<-rep(c("After","Before"),12)
List_stat_Y2$Time <- factor(List_stat_Y2$Time, levels = rev(levels(factor(List_stat_Y2$Time))))
List_stat_Y2$Tukey<-tukey_List_groups_order_Y2$groups

#Plot
List_stat_plot_Y1<-ggplot(List_stat_Y1, aes(x=Time, y=mean, fill=Time))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=Tukey), position=position_dodge(width=0.9), vjust=-0.25)+
  facet_grid(Facility~Treatment)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log MPN/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits = c(0,10))+
  ggtitle("Listeria spp MPN - Summary Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
List_stat_plot_Y1
ggsave("SummaryList_Y1.svg", plot=List_stat_plot_Y1, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryList_Y1.png", plot=List_stat_plot_Y1, device="png", width=8, height=11, units="in", dpi=600)

List_stat_plot_Y2<-ggplot(List_stat_Y2, aes(x=Time, y=mean, fill=Time))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=Tukey), position=position_dodge(width=0.9), vjust=-0.25)+
  facet_grid(Facility~Treatment)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log MPN/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits = c(0,10))+
  ggtitle("Listeria spp. MPN - Summary Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
List_stat_plot_Y2
ggsave("SummaryList_Y2.svg", plot=List_stat_plot_Y2, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryList_Y2.png", plot=List_stat_plot_Y2, device="png", width=8, height=11, units="in", dpi=600)

#No significant difference in before and after samples due to application of SSOPs

#Paired t-test for each Treatment, facility and year
#Make Time column as a factor
ListMPN_Y1$Time<-as.factor(micro_data_Y1$Time)

#Subset data by facility and treatment
F1T1Y1_List<-subset(ListMPN_Y1, Facility=="F1" & Treatment=="T1")
F1T2Y1_List<-subset(ListMPN_Y1, Facility=="F1" & Treatment=="T2")
F1T3Y1_List<-subset(ListMPN_Y1, Facility=="F1" & Treatment=="T3")
F1T4Y1_List<-subset(ListMPN_Y1, Facility=="F1" & Treatment=="T4")

F3T1Y1_List<-subset(ListMPN_Y1, Facility=="F3" & Treatment=="T1")
F3T2Y1_List<-subset(ListMPN_Y1, Facility=="F3" & Treatment=="T2")
F3T3Y1_List<-subset(ListMPN_Y1, Facility=="F3" & Treatment=="T3")
F3T4Y1_List<-subset(ListMPN_Y1, Facility=="F3" & Treatment=="T4")


#Paired t-test by facility and treatment and year
t.test(F1T1Y1_List$logList ~ F1T1Y1_List$Time, paired=TRUE) #p-value = 0.8928
t.test(F1T2Y1_List$logList ~ F1T2Y1_List$Time, paired=TRUE) # p-value = 0.1805
t.test(F1T3Y1_List$logList ~ F1T3Y1_List$Time, paired=TRUE) #p-value = 0.713
t.test(F1T4Y1_List$logList ~ F1T4Y1_List$Time, paired=TRUE) # p-value = 0.1768

t.test(F3T1Y1_List$logList ~ F3T1Y1_List$Time, paired=TRUE) #p-value = 0.6253
t.test(F3T2Y1_List$logList ~ F3T2Y1_List$Time, paired=TRUE) #p-value = 0.1113
t.test(F3T3Y1_List$logList ~ F3T3Y1_List$Time, paired=TRUE) # p-value = 0.002764
t.test(F3T4Y1_List$logList ~ F3T4Y1_List$Time, paired=TRUE) # p-value = 0.04099

t.test(F1T1Y2$logList ~ F1T1Y2$Time, paired=TRUE) #p-value = 0.06103
t.test(F1T2Y2$logList ~ F1T2Y2$Time, paired=TRUE) #p-value = 0.06602
t.test(F1T3Y2$logList ~ F1T3Y2$Time, paired=TRUE) # p-value = 0.02566
t.test(F1T4Y2$logList ~ F1T4Y2$Time, paired=TRUE) #p-value = 0.08736

t.test(F2T1Y2$logList ~ F2T1Y2$Time, paired=TRUE) #p-value = 0.9656
t.test(F2T2Y2$logList ~ F2T2Y2$Time, paired=TRUE) #p-value = 0.7814
t.test(F2T3Y2$logList ~ F2T3Y2$Time, paired=TRUE) # p-value = 0.5067
t.test(F2T4Y2$logList ~ F2T4Y2$Time, paired=TRUE) #p-value = 0.09947

t.test(F3T1Y2$logList ~ F3T1Y2$Time, paired=TRUE) #p-value = 0.09266
t.test(F3T2Y2$logList ~ F3T2Y2$Time, paired=TRUE) #p-value = 0.9341
t.test(F3T3Y2$logList ~ F3T3Y2$Time, paired=TRUE) #p-value = 0.6796
t.test(F3T4Y2$logList ~ F3T4Y2$Time, paired=TRUE) #p-value = 0.0223


# If the question is which treatment was most effective at reducing List overall, subset by treatment and look at each independently
#Paired t-test for treatment effect regardless of facility


#Y1
t.test(micro_data_Y1T1$logList ~ micro_data_Y1T1$Time, paired=TRUE) # p=0.6795
t.test(micro_data_Y1T2$logList ~ micro_data_Y1T2$Time, paired=TRUE) # p=0.02411
t.test(micro_data_Y1T3$logList ~ micro_data_Y1T3$Time, paired=TRUE) # p=0.03695
t.test(micro_data_Y1T4$logList ~ micro_data_Y1T4$Time, paired=TRUE) # p=0.01576

#Y2
t.test(micro_data_Y2T1$logList ~ micro_data_Y2T1$Time, paired=TRUE) # p=0.7633
t.test(micro_data_Y2T2$logList ~ micro_data_Y2T2$Time, paired=TRUE) # p=0.1513
t.test(micro_data_Y2T3$logList ~ micro_data_Y2T3$Time, paired=TRUE) # p=0.4365
t.test(micro_data_Y2T4$logList ~ micro_data_Y2T4$Time, paired=TRUE) # p=0.0005468


# Summary plots 

#Calculate mean and sd by treatment
List_stat<-describeBy(micro_data$logList, list(micro_data$Year,micro_data$Trt_time ), mat = TRUE) #By facility and treatment
List_stat$Treatment<-c(rep("T1", 4), rep("T2",4), rep("T3",4), rep("T4", 4))
List_stat$Time<-rep(c(rep("After",2),rep("Before",2)),4)
List_stat$Year<-rep(c("Y1","Y2"),8)
List_stat$Time <- factor(List_stat$Time, levels = rev(levels(factor(List_stat$Time))))


#Plot
List_stat_plot<-ggplot(List_stat, aes(x=Time, y=mean, fill=Time))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  facet_grid(Year~Treatment)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log MPN/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  #theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12), limits = c(0,12))+
  ggtitle("Listeria spp - Summary")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
List_stat_plot
ggsave("SummaryList.svg", plot=List_stat_plot, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryList.png", plot=List_stat_plot, device="png", width=8, height=11, units="in", dpi=600)


#### FDA: Lm MPN ####

#ANOVA
anova_Lm_Y1<-aov(logLm ~ F_T_t_Y, data=micro_data_Y1)
summary(anova_Lm_Y1)
#p-value  0.0152 *

anova_Lm_Y2<-aov(logLm ~ F_T_t_Y, data=micro_data_Y2)
summary(anova_Lm_Y2)
#p-value 0.131 - No Tukey Test needed


#tukey test
tukey_Lm_Y1<-HSD.test(anova_Lm_Y1, trt="F_T_t_Y") #No sig different groups

# Summary plots
#Calculate mean and sd by treatment
Lm_stat_Y1<-describeBy(micro_data_Y1$logLm, list(micro_data_Y1$F_T_t_Y), mat = TRUE) #By facility and treatment
Lm_stat_Y1$Facility<-c(rep("F1", 8), rep("F2",8), rep("F3",8))
Lm_stat_Y1$Treatment<-rep(c(rep("T1", 2), rep("T2",2), rep("T3",2), rep("T4", 2)),3)
Lm_stat_Y1$Time<-rep(c("After","Before"),12)
Lm_stat_Y1$Time <- factor(Lm_stat_Y1$Time, levels = rev(levels(factor(Lm_stat_Y1$Time))))


Lm_stat_Y1_plot<-ggplot(Lm_stat_Y1, aes(x=Time, y=mean, fill=Time))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  #geom_hline(yintercept=2.55, linetype=2, color="grey20")+
  facet_grid(Facility~Treatment)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log MPN/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits = c(0,10))+
  ggtitle("L. monocytogenes MPN Y1- Summary")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
Lm_stat_Y1_plot
ggsave("SummaryLm Y1.svg", plot=Lm_stat_Y1_plot, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryLm Y1.png", plot=Lm_stat_Y1_plot, device="png", width=8, height=11, units="in", dpi=600)

Lm_stat_Y2<-describeBy(micro_data_Y2$logLm, list(micro_data_Y2$F_T_t_Y), mat = TRUE) #By facility and treatment
Lm_stat_Y2$Facility<-c(rep("F1", 8), rep("F2",8), rep("F3",8))
Lm_stat_Y2$Treatment<-rep(c(rep("T1", 2), rep("T2",2), rep("T3",2), rep("T4", 2)),3)
Lm_stat_Y2$Time<-rep(c("After","Before"),12)
Lm_stat_Y2$Time <- factor(Lm_stat_Y2$Time, levels = rev(levels(factor(Lm_stat_Y2$Time))))


Lm_stat_Y2_plot<-ggplot(Lm_stat_Y2, aes(x=Time, y=mean, fill=Time))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  #geom_hline(yintercept=2.55, linetype=2, color="grey20")+
  facet_grid(Facility~Treatment)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log MPN/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits = c(0,10))+
  ggtitle("L. monocytogenes MPN Y2- Summary")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
Lm_stat_Y2_plot
ggsave("SummaryLm Y2.svg", plot=Lm_stat_Y2_plot, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryLm Y2.png", plot=Lm_stat_Y2_plot, device="png", width=8, height=11, units="in", dpi=600)

#Paired t-test by facility and treatment and year
t.test(F1T1Y1$logLm ~ F1T1Y1$Time, paired=TRUE) #p-value = 0.4825
t.test(F1T2Y1$logLm ~ F1T2Y1$Time, paired=TRUE) # p-value = 0.0246
t.test(F1T3Y1$logLm ~ F1T3Y1$Time, paired=TRUE) #p-value = 0.2189
t.test(F1T4Y1$logLm ~ F1T4Y1$Time, paired=TRUE) # p-value = 0.08132

t.test(F2T1Y1$logLm ~ F2T1Y1$Time, paired=TRUE) #p-value = 0.6213
t.test(F2T2Y1$logLm ~ F2T2Y1$Time, paired=TRUE) #p-value = 0.5222
t.test(F2T3Y1$logLm ~ F2T3Y1$Time, paired=TRUE) #p-value = 0.3739
t.test(F2T4Y1$logLm ~ F2T4Y1$Time, paired=TRUE) #p-value = 0.3739

t.test(F3T1Y1$logLm ~ F3T1Y1$Time, paired=TRUE) #p-value = 0.7103
t.test(F3T2Y1$logLm ~ F3T2Y1$Time, paired=TRUE) #p-value = 0.01812
t.test(F3T3Y1$logLm ~ F3T3Y1$Time, paired=TRUE) #p-value = 0.01741
t.test(F3T4Y1$logLm ~ F3T4Y1$Time, paired=TRUE) #p-value = 0.02332

t.test(F1T1Y2$logLm ~ F1T1Y2$Time, paired=TRUE) # p-value = 0.9988
t.test(F1T2Y2$logLm ~ F1T2Y2$Time, paired=TRUE) #p-value = 0.4582
t.test(F1T3Y2$logLm ~ F1T3Y2$Time, paired=TRUE) #p-value = 0.1354
t.test(F1T4Y2$logLm ~ F1T4Y2$Time, paired=TRUE) #p-value = 0.2858

t.test(F2T1Y2$logLm ~ F2T1Y2$Time, paired=TRUE) #p-value = 0.6965
t.test(F2T2Y2$logLm ~ F2T2Y2$Time, paired=TRUE) #p-value = NA
t.test(F2T3Y2$logLm ~ F2T3Y2$Time, paired=TRUE) #p-value = NA
t.test(F2T4Y2$logLm ~ F2T4Y2$Time, paired=TRUE) #p-value = 0.6406

t.test(F3T1Y2$logLm ~ F3T1Y2$Time, paired=TRUE) #p-value = 0.1821
t.test(F3T2Y2$logLm ~ F3T2Y2$Time, paired=TRUE) #p-value = NA
t.test(F3T3Y2$logLm ~ F3T3Y2$Time, paired=TRUE) #p-value = 0.3739
t.test(F3T4Y2$logLm ~ F3T4Y2$Time, paired=TRUE) #p-value = 0.374


# If the question is which treatment was most effective at reducing Lm overall, subset by treatment and look at each independently
#Paired t-test for treatment effect regardless of facility


#Y1
t.test(micro_data_Y1T1$logLm ~ micro_data_Y1T1$Time, paired=TRUE) # p=0.7356
t.test(micro_data_Y1T2$logLm ~ micro_data_Y1T2$Time, paired=TRUE) # p=0.0448
t.test(micro_data_Y1T3$logLm ~ micro_data_Y1T3$Time, paired=TRUE) # p=0.1107
t.test(micro_data_Y1T4$logLm ~ micro_data_Y1T4$Time, paired=TRUE) # p=0.002682

#Y2
t.test(micro_data_Y2T1$logLm ~ micro_data_Y2T1$Time, paired=TRUE) # p=0.2901
t.test(micro_data_Y2T2$logLm ~ micro_data_Y2T2$Time, paired=TRUE) # p=0.415
t.test(micro_data_Y2T3$logLm ~ micro_data_Y2T3$Time, paired=TRUE) # p=0.08317
t.test(micro_data_Y2T4$logLm ~ micro_data_Y2T4$Time, paired=TRUE) # p=0.2503


# Summary plots 

#Calculate mean and sd by treatment
Lm_stat<-describeBy(micro_data$logLm, list(micro_data$Year,micro_data$Trt_time ), mat = TRUE) #By facility and treatment
Lm_stat$Treatment<-c(rep("T1", 4), rep("T2",4), rep("T3",4), rep("T4", 4))
Lm_stat$Time<-rep(c(rep("After",2),rep("Before",2)),4)
Lm_stat$Year<-rep(c("Y1","Y2"),8)
Lm_stat$Time <- factor(Lm_stat$Time, levels = rev(levels(factor(Lm_stat$Time))))


#Plot
Lm_stat_plot<-ggplot(Lm_stat, aes(x=Time, y=mean, fill=Time))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  facet_grid(Year~Treatment)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log MPN/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  #theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12), limits = c(0,12))+
  ggtitle("L. monocytogenes - Summary")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
Lm_stat_plot
ggsave("SummaryLm.svg", plot=Lm_stat_plot, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryLm.png", plot=Lm_stat_plot, device="png", width=8, height=11, units="in", dpi=600)

#### FDA: Lm presence by enrichment ####
#Upload data
Lm_enrichment<-read_excel("Results_TwoYears.xlsx", sheet=3, col_names = TRUE)

Lm_enrichment$Time <- factor(Lm_enrichment$Time, levels = rev(levels(factor(Lm_enrichment$Time))))

#All Facilities
Lm_presence_plot<-ggplot(Lm_enrichment, aes(x=Time, y=Percentage, fill=Lm))+
  geom_bar(stat='identity', color='black')+
  facet_grid(Facility~Year+Treatment)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  ylab("Percentage of samples") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_fill_manual(values = c('#F8CD9C','#EA7580'))+
  theme(legend.position = "bottom")+
  ggtitle("L. monocytogenes occurrance", subtitle = 'By enrichment method')+theme(plot.title = element_text(hjust = 0.5, face = 'italic'))
Lm_presence_plot
ggsave("Lm_presence.svg", plot=Lm_presence_plot, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("Lm_presence.png", plot=Lm_presence_plot, device="png", width=8, height=11, units="in", dpi=600)

#All Facilities combined
#Upload data
Lm_enrichment_nofac<-read_excel("Results_TwoYears.xlsx", sheet=5, col_names = TRUE)
Lm_enrichment_nofac$Time <- factor(Lm_enrichment_nofac$Time, levels = rev(levels(factor(Lm_enrichment_nofac$Time))))


Lm_presence_plot_nofac<-ggplot(Lm_enrichment_nofac, aes(x=Time, y=Percentage, fill=Lm))+
  geom_bar(stat='identity', color='black')+
  facet_grid(Year~Treatment)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  ylab("Percentage of samples") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_fill_manual(values = c('#F8CD9C','#EA7580'))+
  theme(legend.position = "bottom")+
  ggtitle("L. monocytogenes occurrance", subtitle = 'By enrichment method')+theme(plot.title = element_text(hjust = 0.5, face = 'italic'))
Lm_presence_plot_nofac
ggsave("Lm_presence_noFac.svg", plot=Lm_presence_plot_nofac, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("Lm_presence_noFac.png", plot=Lm_presence_plot_nofac, device="png", width=8, height=11, units="in", dpi=600)



## Treatment comparison: Before v after (all facilities)
# Testing each pair of results (i.e., F1T1Y1Before v F1T1Y1After) only results in significance is Before samples have 5 positive and after samples have 0 positives
# The test in underpowered to see effect of treatments due to small sample size
pos.test<-c(5,0)
tot.test<-c(5,5)
prop.test(pos.test, tot.test)


#Ask JK: is it worth to compare overall treatment differences (regardless of facility) - need difference of 7 samples to get a significant difference
#two-proportion test
#Y1
pos.T1Y1<-c(9,8)
tot.T1Y1<-c(15,15)
prop.test(pos.T1Y1, tot.T1Y1) #X-squared = 6.6928e-33, df = 1, p-value = 1

pos.T2Y1<-c(11,5)
tot.T2Y1<-c(15,15)
prop.test(pos.T2Y1, tot.T2Y1) #X-squared = 3.3482, df = 1, p-value = 0.06728

pos.T3Y1<-c(7,5)
tot.T3Y1<-c(15,15)
prop.test(pos.T3Y1, tot.T3Y1) #X-squared = 0.13889, df = 1, p-value = 0.7094

pos.T4Y1<-c(8,2)
tot.T4Y1<-c(15,15)
prop.test(pos.T4Y1, tot.T4Y1) #X-squared = 3.75, df = 1, p-value = 0.05281

#Y2
pos.T1Y2<-c(4,7)
tot.T1Y2<-c(15,15)
prop.test(pos.T1Y2, tot.T1Y2) #X-squared = 0.57416, df = 1, p-value = 0.4486

pos.T2Y2<-c(6,5)
tot.T2Y2<-c(15,15)
prop.test(pos.T2Y2, tot.T2Y2) #X-squared = 0, df = 1, p-value = 1

pos.T3Y2<-c(8,4)
tot.T3Y2<-c(15,15)
prop.test(pos.T3Y2, tot.T3Y2) #X-squared = 1.25, df = 1, p-value = 0.2636

pos.T4Y2<-c(8,6)
tot.T4Y2<-c(15,15)
prop.test(pos.T4Y2, tot.T4Y2) #X-squared = 0.13393, df = 1, p-value = 0.7144

#Note: none of the treatments applied significantly reduced the presence of Lm

#Both seasons combines
pos.T1<-c(13,15)
tot.T1<-c(30,30)
prop.test(pos.T1, tot.T1) #X-squared = 0.066964, df = 1, p-value = 0.7958

pos.T2<-c(17,10)
tot.T2<-c(30,30)
prop.test(pos.T2, tot.T2) #X-squared = 2.4242, df = 1, p-value = 0.1195

pos.T3<-c(15,9)
tot.T3<-c(30,30)
prop.test(pos.T3, tot.T3) #X-squared = 1.7361, df = 1, p-value = 0.1876

pos.T4<-c(16,8)
tot.T4<-c(30,30)
prop.test(pos.T4, tot.T4) #X-squared = 3.4028, df = 1, p-value = 0.06509



#### ATP measurements ####
#Upload data - brushes
ATP_brushes<-read_excel("Results_TwoYears.xlsx", sheet=5, col_names = TRUE)
ATP_brushes$ATP<-as.numeric(ATP_brushes$ATP)

#Subset by year
ATP_brushes_Y1<-subset(ATP_brushes, Year=="Y1")
ATP_brushes_Y2<-subset(ATP_brushes, Year=="Y2")

#Plot raw data
ATP_brushes_plot_Y1<-ggplot(ATP_brushes_Y1, aes(x=reorder(Time,Order), y=ATP, fill=Time))+
  geom_boxplot()+
  facet_grid(Facility~Site)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("ATP (RLU)") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(limits = c(0,250000))+
  ggtitle("ATP - Brushes Summary Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
ATP_brushes_plot_Y1
ggsave("Boxplot ATP brushes Y1.svg", plot=ATP_brushes_plot_Y1, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("Boxplot ATP brushes Y1.png", plot=ATP_brushes_plot_Y1, device="png", width=8, height=11, units="in", dpi=600)


ATP_brushes_plot_Y2<-ggplot(ATP_brushes_Y2, aes(x=reorder(Time,Order), y=ATP, fill=Time))+
  geom_boxplot()+
  facet_grid(Facility~Site)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("ATP (RLU)") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(limits = c(0,250000))+
  ggtitle("ATP - Brushes Summary Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
ATP_brushes_plot_Y2
ggsave("Boxplot ATP brushes Y2.svg", plot=ATP_brushes_plot_Y2, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("Boxplot ATP brushes Y2.png", plot=ATP_brushes_plot_Y2, device="png", width=8, height=11, units="in", dpi=600)


#Calculate stats by year
ATP_brushes_stat_Y1<-describeBy(ATP_brushes_Y1$ATP, list(ATP_brushes_Y1$Facility, ATP_brushes_Y1$Time, ATP_brushes_Y1$Site), mat = TRUE) #By facility and treatment
ATP_brushes_stat_Y2<-describeBy(ATP_brushes_Y2$ATP, list(ATP_brushes_Y2$Facility, ATP_brushes_Y2$Time, ATP_brushes_Y2$Site), mat = TRUE) #By facility and treatment

ATP_brushes_stat_plot_Y1<-ggplot(ATP_brushes_stat_Y1, aes(x=group2, y=mean, fill=group2))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  facet_grid(group1~group3)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("ATP (RLU)") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  #scale_y_continuous(breaks=c(0,2,4,6,8), limits = c(0,8))+
  ggtitle("ATP Brushes - Summary Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
ATP_brushes_stat_plot_Y1
ggsave("SummaryATP brushes Y1.svg", plot=ATP_brushes_stat_plot_Y1, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryATP brushes Y1.png", plot=ATP_brushes_stat_plot_Y1, device="png", width=8, height=11, units="in", dpi=600)


ATP_brushes_stat_plot_Y2<-ggplot(ATP_brushes_stat_Y2, aes(x=group2, y=mean, fill=group2))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  facet_grid(group1~group3)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("ATP (RLU)") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  #scale_y_continuous(breaks=c(0,2,4,6,8), limits = c(0,8))+
  ggtitle("ATP Brushes - Summary Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
ATP_brushes_stat_plot_Y2
ggsave("SummaryATP brushes Y2.svg", plot=ATP_brushes_stat_plot_Y2, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("SummaryATP brushes Y2.png", plot=ATP_brushes_stat_plot_Y2, device="png", width=8, height=11, units="in", dpi=600)


#ANOVA
anova_ATPbrushes_Y1<-aov(ATP ~ F_S_t, data=ATP_brushes_Y1)
summary(anova_ATPbrushes_Y1) #p=0.0933

anova_ATPbrushes_Y2<-aov(ATP ~ F_S_t, data=ATP_brushes_Y2)
summary(anova_ATPbrushes_Y2) #p=3.59e-07 ***


#tukey test

tukey_ATPbrushes_Y2<-HSD.test(anova_ATPbrushes_Y2, trt="F_S_t") 
#No significant differences between pairs of interest




#Upload ATP data - floors
ATP_floor<-read_excel("Results_TwoYears.xlsx", sheet=4, col_names = TRUE)
ATP_floor$ATP<-as.numeric(ATP_floor$ATP)

#Subset by year
ATP_floor_Y1<-subset(ATP_floor, Year=="Y1")
ATP_floor_Y2<-subset(ATP_floor, Year=="Y2")

#Plot raw data
ATP_floor_plot_Y1<-ggplot(ATP_floor_Y1, aes(x=reorder(Time,Order), y=ATP, fill=Time))+
  geom_boxplot()+
  facet_grid(Facility~Site)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("ATP (RLU)") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(limits = c(0,200000))+
  ggtitle("ATP - floor Summary Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
ATP_floor_plot_Y1
ggsave("Boxplot ATP floor Y1.svg", plot=ATP_floor_plot_Y1, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("Boxplot ATP floor Y1.png", plot=ATP_floor_plot_Y1, device="png", width=8, height=11, units="in", dpi=600)


ATP_floor_plot_Y2<-ggplot(ATP_floor_Y2, aes(x=reorder(Time,Order), y=ATP, fill=Time))+
  geom_boxplot()+
  facet_grid(Facility~Site)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("ATP (RLU)") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(limits = c(0,200000))+
  ggtitle("ATP - floor Summary Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
ATP_floor_plot_Y2
ggsave("Boxplot ATP floor Y2.svg", plot=ATP_floor_plot_Y2, device="svg", width=8, height=11, units="in", dpi=600)
ggsave("Boxplot ATP floor Y2.png", plot=ATP_floor_plot_Y2, device="png", width=8, height=11, units="in", dpi=600)


#ANOVA
anova_ATPfloor_Y1<-aov(ATP ~ F_S_time, data=ATP_floor_Y1)
summary(anova_ATPfloor_Y1) #p=0.0048 *

anova_ATPfloor_Y2<-aov(ATP ~ F_S_time, data=ATP_floor_Y2)
summary(anova_ATPfloor_Y2) #p= 0.0401 *


#tukey test
tukey_ATPfloor_Y1<-HSD.test(anova_ATPfloor_Y1, trt="F_S_time")
tukey_ATPfloor_Y2<-HSD.test(anova_ATPfloor_Y2, trt="F_S_time") 
#No significant differences between groups regardless of significant ANOVA






# Summary plots 
#Calculate mean and sd by treatment

ATP_stat<-describeBy(ATP_Y1$ATP, list(ATP_Y1$Fac_Trt_Time), mat = TRUE) #By facility and treatment
ATP_stat$Facility<-c(rep("P1", 8), rep("P2",8), rep("P3",8))
ATP_stat$Treatment<-rep(c(rep("T1", 2), rep("T2",2), rep("T3",2), rep("T4", 2)),3)
ATP_stat$Time<-rep(c("After","Before"),12)
ATP_stat$Time <- factor(ATP_stat$Time, levels = rev(levels(factor(ATP_stat$Time))))


ATP_stat_plot<-ggplot(ATP_stat, aes(x=Time, y=mean, fill=Time))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  facet_grid(Facility~Treatment)+ theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("ATP (RLU)") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  #scale_y_continuous(breaks=c(0,2,4,6,8), limits = c(0,8))+
  ggtitle("ATP - Summary")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
ATP_stat_plot
ggsave("SummaryATP.svg", plot=ATP_stat_plot, device="svg", width=8, height=11, units="in", dpi=600)

#Subset all before samples
#Cannot include the after samples because ATP was taken before application of the sanitizer

ATP_before<-subset(ATP_Y1, Time =="Before")

#Plot ATP v APC

ATP_before_plot<-ggplot(ATP_Y1, aes(x=ATP, y=logAPC, shape=Site, color=Facility))+
  geom_point(stat='identity', size=3)+facet_grid(Facility~.)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log CFU/swab") + xlab("ATP (RLU)")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_continuous(breaks=c(0,2,4,6,8,10,12), limits = c(0,12))+
  ggtitle("ATP - Before")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
ATP_before_plot
ggsave("ATPvAPC Before.svg", plot=ATP_before_plot, device="svg", width=8, height=11, units="in", dpi=600)


#### STOPPED HERE 4/20/23 ####






#### Nanopore Lm results ####
Nanopore<-read_excel("Nanopore.xlsx", sheet=2, col_names = TRUE)

Nanopore_stat_List<-describeBy(Nanopore$RAListeria, list(Nanopore$Lm_presence), mat = TRUE)
Nanopore_stat_List$Bact<-rep("Listeria spp.",2)
Nanopore_stat_Lm<-describeBy(Nanopore$RALm, list(Nanopore$Lm_presence), mat = TRUE)
Nanopore_stat_Lm$Bact<-rep("L. monocytogenes",2)

Nanopore_stat<-bind_rows(Nanopore_stat_List, Nanopore_stat_Lm)
Nanopore_stat$Order<-c(2,2,1,1)

Nano_plot<-ggplot(Nanopore_stat, aes(x=mean, y=reorder(Bact,Order), fill=group1))+
  geom_bar(stat='identity', color="black", position=position_dodge(0.9))+
  geom_errorbar(aes(xmin=mean-se, xmax=mean+se), width=.2, position=position_dodge(.9))+
  theme(legend.text=element_text(size=7,color='black'), legend.title= element_text(size=7, face='italic')) +
  ylab("") + xlab("Relative Abundance (%)")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text.x = element_text(color='black', size=7), axis.ticks = element_line(color='black'),
        axis.text.y = element_text(color='black', size=7,face='italic')) +
  theme(axis.title = element_text(size=7,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(strip.background= element_blank(), strip.text = element_text(size=7),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Nanopore Listeria - Summary")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_fill_manual(values = c('#F8CD9C','#EA7580'),name="L. monocytogenes")
Nano_plot
ggsave("SummaryNanoporeList.svg", plot=Nano_plot, device="svg", width=10, height=6, units="in", dpi=600)
ggsave("SummaryNanoporeList.png", plot=Nano_plot, device="png", width=10, height=6, units="in", dpi=600)






















#### PDA project ####

#### PDA: Aerobic plate counts ####
#Upload data
PDA_Y1<-read_excel("Results_Y1.xlsx", sheet=3, col_names = TRUE)


#Plot APC results y treatment
PDA_APCplot<-ggplot(PDA_Y1, aes(x=Site, y=logAPC, fill=Week))+
  geom_bar(stat='identity', position = 'dodge', color='black')+
  facet_grid(.~Treatment)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log CFU/swab") +
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(legend.position = 'bottom')+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Aerobic plate count - Long term P2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
PDA_APCplot
ggsave("PDAplot_APC.png", plot=PDA_APCplot, device="png", width=8, height=11, units="in", dpi=600)

#Summary plot
#Calculate mean and sd by week and treatment

PDA_APC_stat<-describeBy(PDA_Y1$logAPC, list(PDA_Y1$Treatment, PDA_Y1$Week), mat = TRUE) #By facility and treatment

PDA_APC_stat_plot<-ggplot(PDA_APC_stat, aes(x=group2, y=mean, fill=group2))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  facet_grid(.~group1)+
  theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log CFU/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Aerobic plate count - Summary")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
PDA_APC_stat_plot
ggsave("PDA_SummaryAPC.png", plot=PDA_APC_stat_plot, device="png", width=11, height=8, units="in", dpi=600)

#### PDA: MPN results ####

PDA_Y1$`Listspp_logMPN/swab`<-as.numeric(PDA_Y1$`Listspp_logMPN/swab`)
PDA_Y1$`Lm_logMPN/swab`<-as.numeric(PDA_Y1$`Lm_logMPN/swab`)

MPN_List_stat_PDA<-describeBy(PDA_Y1$`Listspp_logMPN/swab`, list(PDA_Y1$Treatment, PDA_Y1$Week), mat = TRUE)
MPN_Lm_stat_PDA<-describeBy(PDA_Y1$`Lm_logMPN/swab`, list(PDA_Y1$Treatment, PDA_Y1$Week), mat = TRUE) #By facility and treatment

List_stat_PDA<-ggplot(MPN_List_stat_PDA, aes(x=group2, y=mean, fill=group2))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  facet_grid(.~group1)+
  theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log MPN/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Listeria spp. - Long Term P2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  geom_hline(yintercept = 1.52, color='grey30', linetype="dashed", size=1)+
  guides(fill=guide_legend(title="Time"))
List_stat_PDA
ggsave("SummaryList-PDA.png", plot=List_stat_PDA, device="png", width=11, height=8, units="in", dpi=600)



Lm_stat_PDA<-ggplot(MPN_Lm_stat_PDA, aes(x=group2, y=mean, fill=group2))+
  geom_bar(stat='identity', position = position_dodge2(reverse=TRUE), color='black')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  facet_grid(.~group1)+
  theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=15, face='bold')) +
  ylab("log MPN/swab") + xlab("")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=15), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("L.monocytogenes - Long Term P2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  geom_hline(yintercept = 1.52, color='grey30', linetype="dashed", size=1)+
  guides(fill=guide_legend(title="Time"))
Lm_stat_PDA
ggsave("SummaryLm-PDA.png", plot=Lm_stat_PDA, device="png", width=11, height=8, units="in", dpi=600)
