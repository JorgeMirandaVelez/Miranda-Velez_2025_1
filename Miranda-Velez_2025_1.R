library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)

#OBS: Set working directory to source file location

#GAM analysis of suction cup NO3- measurements----
#Load and wrangle
nitN_SC <- read.csv(file="Miranda-Velez_2025_1_NO3_suction_cups.csv")
nitN_SC$Date <- date(nitN_SC$Date)
nitN_SC$Season <- factor(nitN_SC$Season)
nitN_SC$Block <- factor(nitN_SC$Block)
nitN_SC$Rotation <- factor(nitN_SC$Rotation)
nitN_SC$Tillage <- factor(nitN_SC$Tillage)
nitN_SC$Cover <- factor(nitN_SC$Cover)
str(nitN_SC)

#A few measurements fell below the limit of detection. Assigning value 0
nitN_SC$NitrateN[nitN_SC$NitrateN<0] <- 0

##GAM----
library(mgcv)
library(gratia)

#OBS: run-time ~10 minutes on Intel Core i5-1345U 1600 Mhz, 16 GB RAM.
ptm <- proc.time()
nit_mod <- bam(NitrateN ~ 
                 s(Sweek,bs="cr",k=15)+
                 s(Sweek,Rotation,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Tillage,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Cover,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Season,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Rotation,Tillage,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Rotation,Cover,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Rotation,Season,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Tillage,Cover,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Tillage,Season,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Cover,Season,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Rotation,Tillage,Cover,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Rotation,Tillage,Season,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Rotation,Cover,Season,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Tillage,Cover,Season,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Sweek,Rotation,Tillage,Cover,Season,bs="sz",xt=list(bs="cr"),k=15)+
                 s(Rotation,Tillage,Cover,Season,Block,bs="re"), #Effects by block structured as "random intercepts"
               method="fREML",
               family=tw(link="log"),
               discrete=T,
               data=nitN_SC)
proc.time() - ptm

appraise(nit_mod,method="simulate") #Diagnostic plots

s <- summary(nit_mod) #Model summary

nitN_SC <- nitN_SC %>% mutate(fit=fitted(nit_mod)) #Calculate fitted values

#Extract partial effects
par_eff_GAM <- smooth_estimates(nit_mod) %>% add_confint()
par_eff_GAM$est_res <- nit_mod$family$linkinv(par_eff_GAM$.estimate)
par_eff_GAM$upr <- nit_mod$family$linkinv(par_eff_GAM$.upper_ci)
par_eff_GAM$lwr <- nit_mod$family$linkinv(par_eff_GAM$.lower_ci)

#Extract conditional effects
library(marginaleffects)
cond_eff_GAM <- avg_predictions(model = nit_mod,
                                variables = c("Rotation", "Season"),
                                by = c('Sweek','Cover','Tillage'), 
                                type = "response") %>% data.frame()


#Aggregated annual N-leaching----
#Load and wrangle the interpolated NO3- and leaching calculated elsewhere
#For interpolation script please contact the corresponding author: jorge_mv@agro.au.dk
nit_inter <- read.csv(file="Miranda-Velez_2025_1_NO3_interpolated.csv")
nit_inter$Date <- date(nit_inter$Date)
nit_inter$Block <- factor(nit_inter$Block)
nit_inter$Rotation <- factor(nit_inter$Rotation)
nit_inter$Tillage <- factor(nit_inter$Tillage)
nit_inter$Cover <- factor(nit_inter$Cover)
nit_inter <- filter(nit_inter,is.na(Ci)==F&Date>"2019-07-16")
str(nit_inter)

#Aggregate to annual leaching by sampling seasons
nit_inter$Season <- year(nit_inter$Date-weeks(28)) %>% factor()
nit_agg <- nit_inter %>% group_by(Block,Rotation,Tillage,Cover,Season) %>% summarise(totBD=sum(Bottom_Drainage),totNleach=-1*sum(Nleach)) %>% data.frame()

##GLMMs----
library(glmmTMB)
library(DHARMa)

leach_mod_1 <- glmmTMB(totNleach~Rotation*Tillage*Cover*Season+
                         (1|Rotation:Tillage:Season/Block),
                     family=tweedie,
                     data=nit_agg)

plot(simulateResiduals(leach_mod_1)) #Diagnostic plots

car::Anova(leach_mod_1, type=2)
summary(leach_mod_1)

nit_agg$fit <- fitted(leach_mod_1) #Calculate fitted values

#Coniditonal estimates
library(emmeans)
library(multcomp)

#Cover and Tillage averaged across rotations and seasons
leach_mod_1 %>% emmeans("pairwise"~Cover|Tillage,type="response") %>% cld(Letters=letters,adjust="Tukey")
leach_mod_1 %>% emmeans("pairwise"~Tillage|Cover,type="response") %>% cld(Letters=letters,adjust="Tukey")

#Cover and season averaged across rotations and tillage
leach_mod_1 %>% emmeans(~Season|Cover,type="response") %>% cld(Letters=letters,adjust="Tukey")
leach_mod_1 %>% emmeans(~Cover|Season,type="response") %>% cld(Letters=letters,adjust="Tukey")

#Load weather data for Figure 1----
#This is the weather file used as input for simulated percolation using Daisy
meteo <- read.table(file="Foulum_2002_2024_Allerup_Daily.dwf",skip=21,sep=" ",header=T)[-1,] %>% type.convert(as.is=T) #Allerup corrected for Daisy
meteo$Date <- paste(meteo$Year,meteo$Month,meteo$Day) %>% ymd() %>% date()
meteo <- filter(meteo,Date>"2019-01-01"&Date<"2024-01-01")
meteo$ymonth <- tsibble::yearmonth(meteo$Date)
m_meteo <- meteo %>% group_by(ymonth) %>% summarize(st_date = min(Date), m_temp = mean(AirTemp), t_prec = sum(Precip)) %>% data.frame()
meteo$yweek <- tsibble::yearweek(meteo$Date)
w_meteo <- meteo %>% group_by(yweek) %>% summarize(st_date = min(Date), m_temp = mean(AirTemp), t_prec = sum(Precip)) %>% data.frame()

#Include the bottom drainage predicted by Daisy
d_drain <- nit_inter %>%
  group_by(Date) %>%
  summarise(mean_drain=mean(Bottom_Drainage))

m_drain <- d_drain %>% 
  mutate(yrm=tsibble::yearmonth(Date)) %>%
  group_by(yrm) %>%
  summarise(mdrain=sum(mean_drain)) %>% data.frame()
m_drain$st_date <- date(m_drain$yrm)


#Visualization----
library(patchwork)

#Figure 1
ggplot(m_meteo,aes(x=st_date,y=t_prec))+
  geom_col(fill="steelblue",alpha=0.8)+
  geom_col(data=m_drain,aes(x=st_date,y=mdrain),fill="indianred")+
  geom_line(data=w_meteo,aes(x=st_date,y=m_temp*5),col="black",linewidth=0.8)+
  scale_y_continuous(name="Monthly precipitaion/percolation (mm)",
                     sec.axis = sec_axis( trans=~.*0.2, name="Mean weekly temperature (°C)"))+
  scale_x_continuous(name = "",
                     limits = c(date("2019-01-01"),date("2024-01-01")),
                     breaks = c(date("2019-01-01"),date("2020-01-01"),date("2021-01-01"),date("2022-01-01"),date("2023-01-01"),date("2024-01-01")),
                     labels = c("2019","2020","2021","2022","2023","2024"))+
  theme_bw()+
  theme(axis.text = element_text(size=14),axis.title = element_text(size=14))


#Figure 2
a <- ggplot(filter(nitN_SC,Rotation=="R4"),aes(x=Date,y=NitrateN))+
  geom_point(col="gray30",alpha=0.6)+
  geom_line(aes(y=fit,group=interaction(Season,Block)),col="red3",linewidth=0.8)+
  facet_grid(Tillage~Cover)+
  theme_bw()+
  ggtitle("R4")+
  xlab("")+
  ylab(expression("nit-N (mg N L"^-1*")"))+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11, angle = 0))
b <- ggplot(filter(nitN_SC,Rotation=="R5"),aes(x=Date,y=NitrateN))+
  geom_point(col="gray30",alpha=0.6)+
  geom_line(aes(y=fit,group=interaction(Season,Block)),col="red3",linewidth=0.8)+
  facet_grid(Tillage~Cover)+
  theme_bw()+
  ggtitle("R5")+
  xlab("")+
  ylab(expression("nit-N (mg N L"^-1*")"))+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11, angle = 0))

a+b+plot_layout(nrow=2) #patchwork


#Figure 3
z <- ggplot(data=NULL,aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(data=filter(par_eff_GAM,.smooth=="s(Sweek)"),aes(),fill="gray70",col=NA,alpha=0.2)+
  geom_line(data=filter(par_eff_GAM,.smooth=="s(Sweek)"),aes(),col="black",linewidth=1)+
  theme_bw()+
  ylab(expression("mg N L"^-1))+
  theme(axis.text = element_text(size=14),axis.title = element_text(size=14))

a <- ggplot(data=NULL,aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(data=filter(par_eff_GAM,.smooth=="s(Sweek,Cover)"),aes(fill=Cover),col=NA,alpha=0.2)+
  geom_line(data=filter(par_eff_GAM,.smooth=="s(Sweek,Cover)"),aes(col=Cover),linewidth=1)+
  theme_bw()+
  ylab(expression("mg N L"^-1))+
  scale_color_manual(values = c("springgreen3","olivedrab4","turquoise4"))+
  scale_fill_manual(values = c("springgreen3","olivedrab4","turquoise4"))+
  theme(axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size = 10))

b <- ggplot(data=NULL,aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(data=filter(par_eff_GAM,.smooth=="s(Sweek,Tillage)"),aes(fill=Tillage),col=NA,alpha=0.2)+
  geom_line(data=filter(par_eff_GAM,.smooth=="s(Sweek,Tillage)"),aes(col=Tillage),linewidth=1)+
  theme_bw()+
  ylab(expression("mg N L"^-1))+
  scale_color_manual(values = c("firebrick2","tomato4","orange3"))+
  scale_fill_manual(values = c("firebrick2","tomato4","orange3"))+
  theme(axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size = 10))

c <- ggplot(data=NULL,aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(data=filter(par_eff_GAM,.smooth=="s(Sweek,Season)"),aes(fill=Season),col=NA,alpha=0.2)+
  geom_line(data=filter(par_eff_GAM,.smooth=="s(Sweek,Season)"),aes(col=Season),linewidth=1)+
  theme_bw()+
  ylab(expression("mg N L"^-1))+
  scale_color_manual(values = c("blue4","steelblue","cadetblue2","slategray3"))+
  scale_fill_manual(values = c("blue4","steelblue","cadetblue2","slategray3"))+
  theme(axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size = 10))+
  guides(color = guide_legend(nrow = 2))

d <- ggplot(data=NULL,aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(data=filter(par_eff_GAM,.smooth=="s(Sweek,Rotation)"),aes(fill=Rotation),col=NA,alpha=0.2)+
  geom_line(data=filter(par_eff_GAM,.smooth=="s(Sweek,Rotation)"),aes(col=Rotation),linewidth=1)+
  theme_bw()+
  ylab(expression("mg N L"^-1))+
  scale_color_manual(values = c("black","gray60"))+
  scale_fill_manual(values = c("black","gray60"))+
  theme(axis.text = element_text(size=14),axis.title = element_text(size=14),legend.text = element_text(size = 10))

#patchwork
z+a+b+c+d+guide_area()+plot_annotation(tag_levels = "a")+
  plot_layout(guides = "collect")&
  theme(legend.box = "vertical",
        legend.direction = "horizontal",
        legend.text = element_text(size=11),
        plot.tag.position = c(0, 1),
        plot.tag = element_text(size=16, hjust = -6, vjust = 1)
  )


#Figure 4
A <- ggplot(filter(par_eff_GAM,.smooth=="s(Sweek,Tillage,Cover)"),
            aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(aes(group=interaction(Cover,Tillage),fill=Tillage),col=NA,alpha=0.2)+
  geom_line(aes(col=Tillage),linewidth=1)+
  theme_bw()+
  ylab(expression("mg N L"^-1))+
  facet_wrap(~Cover,ncol=3)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11, angle = 0))+
  scale_color_manual(values = c("firebrick2","tomato4","orange3"))+
  scale_fill_manual(values = c("firebrick2","tomato4","orange3"))

B <- ggplot(filter(par_eff_GAM,.smooth=="s(Sweek,Tillage,Season)"),
            aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(aes(group=interaction(Season,Tillage),fill=Tillage),col=NA,alpha=0.2)+
  geom_line(aes(col=Tillage),linewidth=1)+
  theme_bw()+
  ylab(expression("mg N L"^-1))+
  facet_wrap(~Season,ncol=4)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11, angle = 0))+
  scale_color_manual(values = c("firebrick2","tomato4","orange3"))+
  scale_fill_manual(values = c("firebrick2","tomato4","orange3"))

C <- ggplot(filter(par_eff_GAM,.smooth=="s(Sweek,Cover,Season)"),
            aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(aes(group=interaction(Cover,Season),fill=Cover),col=NA,alpha=0.2)+
  geom_line(aes(col=Cover),linewidth=1)+
  theme_bw()+
  ylab(expression("mg N L"^-1))+
  facet_wrap(~Season,ncol=4)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11, angle = 0))+
  scale_color_manual(values = c("springgreen3","olivedrab4","turquoise4"))+
  scale_fill_manual(values = c("springgreen3","olivedrab4","turquoise4"))

D <- ggplot(filter(par_eff_GAM,.smooth=="s(Sweek,Rotation,Cover)"),
            aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(aes(group=interaction(Cover,Rotation),fill=Cover),col=NA,alpha=0.2)+
  geom_line(aes(col=Cover),linewidth=1)+
  theme_bw()+
  ylab(expression("mg N L"^-1))+
  facet_wrap(~Rotation,ncol=4)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11, angle = 0))+
  scale_color_manual(values = c("springgreen3","olivedrab4","turquoise4"))+
  scale_fill_manual(values = c("springgreen3","olivedrab4","turquoise4"))

E <- ggplot(filter(par_eff_GAM,.smooth=="s(Sweek,Rotation,Tillage)"),
            aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(aes(group=interaction(Tillage,Rotation),fill=Tillage),col=NA,alpha=0.2)+
  geom_line(aes(col=Tillage),linewidth=1)+
  theme_bw()+
  ylab(expression("mg N L"^-1))+
  facet_wrap(~Rotation,ncol=4)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11, angle = 0))+
  scale_color_manual(values = c("firebrick2","tomato4","orange3"))+
  scale_fill_manual(values = c("firebrick2","tomato4","orange3"))+
  scale_y_continuous(position = "left")


G <- ggplot(filter(par_eff_GAM,.smooth=="s(Sweek,Rotation,Season)"),
            aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(aes(group=interaction(Season,Rotation),fill=Rotation),col=NA,alpha=0.2)+
  geom_line(aes(col=Rotation),linewidth=1)+
  theme_bw()+
  ylab(expression("mg N L"^-1))+
  facet_wrap(~Season,ncol=4)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11, angle = 0))+
  scale_color_manual(values = c("black","gray60"))+
  scale_fill_manual(values = c("black","gray60"))

#patchwork
des <- "
AAA#
BBBB
CCCC
DDDD
EEGG
FFGG
"

A+B+C+G+E+D+guide_area()+plot_annotation(tag_levels = "a")+
  plot_layout(nrow=6,guides = "collect",design = des)&
  theme(legend.box = "vertical",
        legend.direction = "horizontal",
        plot.tag.position = c(0, 1),
        plot.tag = element_text(size=14, hjust = 0, vjust = 0))

#Appendix 1: Third-order interactions
ggplot(filter(par_eff_GAM,.smooth=="s(Sweek,Tillage,Cover,Season)"),
       aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(aes(group=interaction(Cover,Tillage,Season),fill=Season),col=NA,alpha=0.2)+
  geom_line(aes(col=Season),linewidth=1)+
  theme_bw()+
  ylab("s(Sweek,Cover,Tillage,Season)")+
  facet_grid(Tillage~Cover)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11, angle = 0))+
  scale_color_manual(values = c("blue4","steelblue","cadetblue2","slategray3"))+
  scale_fill_manual(values = c("blue4","steelblue","cadetblue2","slategray3"))

ggplot(filter(par_eff_GAM,.smooth=="s(Sweek,Rotation,Tillage,Cover)"),
       aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(aes(group=interaction(Cover,Tillage,Rotation),fill=Rotation),col=NA,alpha=0.2)+
  geom_line(aes(col=Rotation),linewidth=1)+
  theme_bw()+
  ylab("s(Sweek,Rotation,Cover,Tillage)")+
  facet_grid(Tillage~Cover)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11, angle = 0))+
  scale_color_manual(values = c("black","gray60"))+
  scale_fill_manual(values = c("black","gray60"))

#Appendix 2: Fourth-order interactions
ggplot(filter(par_eff_GAM,.smooth=="s(Sweek,Rotation,Tillage,Cover,Season)"),
       aes(x=Sweek,y=est_res,ymax=upr,ymin=lwr))+
  geom_ribbon(aes(group=interaction(Rotation,Cover,Tillage,Season),fill=Season),col=NA,alpha=0.2)+
  geom_line(aes(col=Season),linewidth=1)+
  theme_bw()+
  ylab("s(Sweek,Rotation,Cover,Tillage,Season)")+
  facet_grid(Rotation~Tillage~Cover)+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11),
        strip.text.y = element_text(size = 11, angle = 0))+
  scale_color_manual(values = c("blue4","steelblue","cadetblue2","slategray3"))+
  scale_fill_manual(values = c("blue4","steelblue","cadetblue2","slategray3"))


#Figure 5
x <- ggplot(cond_eff_GAM, aes(x=Sweek,y=estimate))+
  geom_ribbon(aes(ymin = conf.low,ymax=conf.high,fill=Tillage),alpha=0.3)+
  geom_line(aes(col=Tillage)) +
  facet_wrap(~Cover,ncol=3,nrow=4)+
  ylab(expression("Estimated nit-N (mg N l"^-1*")"))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11))+
  xlab("Sweek")+
  scale_color_manual(values = c("firebrick2","tomato4","orange3"))+
  scale_fill_manual(values = c("firebrick2","tomato4","orange3"))

y <- ggplot(cond_eff_GAM, aes(x=Sweek,y=estimate))+
  geom_ribbon(aes(ymin = conf.low,ymax=conf.high,fill=Cover),alpha=0.3)+
  geom_line(aes(col=Cover)) +
  facet_wrap(~Tillage,ncol=3,nrow=4)+
  ylab(expression("Estimated nit-N (mg N l"^-1*")"))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 11))+
  xlab("Sweek")+
  scale_color_manual(values = c("springgreen3","olivedrab4","turquoise4"))+
  scale_fill_manual(values = c("springgreen3","olivedrab4","turquoise4"))

#patchwork
x+y+plot_layout(nrow=2)
