library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

model_lst = c("mech_codon_dN","mech_codon_dN_dS","mut_sel_dN","mut_sel_dN_dS")

##for (model in model_lst) {
t1 <- list.files(paste0(model,"/r4s_site_rates/norm_rates"),full.names=T)
info = file.info(t1)
t1 <- t1[info$size != 0]

if (model=="mech_codon_dN" | model=="mech_codon_dN_dS"){
	t2 <- list.files(paste0(model,"/sim_site_rates/merged_output"),
	full.names=T)
	}
	
if (model=="mut_sel_dN" | model=="mut_sel_dN_dS"){
t2 <- list.files(paste0(model,"/sim_site_rates/inferred_rates"),
	full.names=T)
	}
t2 <- t2[info$size != 0]

for (i in 1:length(t1)) 
{
	print(t1[i])
	print(t2[i])

	r1 <- read.table(t1[i],skip=11,sep="\t")
	
	#reformat rate4site output
	r <- r1 %>% separate(V1,into=c("c1","c2","c3","c4","c5","c6","c7","c8"),sep="\\][:blank:]*|[:blank:]*\\[[:blank:]*|\\,[:blank:]*|[:blank:]+",extra="drop") %>%
	  separate(c8,into=c("c8","c9"),sep="\\/") %>% 
	  select(-c(c1,c9))
	colnames(r) <- c("pos","seq","score","qq_int_lower","qq_int_upper","stdev","num_taxa")
	
	if (model=="mech_codon_dN" | model=="mech_codon_dN_dS"){
    s <- read.table(t2[i],sep="\t",header=T,col.names=c("pos","rate_probability","dN","dS"))
	  s$dN.dS <- as.numeric(s$dN)/as.numeric(s$dS) 
  }
	if (model=="mut_sel_dN" | model=="mut_sel_dN_dS"){
	  s <- read.table(t2[i],sep="\t",header=T,col.names=c("pos","dN","dS"))
	  s$dN.dS <- as.numeric(s$dN)/as.numeric(s$dS)
	}
	
	s$r4s_score <- as.numeric(r$score)
	s$num_taxa <- as.numeric(r$num_taxa)
	
	str <- regexpr("b\\d+",t1[i])[1]
	end <- regexpr("_\\d+.txt",t1[i])[1]
	bl <- as.numeric(substr(t1[i],str+1,end-1))
	
	str <- regexpr("_\\d+.txt",t1[i])[1]
	end <- regexpr(".txt",t1[i])[1]
	sim_num <- as.numeric(substr(t1[i],str+1,end-1))
	
	s$branch_len <- rep(bl,length(s$num_taxa))
	s$sim_num <- rep(sim_num,length(s$num_taxa))

	if (model=="mut_sel_dN" | model=="mut_sel_dN_dS" | model=="mech_codon_dN_dS"){
	  cor_test <- cor.test(s$dN.dS,s$r4s_score,method="spearman")
	  cor <- cor_test$estimate
	  p_val <- cor_test$p.value
	}
	
	if (model=="mech_codon_dN"){
	  cor_test <- cor.test(s$dN,s$r4s_score,method="spearman")
	  cor <- cor_test$estimate
	  p_val <- cor_test$p.value
	}
	
	s$cor <- rep(cor,length(s$num_taxa))
	s$p_val <- rep(p_val,length(s$num_taxa))
  
	if (i==1) {
		d <- s
	} else d <- rbind(d, s)
}

d <- d %>% arrange(num_taxa,branch_len)

a <- d %>% filter(sim_num==1)
sig <- rep(" ",length(a$num_taxa))
sig[a$p_val <= 0.05] = rep("*",length(sig[a$p_val <= 0.05]))

if (model == "mech_codon_dN") {
  p1 <- ggplot(a,aes(dN,r4s_score)) + 
    geom_point(size=1,alpha=0.7) + 
    geom_smooth(method=lm) +
    xlab(expression(bold("simulated rate (dN/dS)"))) +
    ylab("rate4site score") +
    theme(axis.text=element_text(size=8),legend.position="none") +
    geom_text(aes(x=0.3,y=5,label=paste0(round(cor,2),sig),size=4)) +
    scale_x_continuous(breaks=seq(0,1,0.5), labels=c("0","0.5","1"), limits = c(0,1)) + 
    scale_y_continuous(breaks=seq(-2,6,2), limits = c(-2,6)) +
    facet_grid(num_taxa ~ branch_len) +
    background_grid(major = 'xy', minor = "none") + 
    panel_border()
  ggsave(paste0("plots/",model,"_r4s_rates_v_sim_rates.png"))
}
if ( model == "mech_codon_dN_dS") {
	p1 <- ggplot(a,aes(dN.dS,r4s_score)) + 
		geom_point(size=1,alpha=0.7) + 
		geom_smooth(method=lm) +
		xlab(expression(bold("simulated rate (dN/dS)"))) +
		ylab("rate4site score") +
		theme(axis.text=element_text(size=8),legend.position="none") +
		geom_text(aes(x=0.3,y=5,label=paste0(round(cor,2),sig),size=4)) +
		scale_x_continuous(breaks=seq(0,1.5,0.5), labels=c("0","0.5","1","1.5"), limits = c(0,1.5)) + 
		scale_y_continuous(breaks=seq(-2,6,2), limits = c(-2,6)) +
		facet_grid(num_taxa ~ branch_len) +
		background_grid(major = 'xy', minor = "none") + 
		panel_border()
	ggsave(paste0("plots/",model,"_r4s_rates_v_sim_rates.png"))
}

b1 <- filter(a,branch_len==0.33,num_taxa==32)
index <- a$branch_len==0.1 & a$num_taxa==32 | a$branch_len==0.33 & a$num_taxa==64 | a$branch_len==0.1 & a$num_taxa==64
b2a <- a[index,]
b2 <- rbind(b1,b2a)
if (model == "mut_sel_dN") {
  p1 <- ggplot(a,aes(dN,r4s_score)) + 
  geom_point(size=1,alpha=0.7) + 
  geom_smooth(method=lm) +
  xlab(expression(bold("predicted rate (dN/dS)"))) +
  ylab("rate4site score") +
  theme(axis.text=element_text(size=8),legend.position="none") +
  geom_text(aes(x=0.3,y=5,label=paste0(round(cor,2),sig),size=4)) +
  scale_x_continuous(breaks=seq(0,1,0.5), labels=c("0","0.5","1.0"), limits = c(0.0,1.0)) + 
  scale_y_continuous(breaks=seq(-2,6,2), limits = c(-2,6)) +
  facet_grid(num_taxa ~ branch_len) +
  background_grid(major = 'xy', minor = "none") + 
  panel_border()
  ggsave(paste0("plots/",model,"_r4s_rates_v_sim_rates.png"))
  
  p1a <- ggplot(b1,aes(dN.dS,r4s_score)) + 
    geom_point(size=3,alpha=0.7) + 
    geom_smooth(method=lm) +
    xlab(expression(bold("predicted rate (dN/dS)"))) +
    ylab("rate4site score") +
    theme(axis.text=element_text(size=8),legend.position="none") +
    geom_text(aes(x=0.1,y=5,label=paste0(round(cor,2),sig)),size=6) +
    scale_x_continuous(breaks=seq(0,1,0.1), 
                       #labels=c("0","0.1","1.0"), 
                       limits = c(0.0,1.0)) + 
    scale_y_continuous(breaks=seq(-2,6,1), limits = c(-2,6)) +
    facet_grid(num_taxa ~ branch_len) +
    background_grid(major = 'xy', minor = "none") + 
    panel_border()+
   theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
  ggsave(paste0("plots/",model,"_one_var_r4s_rates_v_sim_rates.png"))
  
  p1b <- ggplot(b2,aes(dN.dS,r4s_score)) + 
    geom_point(size=3,alpha=0.7) + 
    geom_smooth(method=lm) +
    xlab(expression(bold("predicted rate (dN/dS)"))) +
    ylab("rate4site score") +
    theme(axis.text=element_text(size=8),legend.position="none") +
    geom_text(aes(x=0.1,y=5,label=paste0(round(cor,2),sig)),size=6) +
    scale_x_continuous(breaks=seq(0,1,0.2), 
                       #labels=c("0","0.1","1.0"), 
                       limits = c(0.0,1.0)) + 
    scale_y_continuous(breaks=seq(-2,6,2), limits = c(-2,6)) +
    facet_grid(num_taxa ~ branch_len) +
    background_grid(major = 'xy', minor = "none") + 
    panel_border() +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 17),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  ggsave(paste0("plots/",model,"_two_var_r4s_rates_v_sim_rates.png"))

  }

if (model == "mut_sel_dN_dS") {
	p1 <- ggplot(a,aes(dN.dS,r4s_score)) + 
		geom_point(size=1,alpha=0.7) + 
		geom_smooth(method=lm) +
	  xlab(expression(bold("predicted rate (dN/dS)"))) +
		ylab("rate4site score") +
		theme(axis.text=element_text(size=8),legend.position="none") +
		geom_text(aes(x=0.4,y=5,label=paste0(round(cor,2),sig),size=4)) +
		scale_x_continuous(breaks=seq(0,1.5,0.5), labels=c("0","0.5","1","1.5"), limits = c(0.0,1.5)) + 
		scale_y_continuous(breaks=seq(-2,6,2), limits = c(-2,6)) +
		facet_grid(num_taxa ~ branch_len) +
		background_grid(major = 'xy', minor = "none") + 
		panel_border()
	ggsave(paste0("plots/",model,"_r4s_rates_v_sim_rates.png"))
}

colfunc <- colorRampPalette(c("orange","red","darkred"))
p2 <- ggplot(d,aes(num_taxa,cor,color=factor(branch_len),group=branch_len)) + 
  stat_summary(fun.y = mean,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "pointrange",
               size=0.8) +
  scale_colour_manual(values=colfunc(8)) +
  guides(col = guide_legend(title="Branch Length",reverse = TRUE)) +
  stat_summary(fun.y = mean,geom = "line",size=0.8,aes(color=factor(branch_len)))+
  xlab("Number of Taxa") +
  ylab("Correlation (spearman)") +
  scale_x_continuous(breaks=seq(0,300,50), labels=c("0"," ","100"," ","200"," ","300"), limits = c(0,300)) + 
	scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1)) +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
ggsave(paste0("plots/",model,"_cor_v_num_taxa.png"))

colfunc <- colorRampPalette(c("cyan2","dodgerblue2","navyblue"))
if (model=="mech_codon_dN" | model=="mech_codon_dN_dS") {
  p3 <- ggplot(d,aes(branch_len,cor,colour=factor(num_taxa),group=num_taxa)) + 
    stat_summary(fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange",
                 size=0.8) +
    scale_colour_manual(values=colfunc(4)) +
    guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
    stat_summary(fun.y = mean,geom = "line",size=0.8,aes(color=factor(num_taxa)))+
    xlab("Branch Length") +
    ylab("Correlation (spearman)") +
    scale_x_log10(breaks=round(c(0.001,0.0033,0.01,0.033,0.1),3),
                  limits=c(0.001,0.1)) + 
    scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1))+ 
    theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))
  ggsave(paste0("plots/",model,"_cor_v_branch_len.png"))
}

if (model=="mut_sel_dN" | model=="mut_sel_dN_dS") {
  p3 <- ggplot(d,aes(branch_len,cor,colour=factor(num_taxa),group=num_taxa)) + 
    stat_summary(fun.y = mean,
                 fun.ymin = function(x) mean(x) - sd(x), 
                 fun.ymax = function(x) mean(x) + sd(x), 
                 geom = "pointrange",
                 size=0.8) +
    scale_colour_manual(values=colfunc(4)) +
    guides(col = guide_legend(title="Number of Taxa",reverse = TRUE)) +
    stat_summary(fun.y = mean,geom = "line",size=0.8,aes(color=factor(num_taxa)))+
    xlab("Branch Length") +
    ylab("Correlation (spearman)") +
    scale_x_log10(breaks=round(c(0.001,0.0033,0.01,0.033,0.1,0.33,1.0,3.3),3),
                  #labels=c(0.001,0.003,0.01,0.033,0.1,0.33,1.0,3.3),
                  limits=c(0.001,3.3)) + 
    scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1)) +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 17),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  ggsave(paste0("plots/",model,"_cor_v_branch_len.png"))
}
}

