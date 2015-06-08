library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

t1 <- list.files("r4s_site_rates/final_site_rates",full.names=T)
t2 <- list.files("sim_site_rates/final_site_rates",full.names=T)

for (i in 1:length(t1)) 
{
	str <- regexpr("b\\d+",t1[i])[1]
	end <- regexpr(".txt",t1[i])[1]
	bl <- as.numeric(substr(t1[i],str+1,end-1))
	
	r1 <- read.table(t1[i],skip=11,sep="\t")

	#reformat rate4site output
	r <- r1 %>% separate(V1,into=c("pos","seq","score","none1","qq_int_lower","qq_int_upper","stdev","msa_data"),sep="[:blank:]+|[:blank:]*\\[[:blank:]*|[:blank:]*\\,[:blank:]*|[:blank:]*\\][:blank:]*",extra="drop") %>% 
	separate(msa_data,into=c("num_taxa","none2"),sep="\\/") %>% 
	select(-c(none1,none2))

	s <- read.table(t2[i],sep="\t",header=T,col.names=c("pos","rate_probability","dN","dS"))

	s$r4s_score <- as.numeric(r$score)
	s$num_taxa <- as.numeric(r$num_taxa)
	s$branch_len <- rep(bl,length(s$num_taxa))
	
	cor_test <- cor.test(s$dN,s$r4s_score,method="spearman")
	cor <- cor_test$estimate
	p_val <- cor_test$p.value

	
	s$cor_p_val <- rep(p_val,length(s$num_taxa))
	s$cor <- rep(cor,length(s$num_taxa))
	
	if (i==1) {
		d <- s
	} else d <- rbind(d, s)
}

d <- d %>% arrange(num_taxa,branch_len)
sig <- rep(" ",length(d$num_taxa))
sig[d$cor_p_val <= 0.05] = rep("*",length(sig[d$cor_p_val <= 0.05]))

p1 <- ggplot(d,aes(dN,r4s_score)) + 
	geom_point() + 
	geom_smooth(method=lm) +
	xlab("simulated rate (dN)") +
	ylab("rate4site score") +
	theme(axis.text=element_text(size=8),legend.position="none") +
	geom_text(aes(x=0.3,y=5,label=paste0(round(cor,2),sig),size=8)) +
	scale_x_continuous(breaks=seq(0.1,1.5,0.2), limits = c(0.1,1.5)) + 
	scale_y_continuous(breaks=seq(-1,6,2), limits = c(-1,6)) +
	facet_grid(num_taxa ~ branch_len)

c2 <- cor.test(d$num_taxa,d$cor,method="spearman")
if (c2$p.value <= 0.05) { sig="*"
} else sig=""

p2 <- ggplot(d,aes(num_taxa,cor)) + 
	geom_point() +
	geom_smooth(method=lm) +
	xlab("Number of Taxa") +
	ylab("Correlation (spearman)") +
	scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1)) +
	geom_text(aes(x=30,y=0.95,label=paste0(round(c2$estimate,2),sig)))


c3 <- cor.test(d$branch_len,d$cor,method="spearman")
if (c3$p.value <= 0.05) { sig="*"
} else sig=""
p3 <- ggplot(d,aes(branch_len,cor)) + 
	geom_point() +
	geom_smooth(method=lm) +
	xlab("Branch Length") +
	ylab("Correlation (spearman)") +
	scale_y_continuous(breaks=seq(0,1,0.2), limits = c(0,1)) +
	geom_text(aes(x=0.001,y=0.95,label=paste0(round(c3$estimate,2),sig)))

ggsave("plots/r4s_rates_v_sim_rates.pdf",p1)
ggsave("plots/cor_v_num_taxa.pdf",p2)
ggsave("plots/cor_v_branch_len.pdf",p3)