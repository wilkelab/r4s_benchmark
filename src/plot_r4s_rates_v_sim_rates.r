library(tidyr)
library(ggplot2)
library(dplyr)
library(cowplot)

t1 <- list.files("r4s_site_rates/final_site_rates",full.names=T)
t2 <- list.files("sim_site_rates/final_site_rates",full.names=T)

for (i in 1:length(t1)) 
{
	taxa_num <- gregexpr()("t[:digit:]+",t1[i],value = T)	
	
	t <- read.table(t1[i],skip=11,sep="\t")

	#reformat rate4site output
	t <- t %>% separate(V1,into=c("pos","seq","score","e","qq_int_lower","qq_int_upper","stdev","msa_data"),sep="[:blank:]+|[:blank:]*\\[[:blank:]*|[:blank:]*\\,[:blank:]*|[:blank:]*\\][:blank:]*",extra="drop")
	r4s <- data.frame(pos = as.numeric(t$pos),
			seq = t$seq,
			score = as.numeric(t$score),
			qq_int_lower = as.numeric(t$qq_int_lower),
			qq_int_upper = as.numeric(t$qq_int_upper),
			stdev = as.numeric(t$stdev),
			msa_data = t$msa_data
			)
	
	s <- read.table(t2[i],sep="\t",header=T)

	s$r4s_score <- r4s$score
	cor <- cor(s$dN,s$r4s_score,method="spearman")

	p <- ggplot(s,aes(dN,r4s_score)) + 
	geom_point() + 
	geom_smooth(method=lm) +
	xlab("simulated rate (dN)") +
	ylab("rate4site score") +
	scale_x_continuous(breaks=seq(0.1,1.5,0.1), limits = c(0.1,1.5)) + 
	scale_y_continuous(breaks=seq(-1,6,1), limits = c(-1,6)) +
	annotate("text",x=0.2,y=6,label=paste("rho =",round(cor,2)))
	
}




