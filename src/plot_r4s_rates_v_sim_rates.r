t <- read.table("../test_r4s_output/r4s.res",skip=11,sep="\t")

#reformat rate4site output
t <- t %>% separate(V1,into=c("a","pos","seq","score","e","qq_int_lower","qq_int_upper","stdev","msa_data"),sep="[:blank:]+|[:blank:]*\\[[:blank:]*|[:blank:]*\\,[:blank:]*|[:blank:]*\\][:blank:]*",extra="drop")
r4s <- data.frame(pos = as.numeric(t$pos),
				seq = t$seq,
				score = as.numeric(t$score),
				qq_int_lower = as.numeric(t$qq_int_lower),
				qq_int_upper = as.numeric(t$qq_int_upper),
				stdev = as.numeric(t$stdev),
				msa_data = t$msa_data
				)
				
s <- read.table("../test_sim_site_rates/site_rate_final.txt",sep="\t",header=T)

s$r4s_score <- r4s$score
p <- ggplot(s,aes(dN,r4s_score))
p <- p + geom_point()

