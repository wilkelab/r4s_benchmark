library(tidyr)
library(dplyr)

#read in r4s output
r_bias <- read.csv("other_methods_results/r4s_site_rates/processed_rates/full_results_r4s_bias.csv",row.names=1)
r_nobias <- read.csv("other_methods_results/r4s_site_rates/processed_rates/full_results_r4s_nobias.csv",row.names=1)

#read in other methods output
#files acquired out of github repository sitewise_dnds_mutsel/results/processed_results/
d_bias <- read.csv("other_methods_results/other_methods_dnds/full_results_bias_dataset.csv")
d_nobias <- read.csv("other_methods_results/other_methods_dnds/full_results_nobias_dataset.csv")

##select true1 and true2 values
s_bias <- filter(d_bias,rep==1,ntaxa==128,bl==0.0025,method=="FEL1")
s_nobias <- filter(d_nobias,rep==1,ntaxa==128,bl==0.0025,method=="FEL1")
num_bias <- length(r_bias$num_taxa)/length(s_bias$ntaxa)
num_nobias <- length(r_nobias$num_taxa)/length(s_nobias$ntaxa)

##other methods file format "ntaxa","bl","site","rep","true1","true2","type","dnds","method"
##reformate r4s output to match other methods output csv
ref_bias_r <- data.frame(ntaxa=r_bias$num_taxa,
                         bl=r_bias$bl,
                         site=r_bias$pos,
                         rep=r_bias$rep,
                         true1=rep(s_bias$true1,num_bias),
                         true2=rep(s_bias$true2,num_bias),
                         type=r_bias$type,
                         dnds=r_bias$score,
                         method=rep("Rate4Site",num_bias)
)
ref_nobias_r <- data.frame(ntaxa=r_nobias$num_taxa,
                           bl=r_nobias$bl,
                           site=r_nobias$pos,
                           rep=r_nobias$rep,
                           true1=rep(s_bias$true1,num_nobias),
                           true2=rep(s_bias$true2,num_nobias),
                           type=r_nobias$type,
                           dnds=r_nobias$score,
                           method=rep("Rate4Site",num_nobias)
)

combined_bias <- rbind(d_bias,ref_bias_r)
combined_nobias <- rbind(d_nobias,ref_nobias_r)
t_bias <- combined_bias %>% 
  filter(method!="FEL2_1" & method!="FUBAR2_1") %>%
  arrange(method,bl,ntaxa,rep)
t_nobias <- combined_nobias %>% 
  filter(method!="FEL2_1" & method!="FUBAR2_1") %>%
  arrange(method,bl,ntaxa,rep)

cor_list1 <- c()
cor_list2 <- c()
start <- seq(1,length(t_bias$ntaxa),100)
for (i in start) {
  a <- slice(t_bias,i:(i+99))
  if (all(a$site==c(1:100))) {
    cor_test1 <- cor.test(a$dnds,a$true1,method="spearman",na.action = na.omit)
    cor_test2 <- cor.test(a$dnds,a$true2,method="spearman",na.action = na.omit)
    cor1 <- cor_test1$estimate
    cor2 <- cor_test2$estimate
    #p_val <- cor_test$p.value
    cor_list1 <- append(cor_list1,rep(cor1,length(a$ntaxa)))
    cor_list2 <- append(cor_list2,rep(cor2,length(a$ntaxa)))
  }
  else {
    print("Wrong sites used for correlations")
    break
  }
}
t_bias$cor1 <- cor_list1
t_bias$cor2 <- cor_list2

cor_list1 <- c()
cor_list2 <- c()
start <- seq(1,length(t_nobias$ntaxa),100)
for (i in start) {
  a <- slice(t_nobias,i:(i+99))
  if (all(a$site==c(1:100))) {
    cor_test1 <- cor.test(a$dnds,a$true1,method="spearman",na.action = na.omit)
    cor_test2 <- cor.test(a$dnds,a$true2,method="spearman",na.action = na.omit)
    cor1 <- cor_test1$estimate
    cor2 <- cor_test2$estimate
    #p_val <- cor_test$p.value
    cor_list1 <- append(cor_list1,rep(cor1,length(a$ntaxa)))
    cor_list2 <- append(cor_list2,rep(cor2,length(a$ntaxa)))
  }
  else {
    print("Wrong sites used for correlations")
    break
  }
}
t_nobias$cor1 <- cor_list1
t_nobias$cor2 <- cor_list2

write.csv(t_bias,file="other_methods_results/all_methods_processed_rates/all_methods_bias.csv",quote=F)
write.csv(t_nobias,file="other_methods_results/all_methods_processed_rates/all_methods_nobias.csv",quote=F)
