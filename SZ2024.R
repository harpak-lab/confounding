# Install GenomicSEM package first. Follow instructions in https://github.com/GenomicSEM/GenomicSEM

require(GenomicSEM)

datadir <- "sumstats/"
set <- "sz_rep" # replicate SZ's results

# get # children (both sexes) sumstat file
sz_nc <- paste0(datadir, set, "/ukb_wb.quant.both_n_children.regenie")
sz_nc_new <- gsub(".regenie", "_p.regenie", sz_nc)
ss <- read_delim(sz_nc, delim=" ") %>% mutate(P=10^(-LOG10P))
write_tsv(ss, sz_nc_new)

# list all male sumstat files
regenie_results_all <- list.files(paste0(datadir, set), pattern="*[.]male*", full.names=T)

# get sumstat files with transformed P-vals
regenie_results <- regenie_results_all[!grepl("_p.regenie", regenie_results_all)]

# check if sumstat files with transformed P-vals exist, otherwise create them
for (i in 1:length(regenie_results)){
  outfile <- gsub(".regenie", "_p.regenie", regenie_results[i])
  if (!file.exists(outfile)){
    print(regenie_results[i])
    ss <- read_delim(regenie_results[i], delim=" ") %>% mutate(P=10^(-LOG10P))
    write_tsv(ss, outfile)
  }
}

regenie_results_new <- gsub(".regenie", "_p.regenie", regenie_results)


# specify sumstat files and trait names
trait_files_new <- c(sz_nc_new, regenie_results_new)

trait.names<-c("NC", gsub(".*.male[_,.]", "", gsub("_p.regenie", "", basename(regenie_results_new))))

##########################
# Munge GWAS summary statistics first
##########################
hm3 <- "w_hm3.snplist"

# doing by batch to avoid hitting memory limits
munge(files=trait_files_new[1:10],hm3=hm3,trait.names=trait.names[1:10])
munge(files=trait_files_new[11:20],hm3=hm3,trait.names=trait.names[11:20])
munge(files=trait_files_new[21:31],hm3=hm3,trait.names=trait.names[21:31])

# get prevalences for binary traits
binary_prevs <- phenos_new %>%
  dplyr::filter(genetic_sex==1) %>%
  dplyr::filter(FID %in% keep_ids$eid) %>%
  mutate(Recent_poor_appetite_or_overeating=ifelse(Recent_poor_appetite_or_overeating <0, NA,
                                                   ifelse(Recent_poor_appetite_or_overeating ==1, 0, 1)),
         Maternal_smoking_around_birth =ifelse(Maternal_smoking_around_birth  <0, NA, Maternal_smoking_around_birth  ),
         Ever_taken_cannabis =ifelse(Ever_taken_cannabis  <0, NA,
                                     ifelse(Ever_taken_cannabis ==0, 0, 1)),
         afs10   =ifelse(Age_first_had_sexual_intercourse  < 0, NA,
                         ifelse(Age_first_had_sexual_intercourse < 10, 1, 0)),
         afs13   =ifelse(Age_first_had_sexual_intercourse  < 0, NA,
                         ifelse(Age_first_had_sexual_intercourse < 13, 1, 0)),
         Physically_abused_by_family_as_a_child  =ifelse(Physically_abused_by_family_as_a_child   <0, NA,
                                                         ifelse(Physically_abused_by_family_as_a_child  ==0, 0, 1)),
         Victim_of_sexual_assault    =ifelse(Victim_of_sexual_assault    <0, NA,
                                             ifelse(Victim_of_sexual_assault   ==0, 0, 1))) %>%
  dplyr::select(FID, IID,
                Recent_poor_appetite_or_overeating, Maternal_smoking_around_birth, Ever_smoked, Ever_taken_cannabis,
                Physically_abused_by_family_as_a_child, Victim_of_sexual_assault, afs10) %>%
  mutate(eid=FID) %>%
  left_join(phen2, by="eid") %>%
  mutate(risk_taking    =ifelse(risk_taking    <0, NA,
                                ifelse(risk_taking   ==0, 0, 1))) %>%
  # dplyr::rename(eid=FID) %>%
  dplyr::select(FID, IID, risk_taking, ssb=ssb2, bsb, essb,
                Recent_poor_appetite_or_overeating, Maternal_smoking_around_birth, Ever_smoked, Ever_taken_cannabis,
                Physically_abused_by_family_as_a_child, Victim_of_sexual_assault, afs10, afs13) %>% 
  summarise(across(1:afs10, ~ base::mean(.x,  na.rm=T))) %>% 
  dplyr::select(-c(FID, IID)) %>% 
  pivot_longer(cols=risk_taking:afs10) %>% 
  arrange(tolower(name))

sample.prev<-c(NA, binary_prevs$value, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
population.prev<-sample.prev

# the folder of LD scores
ld <- "eur_w_ld_chr/"

# the folder of LD weights [typically the same as folder of LD scores]
wld <- "eur_w_ld_chr/"

munged_files <- paste0(trait.names, ".sumstats.gz")

# run LDSC
LDSCoutput <- ldsc(traits=munged_files,sample.prev=sample.prev,population.prev=population.prev,
                   ld=ld,wld=wld,trait.names=trait.names)

# Run SZ's model with reversed paths for Fig S14
MODEL0 <- 'NC ~ bsb+risk_taking
          bsb ~ risk_taking'

result <- usermodel(LDSCoutput, estimation = "DWLS", model = MODEL0,  std.lv = TRUE, imp_cov = TRUE)
result$results


MODEL1 <- 'NC ~ risk_taking + bsb
          risk_taking ~ bsb'

result <- usermodel(LDSCoutput, estimation = "DWLS", model = MODEL1,  std.lv = TRUE, imp_cov = TRUE)
result$results


MODEL2 <- 'bsb ~ NC
          risk_taking ~ bsb+NC'

result <- usermodel(LDSCoutput, estimation = "DWLS", model = MODEL2,  std.lv = TRUE, imp_cov = TRUE)
result$results


MODEL3 <- 'NC ~ bsb+risk_taking
          bsb ~~ risk_taking'

result <- usermodel(LDSCoutput, estimation = "DWLS", model = MODEL3,  std.lv = TRUE, imp_cov = TRUE)
result$results


# Run alt models replacing risk-taking with other measures
results_all_b <- data.frame()
for (i in c(2:3,5:31)){
  MODEL <- paste0("NC ~ ", trait.names[i], "+bsb
          bsb ~ ", trait.names[i])
  
  result <- usermodel(LDSCoutput, estimation = "DWLS", model = MODEL,  std.lv = TRUE, imp_cov = TRUE)
  
  results_all_b <- bind_rows(results_all_b, result$results)
  
  
}

results_all_b$trait <- rep(trait.names[c(2:3,5:31)], each=6)
 
# Measure X x BSB (Fig S15)
trait_b <- results_all_b %>% dplyr::filter(op=="~" & lhs=="bsb")
tb <- trait_b %>% 
  dplyr::filter(!grepl("children", rhs)) %>%
  dplyr::filter(!grepl("partners", rhs)) %>%
  dplyr::filter(!grepl("osb_any", rhs)) %>%
  dplyr::filter(!grepl("ssb", trait)) %>%
  dplyr::filter(!grepl("essb", trait)) %>%
  dplyr::filter(!grepl("afs10", trait)) %>%
  dplyr::filter(abs(STD_Genotype)<1) %>% 
  # mutate(rhs=gsub("_", " ", rhs)) %>%
  mutate(trait=gsub("_", " ", trait)) %>%
  mutate(trait=ifelse(trait=="afs13", "First had sexual intercourse before age 13", trait)) %>%
  mutate(trait=ifelse(trait=="Frequency of solariumsunlamp use", "Frequency of solarium [tanning bed]/sunlamp use", trait)) %>%
  mutate(trait=ifelse(trait=="essb", "Song & Zhang's 'Exclusively Same-sex Sexual Behavior' phenotype", trait)) %>%
  mutate(trait=ifelse(trait=="ssb", "Ever had same-sex sexual intercourse", trait)) %>%
  mutate(trait=ifelse(trait=="sb any", "Ever had sexual intercourse", trait)) %>%
  mutate(trait=ifelse(trait=="risk taking", "Risk-taking behavior", trait)) %>%
  mutate(trait=ifelse(p_value<0.05, paste0(trait, "*"), trait)) %>%
  mutate(ci_lower=STD_Genotype-1.96*STD_Genotype_SE, 
         ci_upper=STD_Genotype+1.96*STD_Genotype_SE) %>%
  # mutate(highlight=ifelse(grepl("risk taking", trait), "a", ifelse(grepl("[*]", trait), "c", "b"))) %>%
  mutate(highlight=ifelse(grepl("Risk-taking behavior", trait), "a", "b")) %>%
  ggplot(aes(x=STD_Genotype, y=reorder(trait, STD_Genotype), colour=highlight))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey60")+
  geom_point()+
  scale_colour_manual(values=c("#008a9c", "blue"))+
  geom_text(aes(label=trait), nudge_y=0.4, size=5)+
  # geom_errorbar(aes(xmax=STD_Genotype+STD_Genotype_SE, xmin=STD_Genotype-STD_Genotype_SE), width=0.25)+
  geom_errorbar(aes(xmax=ci_upper, xmin=ci_lower), width=0.25)+
  # scale_x_continuous(breaks=round(seq(-0.3, 0.4, by=0.1), digits=2), limits=c(-0.31,0.29))+
  scale_y_discrete(expand=c(0.01, 0))+
  # xlab("Standardized Genetic Correlation with Song & Zhang's 'Bisexual Sexual Behavior' phenotype (+/- SE)")+
  xlab("Partial genetic correlation \n between BSB and Measure X in males")+
  # ylab(NULL)+
  ylab("Measure X")+
  # ylab(NULL)+
  theme_classic()+
  coord_cartesian(clip = 'off') +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        legend.position="none")

# Measure X x # children (Fig S16)
trait_nc2 <- results_all_b %>% dplyr::filter(op=="~" & lhs=="NC" & rhs!="bsb")
tnc2 <- trait_nc2 %>% 
  dplyr::filter(!grepl("children", trait)) %>%
  dplyr::filter(!grepl("partners", trait)) %>%
  dplyr::filter(!grepl("osb_any", trait)) %>%
  dplyr::filter(!grepl("ssb", trait)) %>%
  dplyr::filter(!grepl("essb", trait)) %>%
  dplyr::filter(!grepl("afs10", trait)) %>%
  dplyr::filter(abs(STD_Genotype)<1) %>% 
  mutate(trait=gsub("_", " ", trait)) %>%
  # mutate(trait=ifelse(trait=="afs10", "First had sexual intercourse before age 10", trait)) %>%
  mutate(trait=ifelse(trait=="afs13", "First had sexual intercourse before age 13", trait)) %>%
  mutate(trait=ifelse(trait=="Frequency of solariumsunlamp use", "Frequency of solarium [tanning bed]/sunlamp use", trait)) %>%
  mutate(trait=ifelse(trait=="essb", "Song & Zhang's 'Exclusively Same-sex Sexual Behavior' phenotype", trait)) %>%
  mutate(trait=ifelse(trait=="ssb", "Ever had same-sex sexual intercourse", trait)) %>%
  mutate(trait=ifelse(trait=="sb any", "Ever had sexual intercourse", trait)) %>%
  mutate(trait=ifelse(trait=="risk taking", "Risk-taking behavior", trait)) %>%
  mutate(trait=ifelse(p_value<0.05, paste0(trait, "*"), trait)) %>%
  mutate(ci_lower=STD_Genotype-1.96*STD_Genotype_SE, 
         ci_upper=STD_Genotype+1.96*STD_Genotype_SE) %>%
  # mutate(highlight=ifelse(grepl("risk taking", trait), "a", ifelse(grepl("[*]", trait), "c", "b"))) %>%
  mutate(highlight=ifelse(grepl("Risk-taking behavior", trait), "a", "b")) %>%
  
  ggplot(aes(x=STD_Genotype, y=reorder(trait, STD_Genotype), colour=highlight))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey60")+
  geom_point()+
  scale_colour_manual(values=c("#008a9c", "blue"))+
  geom_text(aes(label=trait), nudge_y=0.4, size=5)+
  # geom_errorbar(aes(xmax=STD_Genotype+STD_Genotype_SE, xmin=STD_Genotype-STD_Genotype_SE), width=0.25)+
  geom_errorbar(aes(xmax=ci_upper, xmin=ci_lower), width=0.25)+
  # scale_x_continuous(breaks=round(seq(-0.3, 0.4, by=0.1), digits=2), limits=c(-0.31,0.29))+
  scale_y_discrete(expand=c(0.01, 0))+
  xlab("Partial genetic correlation \n between Measure X in males and number of children")+
  # ylab(NULL)+
  ylab("Measure X")+
  coord_cartesian(clip = 'off') +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        legend.position="none")


# BSB x # children (Fig 3b)
bsb_nc2 <- results_all_b %>% dplyr::filter(op=="~" & lhs=="NC" & rhs=="bsb")
bnc2 <- bsb_nc2 %>% 
  dplyr::filter(!grepl("children", trait)) %>%
  dplyr::filter(!grepl("partners", trait)) %>%
  dplyr::filter(!grepl("osb_any", trait)) %>%
  dplyr::filter(!grepl("ssb", trait)) %>%
  dplyr::filter(!grepl("essb", trait)) %>%
  dplyr::filter(!grepl("afs10", trait)) %>%
  dplyr::filter(abs(STD_Genotype)<1) %>% 
  mutate(trait=gsub("_", " ", trait)) %>%
  # mutate(trait=ifelse(trait=="afs10", "First had sexual intercourse before age 10", trait)) %>%
  mutate(trait=ifelse(trait=="afs13", "First had sexual intercourse before age 13", trait)) %>%
  mutate(trait=ifelse(trait=="Frequency of solariumsunlamp use", "Frequency of solarium [tanning bed]/sunlamp use", trait)) %>%
  mutate(trait=ifelse(trait=="essb", "Song & Zhang's 'Exclusively Same-sex Sexual Behavior' phenotype", trait)) %>%
  mutate(trait=ifelse(trait=="ssb", "Ever had same-sex sexual intercourse", trait)) %>%
  mutate(trait=ifelse(trait=="sb any", "Ever had sexual intercourse", trait)) %>%
  mutate(trait=ifelse(trait=="risk taking", "Risk-taking behavior", trait)) %>%
  mutate(trait=ifelse(p_value<0.05, paste0(trait, "*"), trait)) %>%
  mutate(ci_lower=STD_Genotype-1.96*STD_Genotype_SE, 
         ci_upper=STD_Genotype+1.96*STD_Genotype_SE) %>%
  # mutate(highlight=ifelse(grepl("risk taking", trait), "a", ifelse(grepl("[*]", trait), "c", "b"))) %>%
  mutate(highlight=ifelse(grepl("Risk-taking behavior", trait), "a", "b")) %>%
  ggplot(aes(x=STD_Genotype, y=reorder(trait, STD_Genotype), colour=highlight))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey60")+
  geom_point()+
  scale_colour_manual(values=c("#008a9c", "blue"))+
  geom_text(aes(label=trait), nudge_y=0.4, size=5)+
  # geom_errorbar(aes(xmax=STD_Genotype+STD_Genotype_SE, xmin=STD_Genotype-STD_Genotype_SE), width=0.25)+
  geom_errorbar(aes(xmax=ci_upper, xmin=ci_lower), width=0.25)+
  # scale_x_continuous(breaks=round(seq(-0.3, 0.4, by=0.1), digits=2), limits=c(-0.31,0.29))+
  scale_y_discrete(expand=c(0.01, 0))+
  xlab("Partial genetic correlation \n between BSB in males and number of children")+
  # ylab(NULL)+
  ylab("Measure X \n (measure for which genetic correlations are being adjusted)")+
  coord_cartesian(clip = 'off') +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=12),
        legend.position="none")

# Fig 3c
fig3c_data <- phenos_all %>% 
  dplyr::filter(genetic_sex==1 & IID %in% keep_ids$eid & !is.na(bsb) & !is.na(Age_first_had_sexual_intercourse)) %>% 
  mutate(Age_first_had_sexual_intercourse=ifelse(Age_first_had_sexual_intercourse<10, "<10",
                                                 ifelse(Age_first_had_sexual_intercourse<13, "10-12",
                                                        ifelse(Age_first_had_sexual_intercourse<16, "13-15",
                                                               ifelse(Age_first_had_sexual_intercourse<21, "16-20",
                                                                      ifelse(Age_first_had_sexual_intercourse<30, "21-29",
                                                                             ifelse(Age_first_had_sexual_intercourse<40, "30-39",
                                                                                    ifelse(Age_first_had_sexual_intercourse<50, "40-49", "50+")))))))) %>%
  group_by(Age_first_had_sexual_intercourse, bsb) %>% 
  count %>% 
  group_by(Age_first_had_sexual_intercourse) %>% 
  mutate(ntot=sum(n), prop=n/ntot) %>% 
  dplyr::filter(bsb==1) %>% 
  mutate(Age_first_had_sexual_intercourse=as.factor(Age_first_had_sexual_intercourse))

fig3c <- fig3c_data %>% 
  ggplot(aes(x=forcats::fct_rev(Age_first_had_sexual_intercourse), y=prop))+
  geom_col()+
  coord_flip()+
  geom_hline(yintercept=0.0218, linetype="dashed", colour="grey60")+
  annotate("segment", yend = 0.0218, x = 4.5, xend = 4.5, y = 0.075, colour="grey60",
           arrow = arrow(type = "closed", length = unit(0.02, "npc")))+
  annotate("text", y = 0.08, x=4.5, colour="black", hjust=0, size=6, label="Fraction of all males in sample\n classified as BSB (0.022)")+
  geom_text(aes(label=paste0("N=", ntot)), size=5, angle=0, hjust= -0.1)+
  ylab("Fraction of males in age group classified as BSB")+
  xlab("Age first had sexual intercourse (self-reported)")+
  theme_classic()+
  # scale_x_continuous(breaks=seq(0,60,by=5), expand=c(0,0))+
  scale_x_discrete(expand=c(0.1,0))+
  scale_y_continuous(expand=c(0,0),  position="right", breaks=seq(0,0.25, by=0.05), limits=c(0,0.28))+
  theme(
    axis.title.x=element_text(size=16),
    axis.text.x=element_text(size=14),
    axis.title.y=element_text(size=16),
    axis.text.y=element_text(size=14),
    legend.position="none")

