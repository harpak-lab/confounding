##########################
# Install GenomicSEM package first. Follow instructions in https://github.com/GenomicSEM/GenomicSEM
##########################

require(GenomicSEM)

# GenomicSEM functions used in this script requires two additional downloads:
# - a file of SNPs with A1, A2 and rsID used to allign alleles across trait (typically from HapMap3)
# - an (unzipped) folder containing pre-computed LD scores & weights:
# See https://github.com/GenomicSEM/GenomicSEM for instructions on where/how to download these files

# specify the HapMap3 SNP list file:
hm3 <- "w_hm3.snplist"

# specify the directory containing LD scores:
ld <- "eur_w_ld_chr/"

# specify the directory containing LD weights [typically the same as folder of LD scores]:
wld <- "eur_w_ld_chr/"

##########################
# Get file names for GWAS summary statistics used in this analysis
# these files are available via Zenodo
##########################
# GWAS sumstat file for # children (both sexes) 
sz_nc <- "ukb_wb.quant.both_n_children_p.regenie"

# get list of all male sumstat files
regenie_results_male <- list.files(".", pattern="*[.]male.*_p.regenie", full.names=T)

# combine the above into a single list of files (order is important)
trait_files <- c(sz_nc, regenie_results_male)

# specify the trait names based on the file names of the above sumstat files
trait.names<-c("NC", gsub(".*.male[_,.]", "", gsub("_p.regenie", "", basename(regenie_results_male))))

##########################
# Munge GWAS summary statistics
##########################
munge(files=trait_files, hm3=hm3, trait.names)

# munged files will be created in the same directory as the initial sumstat files with the format
# "trait.sumstats.gz"--we need a list of these files for using the LDSC function below
munged_files <- paste0(trait.names, ".sumstats.gz")

# prevalences of binary traits to be analyzed were calculated from raw UKB phenotype data 
#    name                                     value
#  1 afs10                                  0.00173
#  2 afs13                                  0.0128 
#  3 bsb                                    0.0199 
#  4 essb                                   0.0103 
#  5 Ever_smoked                            0.653  
#  6 Ever_taken_cannabis                    0.255  
#  7 Maternal_smoking_around_birth          0.302  
#  8 Physically_abused_by_family_as_a_child 0.212  
#  9 Recent_poor_appetite_or_overeating     0.137  
# 10 risk_taking                            0.350  
# 11 ssb                                    0.0430 
# 12 Victim_of_sexual_assault               0.0764 

# add these sample prevalences to a list for the traits analyzed (order is important; must align with the order of measures
# as they appear in `trait.names`

sample.prev <- c(NA, 0.00173, 0.0128, 0.0199, 0.0103, 0.653, 0.255, 0.302, 0.212, 0.137, 0.350, 0.043, 0.0764, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

population.prev <- sample.prev

##########################
# run LDSC
##########################
LDSCoutput <- ldsc(traits=munged_files, 
                   sample.prev=sample.prev,
                   population.prev=population.prev,
                   ld=ld,
                   wld=wld,
                   trait.names=trait.names)

##########################
# Fit GenomicSEM models
##########################

# Song & Zhang's postulated model
MODEL0 <- 'NC ~ risk_taking + bsb
          risk_taking ~ bsb'

result <- usermodel(LDSCoutput, estimation = "DWLS", model = MODEL0,  std.lv = TRUE, imp_cov = TRUE)
result$results

# Run Song & Zhang's model with reversed paths for Fig S14
MODEL1 <- 'NC ~ bsb+risk_taking
          bsb ~ risk_taking'

result <- usermodel(LDSCoutput, estimation = "DWLS", model = MODEL1,  std.lv = TRUE, imp_cov = TRUE)
result$results


# Run alternative models where we replace risk-taking in Song & Zhang's postulated model with each of the other measures considered
results_all_b <- data.frame()
for (i in c(2:3,5:31)){
  MODEL <- paste0("NC ~ ", trait.names[i], "+bsb
          bsb ~ ", trait.names[i])
  
  result <- usermodel(LDSCoutput, estimation = "DWLS", model = MODEL,  std.lv = TRUE, imp_cov = TRUE)
  
  results_all_b <- bind_rows(results_all_b, result$results)
  
  
}

# assign trait name to table of results
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
fig3b_data <- results_all_b %>% dplyr::filter(op=="~" & lhs=="NC" & rhs=="bsb")
bnc2 <- fig3b_data %>% 
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
  mutate(highlight=ifelse(grepl("Risk-taking behavior", trait), "a", "b")) %>%
  ggplot(aes(x=STD_Genotype, y=reorder(trait, STD_Genotype), colour=highlight))+
  geom_vline(xintercept=0, linetype="dashed", colour="grey60")+
  geom_point()+
  scale_colour_manual(values=c("#008a9c", "blue"))+
  geom_text(aes(label=trait), nudge_y=0.4, size=5)+
  geom_errorbar(aes(xmax=ci_upper, xmin=ci_lower), width=0.25)+
  scale_y_discrete(expand=c(0.01, 0))+
  xlab("Partial genetic correlation \n between BSB in males and number of children")+
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
# note that this figure is generated using individual phenotype data from UK Biobank.
# Summary data for replicating this figure is provided in the file fig3c_data.tsv

fig3c_data <- read_tsv("fig3c_data.tsv")

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

