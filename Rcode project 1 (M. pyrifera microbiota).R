#Load packages
library(microDecon)
library(readxl)
library(readr)
library(vegan)
library(seqinr)
library(dplyr)
library(phyloseq)
###Read file to work
taxa.silva <- read_excel("taxa.silva.xlsx")
###Clean data
clean.16S <- decon(data = as.data.frame(taxa.silva),numb.blanks=,numb.ind=66,taxa=F) #<- Adjust according to the data
clean.16S<-clean.16S$decon.table
clean.16S<-clean.16S[,-2]
dim(clean.16S)

#Data Preparation for Analysis
namesTax <-c("OTU_ID")
Tax <- clean.16S[,namesTax]
Bac2 <- t(clean.16S[,!colnames(clean.16S)%in%namesTax])
colnames(Bac2)<-Tax
Bac2 <- Bac2[order(row.names(Bac2)),]
dim(Bac2)

#Data Filtering
N <- apply(I(Bac2>0), 2, sum)
range(N)
rarefaction <- Bac2[, N>0] 
dim(rarefaction)

#Singleton Removal and Additional Filter
Bac2 <- Bac2[, N>1]
dim(Bac2)
df <- Bac2[,apply( Bac2, 2, function(x) sum( x>0)) > 1] 
dim(df)

#Rarefaction curves
S <- specnumber(rarefaction)
(raremax <- min(rowSums(rarefaction)))
Srare <- rarefy(rarefaction, raremax)
plot(S, Srare, xlab = "Observed No. of ASV", ylab = "Rarefied No. of ASV")
abline(0, 1)
rarecurve(rarefaction, step = 20, col = "#404788FF", main="Rarefaction Curves 16S", ylab="ASV Observed", cex=0.8)

names=colnames(df)
#Read file taxa.silva.txt 
table<-read.table("taxa.silva.txt", header = T, sep="\t")
table<-table[,1:2]
names=colnames(df)
#Filter table to keep only rows with ASVs present in df
filtered_data <- filter(table, ASV %in% names)
write.fasta(sequences =as.list(filtered_data[,1]),names= filtered_data[,2],file.out = "silva.fasta")
##Read FASTA silva.fasta
silva<-read.fasta(file="silva.fasta", as.string = TRUE, forceDNAtolower = TRUE)

silva.2<-as.SeqFastadna(silva)
is.SeqFastadna(silva.2)
summary(silva.2[[2]])

#Make relative
library(funrar)
##Convert df to relative values
rel_table<-make_relative(df)
rel_table<-t(rel_table)
#Save file
write.csv(rel_table,file="rel_table16S")
df<-t(df)
#Save file
write.csv(df,file="table16S")

#Preparation for Further Analysis
table<-read.table("taxa.silva.txt", header = T, sep="\t")
Reads <- read_xlsx("table16S.xlsx")
colnames(Reads)
Reads<-Reads[,-1]
tax<-table[,c(xx, xx:xx)]

d <- read_xlsx("table16S.xlsx")
dim(d)
#Dataframe containing the environmental variables
env<-read.csv("Env.csv")
#Select the specific columns from the dataframe d to store in tax
tax<-d[,c(xx,xx:xx)]
##Select the columns from dataframe d with abundance readings for each ASV
reads<-d[,c(xx:xx)]
colnames(reads)

#Assign the reads object (ASV abundance table) to the asv_mat variable
asv_mat<- reads
#Assign the tax object (taxonomy table) to the tax_mat variable
tax_mat<- tax
#Assign the env object (environmental or sample information) to the samples_df variable
samples_df <- env

#convert to dataframe
asv_mat <- as.data.frame(asv_mat)
tax_mat <- as.data.frame(tax_mat)

# Setting row names and converting data frames to matrices
row.names(asv_mat) <- asv_mat$ASV
asv_mat <- asv_mat%>% select (-ASV) 
row.names(tax_mat) <- tax_mat$ASV

tax_mat <- tax_mat%>% select (-ASV) 
row.names(samples_df) <- samples_df$ID
samples_df <- samples_df %>% select (-ID) 

asv_mat <- as.matrix(asv_mat)
tax_mat <- as.matrix(tax_mat)
#Creating Objects for phyloseq
ASV = otu_table(asv_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
#Create phyloseq Object
Bac_16S <- phyloseq(ASV, TAX, samples)
#Subset
Bac_16S<-Bac_16S %>% subset_taxa(Kingdom !="Archaea")     
Bac_16S<-Bac_16S %>% subset_taxa(Order != "Chloroplast")     
Bac_16S<-Bac_16S %>% subset_taxa(Family !="Mitochondria")     

#############-------------Alpha diversty-----------------####################
library(ggplot2)
library(phyloseq)
library(dplyr)
library(patchwork)
#Subset data
subset<-subset_samples(Bac_16S,  Type != "xx"  ) # <- Exclude sea water
sample_sums(subset)
#Rarefy data 
rarefied<-rarefy_even_depth(subset, sample.size = min(sample_sums(subset)),
                            replace = FALSE,verbose = TRUE, rngseed = TRUE)
sample_sums(rarefied)
otu_table(rarefied)

#calculated alpha diversity (Shannon)) to make basic plot with ggplot
p<-plot_richness(rarefied, x="Type", measures=c("Shannon"))+
  geom_boxplot(alpha=1)+
  theme(strip.text.x= element_text(size=18),
        axis.text.y.left = element_text(size = 18),
        axis.title.x = element_text(size = 12),
        strip.text.y = element_text(size=18),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15))+ xlab("")+ 
  facet_grid(Date~factor(Site, levels=c('Bahia_Buzo', 'San_Gregorio', 'Buque_Quemado')))+
  stat_compare_means(method="wilcox.test",label.y.npc = 0.95, label.x = 1.3, position = position_nudge(y = 0.2), size = 5)+
  scale_y_continuous(limits = c(2, 6.5), expand = expansion(mult = c(0, 0.1)))

newSTorder = c("Apical","Sporophylls")
p$data$Type <- as.character(p$data$Type)
p$data$Type <- factor(p$data$Type, levels=newSTorder)
print(p)

#Sort by Type
p$data$Type <- factor(p$data$Type, levels = c("Apical", "Sporophylls"))

p1<-plot_richness(rarefied, x = "Date", measures = c("Shannon")) +
  geom_boxplot(aes(fill = Type), alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~Type, nrow = 1) +
  scale_fill_manual(values = c("Apical" = "#C0FF3E", "Sporophylls" = "#698B22")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size = 18),
    strip.text = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.position = "none"
  ) +
  ylab("Alpha diversity measure")+
  stat_compare_means(method="wilcox.test",label.x.npc = 0,position = position_nudge(y = 0.15))+
  labs(title = "")+
  theme(plot.title = element_text(size = 20, face = "bold"))
print(p1)

#Sort by Site
my_comparisons <- list(c("Bahia_Buzo", "San_Gregorio"), 
                       c("San_Gregorio", "Buque_Quemado"),
                       c("Bahia_Buzo", "Buque_Quemado"))
p2<-plot_richness(rarefied, x="Site", measures=c("Shannon"))+
  geom_boxplot(alpha=1)+
  theme(strip.text.x= element_text(size=18),
        axis.text.y.left = element_text(size = 18),
        axis.title.x = element_text(size = 12),
        strip.text.y = element_text(size=18),
        axis.text.x = element_blank(),,
        axis.title.y = element_blank(),)+ xlab("")+ 
  facet_grid(Date~factor(Type, levels=c('Apical', 'Sporophylls')))+
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test",
                     label = "p.signif",
                     label.y.npc = c(0.95, 0.88, 0.81),
                     vjust = 0.1,
                     size = 5) + 
  scale_y_continuous(limits = c(0, 8), expand = expansion(mult = c(0, 0.15)))

print(p2)

sample_data(rarefied)$Site <- factor(sample_data(rarefied)$Site, levels = c("Bahia_Buzo", "San_Gregorio", "Buque_Quemado"))
p3<- plot_richness(rarefied, x = "Date", measures = c("Shannon")) +
  geom_boxplot(aes(fill = Site), alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~Site, nrow = 1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    strip.text = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.position = "none"
  ) +
  ylab("")+
  stat_compare_means(method = "wilcox.test",label.y.npc = 1) +
  labs(title = "")+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  scale_y_continuous(limits = c(2, 6.5), expand = expansion(mult = c(0, 0.15)))
print(p3)
#Chao1
p4<-plot_richness(rarefied, x="Type", measures=c("Chao1"))+
  geom_boxplot(alpha=1)+
  theme(strip.text.x= element_text(size=18),
        axis.text.y.left = element_text(size = 18),
        axis.title.x = element_text(size = 12),
        strip.text.y = element_text(size=18),
        axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 15))+ xlab("")+ 
  facet_grid(Date~factor(Site, levels=c('Bahia_Buzo', 'San_Gregorio', 'Buque_Quemado')))+
  stat_compare_means(method="wilcox.test",label.y.npc = 0.95,label.x = 1.3,position = position_nudge(y = 0.15), size = 5)+
  scale_y_continuous(limits = c(0, 1800), expand = expansion(mult = c(0, 0.1)))

newSTorder = c("Apical","Sporophylls")
p$data$Type <- as.character(p$data$Type)
p$data$Type <- factor(p$data$Type, levels=newSTorder)
print(p4)

#Sort by Type
p$data$Type <- factor(p$data$Type, levels = c("Apical", "Sporophylls"))

p5<-plot_richness(rarefied, x = "Date", measures = c("Chao1")) +
  geom_boxplot(aes(fill = Type), alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~Type, nrow = 1) +
  scale_fill_manual(values = c("Apical" = "#C0FF3E", "Sporophylls" = "#698B22")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,size = 18),
    strip.text = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.position = "none"
  ) +
  ylab("Alpha diversity measure")+
  stat_compare_means(method="wilcox.test",label.x.npc = 0,position = position_nudge(y = 0.15))+
  labs(title = "")+
  theme(plot.title = element_text(size = 20, face = "bold"))
print(p5)

#Sort by Site
my_comparisons <- list(c("Bahia_Buzo", "San_Gregorio"), 
                       c("San_Gregorio", "Buque_Quemado"),
                       c("Bahia_Buzo", "Buque_Quemado"))
p6<-plot_richness(rarefied, x="Site", measures=c("Chao1"))+
  geom_boxplot(alpha=1)+
  theme(strip.text.x= element_text(size=18),
        axis.text.y.left = element_text(size = 15),
        axis.title.x = element_text(size = 12),
        strip.text.y = element_text(size=18),
        axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_blank())+ xlab("")+ 
  facet_grid(Date~factor(Type, levels=c('Apical', 'Sporophylls')))+
  stat_compare_means(comparisons = my_comparisons, 
                     method = "wilcox.test",
                     label = "p.signif",
                     label.y.npc = c(0.95, 0.88, 0.81),
                     , size = 5) + 
  scale_y_continuous(limits = c(0, 2000), expand = expansion(mult = c(0, 0.15)))

print(p6)

unique(p$data$variable)
sample_data(rarefied)$Site <- factor(sample_data(rarefied)$Site, levels = c("Bahia_Buzo", "San_Gregorio", "Buque_Quemado"))
p7<- plot_richness(rarefied, x = "Date", measures = c("Chao1")) +
  geom_boxplot(aes(fill = Site), alpha = 0.7, outlier.shape = NA) +
  facet_wrap(~Site, nrow = 1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1),
    strip.text = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.position = "none"
  ) +
  ylab("")+
  stat_compare_means(method = "wilcox.test",label.y.npc = 1) +
  labs(title = "")+
  theme(plot.title = element_text(size = 20, face = "bold"))+
  scale_y_continuous(limits = c(2, 1600), expand = expansion(mult = c(0, 0.15)))
print(p7)
#Figura Sup.
p1 <- p1 + theme(axis.text.x = element_blank())
p3 <- p3 + theme(axis.text.x = element_blank())

p5 <- p5 + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.03, 0.95))  # (x, y) entre 0 y 1)
p7 <- p7 + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.01, 0.95)) 
p1 <- p1 + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.03, 0.98)) 
p3 <- p3 + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.01, 0.98)) 

final_plot <- (p1 | p3)/ (p5 | p7)
print(final_plot)
#Fig. paper
p <- p + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.03, 0.95))  # (x, y) entre 0 y 1)
p2 <- p2 + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.01, 0.95)) 
p4 <- p4 + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.03, 0.98)) 
p6 <- p6 + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.01, 0.98)) 
plot <- (p | p2) / (p4 | p6)
plot
#############-------------Figures-----------------####################
###Phylum----
phylumabundance <- Bac_16S %>%
  tax_glom(taxrank = "Phylum")%>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Phylum) 
head(phylumabundance)

all <- phylumabundance %>%
  filter( Type != "Sea_Water") %>%
  select(Phylum, Type, Date, Abundance) %>%
  group_by(Phylum, Type, Date) %>%
  summarize(avg_abundance = mean(Abundance), .groups = "drop") %>%
  filter(avg_abundance >= 0.05)

head(all)

my_theme <- theme(
  axis.text.y.left = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.x = element_text(size = 18, vjust = 1, hjust = 1, angle = 45),
  title = element_text(size = 16),
  legend.text = element_text(size = 16),
  strip.text.x = element_text(size = 15)
)
colors=c(Actinobacteriota="#E7EBFA", Bacteroidota="#8470FF",Campylobacterota="pink",	
         Planctomycetota="#EE7600",	Proteobacteria="#CD1076",Verrucomicrobiota="#FF0000")

p5<-ggplot(all) +
  geom_col(aes(
    x = factor(Date, levels = c("August", "March")),
    y = avg_abundance,
    fill = Phylum
  ),
  position = "fill",
  show.legend = FALSE,
  width = 0.8,
  colour = "black"
  ) +
  ylab("Relative abundance (%)") +
  xlab("") +
  facet_grid(. ~ factor(Type, levels = c("xx", "xx"))) +#<- Sample type name
  theme(
    axis.text.y.left = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, vjust = 1, hjust = 1, angle = 45),
    title = element_text(size = 12),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 15)
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors, drop = FALSE) +
  my_theme

p5<-ggplot(all) +
  geom_col(aes(
    x = factor(Date, levels = c("August", "March")),
    y = avg_abundance,
    fill = Phylum
  ),
  position = "fill",
  show.legend = FALSE,
  width = 0.8,
  colour = "black"
  ) +
  ylab("Relative abundance (%)") +
  xlab("") +
  facet_grid(. ~ factor(Type, levels = c("xx", "xx"))) +#<- Sample type name
  theme(
    axis.text.y.left = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, vjust = 1, hjust = 1, angle = 45),
    title = element_text(size = 12),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 15)
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors, drop = FALSE) +
  my_theme

all_sitio <- phylumabundance %>%
  filter(Type != "Sea_Water") %>%
  select(Phylum, Site, Date, Abundance) %>%
  group_by(Phylum, Site, Date) %>%
  summarize(avg_abundance = mean(Abundance), .groups = "drop") %>%
  filter(avg_abundance >= 0.05)
p6<-ggplot(filter(all_sitio)) +
  geom_col(mapping = aes(
    x = factor(Date, levels = c("August", "March")),
    y = avg_abundance,
    fill = Phylum
  ),
  position = "fill",
  show.legend = TRUE,
  width = 0.8,
  colour = "black"
  ) +
  ylab("") +
  xlab("") +
  facet_grid(. ~ factor(Site, levels = c("xx", "xx", "xx"))) + #<- Sample site name
  theme(
    axis.text.y.left = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 18, vjust = 1, hjust = 1, angle = 45),
    title = element_text(size = 12),
    legend.text = element_text(size = 18),
    strip.text.x = element_text(size = 15)
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors) +
  my_theme
p5 <- p5 + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold"),plot.tag.position = c(0.03, 0.98))
p6<- p6 + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.03, 0.98))
filo_plot <- (p5 | p6) 
filo_plot

#Phylum water
all_water <- phylumabundance %>%
  filter( Type != "Apical", Type != "Sporophylls") %>%
  select(Phylum, Type, Date, Abundance) %>%
  group_by(Phylum, Type, Date) %>%
  summarize(avg_abundance = mean(Abundance), .groups = "drop") %>%
  filter(avg_abundance >= 0.05)

head(all_water)

my_theme <- theme(
  axis.text.y.left = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.x = element_text(size = 18, vjust = 1, hjust = 1, angle = 45),
  title = element_text(size = 16),
  legend.text = element_text(size = 16),
  strip.text.x = element_text(size = 15)
)

colors=c(Actinobacteriota="#E7EBFA", Bacteroidota="#8470FF",Crenarchaeota="pink",Cyanobacteria="#00CED1",	
         Planctomycetota="#EE7600",	Proteobacteria="#CD1076",Verrucomicrobiota="#FF0000")

# Gráfico final
w5<-ggplot(all_water) +
  geom_col(aes(
    x = factor(Date, levels = c("xx", "xx")),#<- Sample date name
    y = avg_abundance,
    fill = Phylum
  ),
  position = "fill",
  show.legend = FALSE,
  width = 0.8,
  colour = "black"
  ) +
  ylab("Relative abundance (%)") +
  xlab("") +
  facet_grid(. ~ factor(Type, levels = c( "Sea_Water"))) +
  theme(
    axis.text.y.left = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, vjust = 1, hjust = 1, angle = 45),
    title = element_text(size = 12),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 15)
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors, drop = FALSE) +
  my_theme

all_sitio_w <- phylumabundance %>%
  filter(Type != "Apical", Type != "Sporophylls") %>%
  select(Phylum, Site, Date, Abundance) %>%
  group_by(Phylum, Site, Date) %>%
  summarize(avg_abundance = mean(Abundance), .groups = "drop") %>%
  filter(avg_abundance >= 0.05)
w6<-ggplot(filter(all_sitio_w)) +
  geom_col(mapping = aes(
    x = factor(Date, levels = c("xx", "xx")),
    y = avg_abundance,
    fill = Phylum
  ),
  position = "fill",
  show.legend = TRUE,
  width = 0.8,
  colour = "black"
  ) +
  ylab("") +
  xlab("") +
  facet_grid(. ~ factor(Site, levels = c("xx", "xx", "xx"))) +  # facet por sitio
  theme(
    axis.text.y.left = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 18, vjust = 1, hjust = 1, angle = 45),
    title = element_text(size = 12),
    legend.text = element_text(size = 18),
    strip.text.x = element_text(size = 15)
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors) +
  my_theme
p5 <- p5 + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold"),plot.tag.position = c(0.03, 0.98))
p6<- p6 + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.03, 0.98))
w5 <- w5 + labs(tag = "C") + theme(plot.tag = element_text(size = 20, face = "bold"),plot.tag.position = c(0.03, 0.98))
w6<- w6 + labs(tag = "D") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.03, 0.98))

WSP<-(p5 | p6) / (w5 | w6)
WSP
###Class----

abundance_class <- rarefied %>% # <- Same flow applicable for Order
  tax_glom(taxrank = "Class")%>%                     # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Class) 
head(abundance_class)


all_class <- abundance_class %>%
  filter(Type != "Sea_Water") %>%
  select(Class,Type,Date,Abundance) %>%
  group_by(Type, Date, Class) %>%
  summarize(
    avg_abundance = mean(Abundance)) %>%
  filter(avg_abundance >= 0.05)
head(all_class)
write.csv(all_class,file="class.csv")

unique(all_class$Class)

colors4 <- c(Cyanobacteriia="#00CED1" ,Bacteroidia="navy", Verrucomicrobiae="#FF0000" ,Planctomycetes="#CAACCB",Campylobacteria="pink", 
             Gammaproteobacteria="#882E72",Nitrososphaeria="#777777",Alphaproteobacteria="#008B00", 
             Actinobacteria="lightgoldenrod1" )

p7<-ggplot(all_class) +
  geom_col(aes(
    x = factor(Date, levels = c("xx", "xx")),#<-  date name
    y = avg_abundance,
    fill = Class
  ),
  position = "fill",
  show.legend = FALSE,
  width = 0.8,
  colour = "black"
  ) +
  ylab("Relative abundance (%)") +
  xlab("") +
  facet_grid(. ~ factor(Type, levels = c("xx", "xx"))) +#<- Sample type name
  my_theme+
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors4, drop = FALSE) 


all_class_sitio <- abundance_class %>%
  filter(Type != "Sea_Water") %>%
  select(Class,Site,Date,Abundance) %>%
  group_by(Site, Date, Class) %>%
  summarize(
    avg_abundance = mean(Abundance)) %>%
  filter(avg_abundance >= 0.05)
head(all_class)

p8<-ggplot(filter(all_class_sitio )) +
  geom_col(mapping = aes(
    x = factor(Date, levels = c("August", "March")),
    y = avg_abundance,
    fill = Class
  ),
  position = "fill",
  show.legend = TRUE,
  width = 0.8,
  colour = "black"
  ) +
  ylab("") +
  xlab("") +
  facet_grid(. ~ factor(Site, levels = c("xx", "xx", "xx"))) + #<- Sample site name
  my_theme +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors4)

p7 <- p7 + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold"),plot.tag.position = c(0.03, 0.98))
p8<- p8 + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.03, 0.98))
class_plot <- (p7 | p8) 
class_plot


#AMP figure----
library(ampvis2)
obj<-rarefied
otutable <- data.frame(phyloseq::otu_table(obj)@.Data,
                       phyloseq::tax_table(obj)@.Data,
                       check.names = FALSE)

otutable$Species = otutable$Genus

av2 <- amp_load(otutable, env)
av2_est <- subset_samples(rarefied, Type == "xx" | Type == "xx")#<- Sample type name
obj3<-av2_est
otutable3 <- data.frame(phyloseq::otu_table(obj3)@.Data,
                        phyloseq::tax_table(obj3)@.Data,
                        check.names = FALSE
)

otutable3$Species = otutable3$Genus

av2_est <- amp_load(otutable3, env)
av2_est$tax
amp_heatmap(av2_est,tax_show = 15,
            group_by = "Type", facet_by = c("Site", "Date"),
            tax_aggregate = "Genus",tax_add = "Order",
            color_vector = c("#5666B6","whitesmoke","#B2182B"))+ facet_grid(Date~factor(Site, levels=c('xx', 'xx', 'xx')))+#<- Sample site name
  theme(axis.text.y = element_text(size = 10),title = element_text(size = 12)) +ggtitle("")+
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1))+theme(strip.text.x = element_text(size = 14))

env$Type <- factor(env$Type, levels = c("xx", "xx"))#<- Sample type name
av2_est <- amp_load(otutable3, env)

p11<-amp_heatmap(av2_est,
                 tax_show = 15,
                 group_by = "Date",
                 facet_by = "Type",
                 tax_aggregate = "Genus",
                 tax_add = "Order",
                 color_vector = c("#5666B6", "whitesmoke", "#B2182B")) +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust = 1),
    strip.text.x = element_text(size = 14),
    title = element_text(size = 12)
  )

env$Site <- factor(env$Site, levels = c('xx', 'xx', 'xx'))#<- Sample site name
av2_site <- amp_load(otutable3, env)

p12<-amp_heatmap(av2_site,
                 tax_show = 15,
                 group_by = "Date",
                 facet_by = "Site",
                 tax_aggregate = "Genus",
                 tax_add = "Order",
                 color_vector = c("#5666B6", "whitesmoke", "#B2182B")) +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust = 1),
    strip.text.x = element_text(size = 14),
    title = element_text(size = 12)
  )

library(patchwork)
p11 <- p11 + labs(tag = "A") + theme(plot.tag = element_text(size = 20, face = "bold"),plot.tag.position = c(0.03, 0.98))
p12<- p12 + labs(tag = "B") + theme(plot.tag = element_text(size = 20, face = "bold"), plot.tag.position = c(0.03, 0.98))
amp_plot <- (p11 | p12) 
amp_plot

#############-------------Beta diversty-----------------####################
#ANOSIM-NMDS 
#-----Site
library(phyloseq)
library(vegan)

# filter
subset <- subset_samples(Bac_16S, Type != "Sea_Water")
sample_sums(subset)

# Rarefaction
set.seed(123) 
subset_rarefied <- rarefy_even_depth(subset, rngseed = 123, verbose = FALSE)
sample_sums(subset_rarefied)

# matrix
otu_mat <- as(otu_table(subset_rarefied), "matrix")

if (taxa_are_rows(subset_rarefied)) {
  otu_mat <- t(otu_mat)
}

# Hellinger
otu_hellinger <- decostand(otu_mat, method = "hellinger") # Select method

# Bray-Curtis post-Hellinger
dist_hell_bray <- vegdist(otu_hellinger, method = "bray")

meta_df <- as(sample_data(subset_rarefied), "data.frame")

env_august <- meta_df[meta_df$Date == "August", ] # <- Date of interest
env_march <- meta_df[meta_df$Date == "March", ]
env_august$ID <- rownames(env_august)
env_march$ID <- rownames(env_march)

# Subset August
samples_august <- env_august$ID

mat_dist_hell_bray <- as.matrix(dist_hell_bray)
sub_mat_august_bray <- mat_dist_hell_bray[samples_august, samples_august]
dist_august_bray <- as.dist(sub_mat_august_bray)

# Subset March
samples_march <- env_march$ID

mat_dist_hell_bray_m <- as.matrix(dist_hell_bray)
sub_mat_march_bray <- mat_dist_hell_bray_m[samples_march, samples_march]
dist_march_bray <- as.dist(sub_mat_march_bray)

# ANOSIM Site-Agosto 
anosim_august_site_bray <- anosim(dist_august_bray, env_august$Site, permutations = 9999)

# ANOSIM Site-Marzo 
anosim_march_site_bray <- anosim(dist_march_bray, env_march$Site, permutations = 9999)

# NMDS
nmds_august_bray <- metaMDS(dist_august_bray, k = 2, trymax = 100)
nmds_march_bray  <- metaMDS(dist_march_bray, k = 2, trymax = 100)

site_colors <- c("xx" = "#1b9e77", "xx" = "#d95f02", "xx" = "#7570b3")
site_shapes <- c("xx" = 17, "xx" = 19, "xx" = 15)

r_august_site_bray <- round(anosim_august_site_bray$statistic, 3)
p_august_site_bray <- format.pval(anosim_august_site_bray$signif, digits = 3, eps = .001)

r_march_site_bray <- round(anosim_march_site_bray$statistic, 3)
p_march_site_bray <- format.pval(anosim_march_site_bray$signif, digits = 3, eps = .001)

# --- NMDS August-Site ---
plot(nmds_august_bray$points, type = "n", main = "xx",
     xlab = "NMDS1", ylab = "NMDS2")

for (s in unique(env_august$Site)) {
  idx <- which(env_august$Site == s)
  points(nmds_august_bray$points[idx, ], 
         col = site_colors[s], pch = site_shapes[s], cex = 1.2)
}

legend(x = -0.6, y = -0.25, legend = names(site_colors), 
       col = site_colors, pch = site_shapes, bty = "n")

text(x = -0.6, y = 0.25, 
     labels = bquote(R == .(r_august_site_bray) ~ ", " ~ p == .(p_august_site_bray)),
     cex = 1, font = 2, pos = 4)

# --- NMDS March-Site ---
plot(nmds_march_bray$points, type = "n", main = "xx",
     xlab = "NMDS1", ylab = "NMDS2")

for (s in unique(env_march$Site)) {
  idx <- which(env_march$Site == s)
  points(nmds_march_bray$points[idx, ], 
         col = site_colors[s], pch = site_shapes[s], cex = 1.2)
}

legend("topright", legend = names(site_colors), 
       col = site_colors, pch = site_shapes, bty = "n")

text(x = 0.4, y = -0.3, 
     labels = bquote(R == .(r_march_site_bray) ~ ", " ~ p == .(p_march_site_bray)),
     cex = 1, font = 2, pos = 3)

#-----Type
dist_hell_bray <- vegdist(otu_hellinger, method = "bray")

mat_dist_hell_bray <- as.matrix(dist_hell_bray)

samples_august <- env_august$ID  
sub_mat_august <- mat_dist_hell_bray[samples_august, samples_august]

dist_august_bray <- as.dist(sub_mat_august)

samples_march <- env_march$ID
sub_mat_march <- mat_dist_hell_bray[samples_march, samples_march]
dist_march_bray <- as.dist(sub_mat_march)

# ANOSIM
anosim_august_bray_type <- anosim(dist_august_bray, env_august$Type, permutations = 9999)
anosim_march_bray_type  <- anosim(dist_march_bray,  env_march$Type,  permutations = 9999)

r_august_bray_type <- round(anosim_august_bray_type$statistic, 3)
p_august_bray_type <- format.pval(anosim_august_bray_type$signif, digits = 3, eps = .001)

r_march_bray_type <- round(anosim_march_bray_type$statistic, 3)
p_march_bray_type <- format.pval(anosim_march_bray_type$signif, digits = 3, eps = .001)

# NMDS
nmds_august_bray_type <- metaMDS(dist_august_bray, k = 2, trymax = 100)
nmds_march_bray_type  <- metaMDS(dist_march_bray, k = 2, trymax = 100)

type_colors <- c("xx" = "#FF8C00", "xx" = "#8B4500")#<- Sample type name
type_shapes <- c("xx" = 19, "xx" = 15)

# --- NMDS August ---
plot(nmds_august_bray_type$points, type = "n", main = "xx",
     xlab = "NMDS1", ylab = "NMDS2")

for (t in unique(env_august$Type)) {
  idx <- which(env_august$Type == t)
  points(nmds_august_bray_type$points[idx, ], 
         col = type_colors[t], pch = type_shapes[t], cex = 1.2)
}

legend("topright", legend = names(type_colors), 
       col = type_colors, pch = type_shapes, bty = "n")

text(x = -0.6, y = -0.3, 
     labels = bquote(R == .(r_august_bray_type) ~ ", " ~ p == .(p_august_bray_type)),
     cex = 1, font = 2, pos = 4)

# --- NMDS March ---
plot(nmds_march_bray_type$points, type = "n", main = "xx",
     xlab = "NMDS1", ylab = "NMDS2")

for (t in unique(env_march$Type)) {
  idx <- which(env_march$Type == t)
  points(nmds_march_bray_type$points[idx, ], 
         col = type_colors[t], pch = type_shapes[t], cex = 1.2)
}

legend("topright", legend = names(type_colors), 
       col = type_colors, pch = type_shapes, bty = "n")

text(x = 0.35, y = -0.2, 
     labels = bquote(R == .(r_march_bray_type) ~ ", " ~ p == .(p_march_bray_type)),
     cex = 1, font = 2, pos = 3)

# Configurar layout de 2x2
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
# ------------------- A. NMDS August por Site (Bray) -------------------
plot(nmds_august_bray$points, type = "n", 
     main = "xx",
     xlab = "NMDS1", ylab = "NMDS2",
     cex.main = 1.5, cex.lab = 1.4, cex.axis = 1.2)

for (s in unique(env_august$Site)) {
  idx <- which(env_august$Site == s)
  points(nmds_august_bray$points[idx, ], 
         col = site_colors[s], pch = site_shapes[s], cex = 1.2)
}
legend(x = -0.6, y = 0.1, legend = names(site_colors), 
       col = site_colors, pch = site_shapes, bty = "n", cex = 1.3)

text(x = 0, y = 0.4, 
     labels = bquote(R == .(r_august_site_bray) ~ ", " ~ p == .(p_august_site_bray)),
     cex = 1.7, font = 2, pos = 4)

mtext("A", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)

# ------------------- B. NMDS March por Site (Bray) -------------------
plot(nmds_march_bray$points, type = "n", 
     main = "xx",
     xlab = "NMDS1", ylab = "NMDS2",
     cex.main = 1.5, cex.lab = 1.4, cex.axis = 1.2)

for (s in unique(env_march$Site)) {
  idx <- which(env_march$Site == s)
  points(nmds_march_bray$points[idx, ], 
         col = site_colors[s], pch = site_shapes[s], cex = 1.2)
}
legend(x = -0.5, y = -0.9, legend = names(site_colors), 
       col = site_colors, pch = site_shapes, bty = "n")

text(x = 0.3, y = 0.22, 
     labels = bquote(R == .(r_march_site_bray) ~ ", " ~ p == .(p_march_site_bray)),
     cex = 1.7, font = 2, pos = 3)

mtext("B", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)
# ------------------- C. NMDS August por Type (Bray) -------------------
plot(nmds_august_bray_type$points, type = "n", 
     main = "xx",
     xlab = "NMDS1", ylab = "NMDS2",
     cex.main = 1.5, cex.lab = 1.4, cex.axis = 1.2)

for (t in unique(env_august$Type)) {
  idx <- which(env_august$Type == t)
  points(nmds_august_bray_type$points[idx, ], 
         col = type_colors[t], pch = type_shapes[t], cex = 1.2)
}
legend(x = -0.6, y = 0, legend = names(type_colors), 
       col = type_colors, pch = type_shapes, bty = "n", cex = 1.3)

text(x = 0.01, y = 0.4, 
     labels = bquote(R == .(r_august_bray_type) ~ ", " ~ p == .(p_august_bray_type)),
     cex = 1.7, font = 2, pos = 4)

mtext("C", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)
# ------------------- D. NMDS March por Type (Bray) -------------------
plot(nmds_march_bray_type$points, type = "n", 
     main = "xx",
     xlab = "NMDS1", ylab = "NMDS2",
     cex.main = 1.5, cex.lab = 1.4, cex.axis = 1.2)

for (t in unique(env_march$Type)) {
  idx <- which(env_march$Type == t)
  points(nmds_march_bray_type$points[idx, ], 
         col = type_colors[t], pch = type_shapes[t], cex = 1.2)
}
legend(x = -0.8, y = 0.9, legend = names(type_colors), 
       col = type_colors, pch = type_shapes, bty = "n")

text(x = 0.35, y = -0.25, 
     labels = bquote(R == .(r_march_bray_type) ~ ", " ~ p == .(p_march_bray_type)),
     cex = 1.7, font = 2, pos = 3)

mtext("D", side = 3, line = 1, adj = 0, cex = 1.5, font = 2)

#############-------------GLM-----------------####################
#GLM Genus blades--------
library(betareg)
library(emmeans)
library(statmod)

av2_est <- subset_samples(rarefied, Type == "xx" | Type == "xx")#<- Sample type name
genusabundance <- av2_est %>%
  tax_glom(taxrank = "Genus")%>%                     # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) 
head(genusabundance)

test <- genusabundance %>%
  select(Genus,Site,Type,Date,Sample, Abundance) %>%
  group_by(Genus, Site,Type, Date,Sample) %>%
  summarize(
    avg_abundance = mean(Abundance)) 
#Run the models for the most abundant genera according to to the Heat map
Psychromonas<-test[test$Genus=="Psychromonas",]
Psychromonas$Type <- factor(Psychromonas$Type)
Psychromonas$Site <- factor(Psychromonas$Site)
Psychromonas$Date <- factor(Psychromonas$Date)
Mod1 <- betareg(avg_abundance  ~Type+Site+Date, data=Psychromonas, na.action=na.omit, link = "logit")
AIC(Mod1)
summary(emmeans(Mod1, pairwise ~ Type, adjust="holm"))
summary(emmeans(Mod1, pairwise ~ Site, adjust="holm"))
summary(emmeans(Mod1, pairwise ~ Date, adjust="holm"))

Persicirhabdus<-test[test$Genus=="Persicirhabdus",]
Persicirhabdus$Type <- factor(Persicirhabdus$Type)
Persicirhabdus$Site <- factor(Persicirhabdus$Site)
Persicirhabdus$Date <- factor(Persicirhabdus$Date)
Mod2 <- betareg(avg_abundance  ~Type+Site+Date, data=Persicirhabdus, na.action=na.omit, link = "logit")
AIC(Mod2)
summary(emmeans(Mod2, pairwise ~ Type, adjust="holm"))
summary(emmeans(Mod2, pairwise ~ Site, adjust="holm"))
summary(emmeans(Mod2, pairwise ~ Date, adjust="holm"))

Rubritalea<-test[test$Genus=="Rubritalea",]
Rubritalea$Type <- factor(Rubritalea$Type)
Rubritalea$Site <- factor(Rubritalea$Site)
Rubritalea$Date <- factor(Rubritalea$Date)
Mod3 <- betareg(avg_abundance  ~Type+Site+Date, data=Rubritalea, na.action=na.omit, link = "logit")
AIC(Mod3)
summary(emmeans(Mod3, pairwise ~ Type, adjust="holm"))
summary(emmeans(Mod3, pairwise ~ Site, adjust="holm"))
summary(emmeans(Mod3, pairwise ~ Date, adjust="holm"))

Blastopirellula<-test[test$Genus=="Blastopirellula",]
Blastopirellula$Type <- factor(Blastopirellula$Type)
Blastopirellula$Site <- factor(Blastopirellula$Site)
Blastopirellula$Date <- factor(Blastopirellula$Date)
Mod4 <- betareg(avg_abundance  ~Type+Site+Date, data=Blastopirellula, na.action=na.omit, link = "logit")
AIC(Mod4)
summary(emmeans(Mod4, pairwise ~ Type, adjust="holm"))
summary(emmeans(Mod4, pairwise ~ Site, adjust="holm"))
summary(emmeans(Mod4, pairwise ~ Date, adjust="holm"))

Arenicella<-test[test$Genus=="Arenicella",]
Arenicella$Type <- factor(Arenicella$Type)
Arenicella$Site <- factor(Arenicella$Site)
Arenicella$Date <- factor(Arenicella$Date)
Mod5<-betareg(avg_abundance ~Type+Site+Date, data=Arenicella,link = "logit",na.action = "na.omit")
AIC(Mod5)
summary(emmeans(Mod5, pairwise ~ Site, adjust="holm"))
summary(emmeans(Mod5, pairwise ~ Type, adjust="holm"))
summary(emmeans(Mod5, pairwise ~ Date, adjust="holm"))

Aliivibrio<-test[test$Genus=="Aliivibrio",]
Aliivibrio$Type <- factor(Aliivibrio$Type)
Aliivibrio$Site <- factor(Aliivibrio$Site)
Aliivibrio$Date <- factor(Aliivibrio$Date)
Mod6<-betareg(avg_abundance ~Type+Site+Date, data=Aliivibrio,link = "logit")
AIC(Mod6)
summary(lsmeans(Mod6, pairwise ~ Site, adjust="holm"))
summary(lsmeans(Mod6, pairwise ~ Type, adjust="holm"))
summary(lsmeans(Mod6, pairwise ~ Date, adjust="holm"))

Cocleimonas<-test[test$Genus=="Cocleimonas",]
Cocleimonas$Type <- factor(Cocleimonas$Type)
Cocleimonas$Site <- factor(Cocleimonas$Site)
Cocleimonas$Date <- factor(Cocleimonas$Date)
Mod7 <- betareg(avg_abundance  ~Type+Site+Date, data=Cocleimonas, na.action=na.omit, link = "logit")
AIC(Mod7)
summary(emmeans(Mod7, pairwise ~ Type, adjust="holm"))
summary(emmeans(Mod7, pairwise ~ Site, adjust="holm"))
summary(emmeans(Mod7, pairwise ~ Date, adjust="holm"))

Thalassotalea<-test[test$Genus=="Thalassotalea",]
Thalassotalea$Type <- factor(Thalassotalea$Type)
Thalassotalea$Site <- factor(Thalassotalea$Site)
Thalassotalea$Date <- factor(Thalassotalea$Date)
Mod8 <- betareg(avg_abundance  ~Type+Site+Date, data=Thalassotalea, na.action=na.omit, link = "logit")
AIC(Mod8)
summary(emmeans(Mod8, pairwise ~ Type, adjust="holm"))
summary(emmeans(Mod8, pairwise ~ Site, adjust="holm"))
summary(emmeans(Mod8, pairwise ~ Date, adjust="holm"))

Wenyingzhuangia<-test[test$Genus=="Wenyingzhuangia",]
Wenyingzhuangia$Type <- factor(Wenyingzhuangia$Type)
Wenyingzhuangia$Site <- factor(Wenyingzhuangia$Site)
Wenyingzhuangia$Date <- factor(Wenyingzhuangia$Date)
Mod9 <- betareg(avg_abundance  ~Type+Site+Date, data=Wenyingzhuangia, na.action=na.omit, link = "logit")
AIC(Mod9)
summary(emmeans(Mod9, pairwise ~ Type, adjust="holm"))
summary(emmeans(Mod9, pairwise ~ Site, adjust="holm"))
summary(emmeans(Mod9, pairwise ~ Date, adjust="holm"))

Granulosicoccus<-test[test$Genus=="Granulosicoccus",]
Granulosicoccus$Type <- factor(Granulosicoccus$Type)
Granulosicoccus$Site <- factor(Granulosicoccus$Site)
Granulosicoccus$Date <- factor(Granulosicoccus$Date)
Mod10 <- betareg(avg_abundance  ~Type+Site+Date, data=Granulosicoccus, na.action=na.omit, link = "logit")
AIC(Mod10)
summary(emmeans(Mod10, pairwise ~ Type, adjust="holm"))
summary(emmeans(Mod10, pairwise ~ Site, adjust="holm"))
summary(emmeans(Mod10, pairwise ~ Date, adjust="holm"))

Psychrobium<-test[test$Genus=="Psychrobium",]
Psychrobium$Type <- factor(Psychrobium$Type)
Psychrobium$Site <- factor(Psychrobium$Site)
Psychrobium$Date <- factor(Psychrobium$Date)
Mod11 <- betareg(avg_abundance  ~Type+Site+Date, data=Psychrobium, na.action=na.omit, link = "logit")
AIC(Mod11)
summary(emmeans(Mod11, pairwise ~ Type, adjust="holm"))
summary(emmeans(Mod11, pairwise ~ Site, adjust="holm"))
summary(emmeans(Mod11, pairwise ~ Date, adjust="holm"))

Leucothrix<-test[test$Genus=="Leucothrix",]
Leucothrix$Type <- factor(Leucothrix$Type)
Leucothrix$Site <- factor(Leucothrix$Site)
Leucothrix$Date <- factor(Leucothrix$Date)
Mod12 <- betareg(avg_abundance  ~Type+Site+Date, data=Leucothrix, na.action=na.omit, link = "logit")
AIC(Mod12)
summary(emmeans(Mod12, pairwise ~ Type, adjust="holm"))
summary(emmeans(Mod12, pairwise ~ Site, adjust="holm"))
summary(emmeans(Mod12, pairwise ~ Date, adjust="holm"))

Algitalea<-test[test$Genus=="Algitalea",]
Algitalea$Type <- factor(Algitalea$Type)
Algitalea$Site <- factor(Algitalea$Site)
Algitalea$Date <- factor(Algitalea$Date)
Mod13 <- betareg(avg_abundance  ~Type+Site+Date, data=Algitalea, na.action=na.omit, link = "logit")
AIC(Mod13)
summary(emmeans(Mod13, pairwise ~ Type, adjust="holm"))
summary(emmeans(Mod13, pairwise ~ Site, adjust="holm"))
summary(emmeans(Mod13, pairwise ~ Date, adjust="holm"))

#HeatMap GLMS-----
library(emmeans)
library(dplyr)

result_lists <- list()
generos <- c("Psychromonas", "Persicirhabdus", "Rubritalea","Blastopirellula","Arenicella",
            "Aliivibrio","Cocleimonas","Thalassotalea","Wenyingzhuangia","Granulosicoccus",
            "Psychrobium","Leucothrix","Algitalea")
for (i in 1:13) {
  modelo <- get(paste0("Mod", i)) 
  genero <- generos[i]
  
tukey_site <- summary(emmeans(modelo, pairwise ~ Site, adjust = "holm")$contrasts)
tukey_site <- tukey_site %>%
    mutate(Genero = genero,
           Variable = "Site")

tukey_date <- summary(emmeans(modelo, pairwise ~ Date, adjust = "holm")$contrasts)
tukey_date <- tukey_date %>%
    mutate(Genero = genero,
           Variable = "Date")
  
tukey_type <- summary(emmeans(modelo, pairwise ~ Type, adjust = "holm")$contrasts)
tukey_type <- tukey_type %>%
    mutate(Genero = genero,
           Variable = "Type")

result_lists [[i]] <- bind_rows(tukey_site, tukey_date, tukey_type)
}

glm_results <- bind_rows(result_lists)

glm_results <- glm_results %>%
  select(Genero, Variable, contrast, estimate, SE, p.value)

glm_results <- glm_results %>%
  mutate(Significativo = ifelse(p.value < 0.05, "Yes", "No"))

ggplot(glm_results, aes(x = contrast, y = Genero, fill = Significativo)) +
  geom_tile(color = "white") +
  geom_text(aes(
    label = ifelse(p.value < 0, "0", round(p.value, 3))
  ), size = 4) +  
  facet_wrap(~ Variable, scales = "free_x") +
  scale_fill_manual(values = c("Yes" = "#d73027", "No" = "#f0f0f0")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
    axis.text.y = element_text(size = 14),                         
    strip.text = element_text(size = 16)           
  ) +
  labs(
    title = "",
    x = "",
    y = "",
    fill = "Significance"
  )

orden_tax <- data.frame(
  Genero = c("Psychromonas", "Persicirhabdus", "Rubritalea","Blastopirellula","Arenicella","Aliivibrio","Cocleimonas",
             "Thalassotalea","Wenyingzhuangia","Granulosicoccus","Psychrobium","Leucothrix","Algitalea" ),
  
  Orden = c("Enterobacterales", "Verrucomicrobiales", "Verrucomicrobiales","Pirellulales","Arenicellales","Enterobacterales","Thiotrichales",
            "Enterobacterales","Flavobacteriales","Granulosicoccales","Enterobacterales","Thiotrichales","Flavobacteriales"))
             


glm_results <- left_join(glm_results, orden_tax, by = "Genero")

glm_results <- glm_results %>%
  mutate(Genero_Orden = paste(Orden, Genero, sep = " | "))


glm_results <- glm_results %>%
  arrange(Orden, Genero) %>%
  mutate(
    Genero_Orden = paste(Orden, Genero, sep = " | "),
    Genero_Orden_label = paste0(Orden, " ~ '|' ~ italic('", Genero, "')"),
    Genero_Orden_label = factor(Genero_Orden_label, levels = unique(Genero_Orden_label))
  )
ggplot(glm_results, aes(x = contrast, y = Genero_Orden_label, fill = Significativo)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(p.value < 0, "0", round(p.value, 3))), size = 4) +
  facet_wrap(~ Variable, scales = "free_x") +
  scale_fill_manual(values = c("Yes" = "#d73027", "No" = "#f0f0f0")) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    strip.text = element_text(size = 16)
  ) +
  labs(
    title = "",
    x = "",
    y = "Order | Genus",
    fill = "Significance"
  ) +
  scale_y_discrete(labels = function(x) parse(text = x))

#############-------------CORE-----------------####################

library(dplyr)
library(vegan)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(ggplot2)
#Subset data
Bac_16S_filt <- subset_samples(Bac_16S, Type != "Sea_Water")

#CORE Type----
#Subset data
apicales <- subset_samples(Bac_16S_filt, Type == "xx")#<- Sample type name
sporophylls <- subset_samples(Bac_16S_filt, Type == "xx")
#Make it relative
apicales_rel <- transform_sample_counts(apicales, function(x) x / sum(x))
sporophylls_rel <- transform_sample_counts(sporophylls, function(x) x / sum(x))

otu_table(apicales_rel) <- otu_table(apicales_rel, taxa_are_rows = TRUE)
otu_table(sporophylls_rel) <- otu_table(sporophylls_rel, taxa_are_rows = TRUE)
#Pick the core 
core_apicales <- core(apicales_rel, detection = 0.001, prevalence = 0.90)
core_sporophylls <- core(sporophylls_rel, detection = 0.001, prevalence = 0.90)

asv_apicales <- taxa_names(core_apicales)
asv_sporophylls <- taxa_names(core_sporophylls)
#Assign to list to manipulate data
asv_apicales_list <- c(xx, xx,.... ) #<-ID ASVs
#Save list in .csv
write.csv(asv_apicales_lista , "core_apical.csv", row.names = TRUE)
#Assign to list to manipulate data
asv_sporophylls_list <- c(xx, xx,....)#<-ID ASVs
#Save list in .csv
write.csv(asv_sporophylls_list, "core_sporophylls.csv", row.names = TRUE)
#Read the .csv modified with the names of bacterial genera
core_apical <- read.csv("core_apical.csv")
core_sporophylls <- read.csv("core_sporophylls.csv")
#Subset data
asv_apical <- unique(core_apical$ASV)
asv_sporophylls <- unique(core_sporophylls$ASV)

asv_list <- list(
  Sporophylls = asv_sporophylls,
  Apical = asv_apical)
#Subset data
core_apical$Origen <- "xx"#<- Sample type name
core_sporophylls$Origen <- "xx"#<- Sample type name
core_total <- bind_rows(core_apical, core_sporophylls)

asv_comunes <- intersect(asv_apical, asv_sporophylls)
asv_solo_apical <- setdiff(asv_apical, asv_sporophylls)
asv_solo_sporophylls <- setdiff(asv_sporophylls, asv_apical)

# Comun ASVs
gen_inter <- core_total %>%
  filter(ASV %in% asv_comunes) %>%
  distinct(ASV, .keep_all = TRUE) %>%
  group_by(Genus) %>%
  summarise(n = n()) %>%
  mutate(texto = paste0(Genus, " (", n,")"))

# only ASVs apical
gen_apical <- core_apical %>%
  filter(ASV %in% asv_solo_apical) %>%
  group_by(Genus) %>%
  summarise(n = n()) %>%
  mutate(texto = paste0(Genus, " (", n,")"))

# only ASVs sporophylls
gen_sporophylls<- core_sporophylls %>%
  filter(ASV %in% asv_solo_sporophylls) %>%
  group_by(Genus) %>%
  summarise(n = n()) %>%
  mutate(texto = paste0(Genus, " (", n,")"))


texto_apical <- paste(gen_apical$texto, collapse = "\n")
texto_esporofila <- paste(gen_sporophylls$texto, collapse = "\n")
texto_inter <- paste(gen_inter$texto, collapse = "\n")

#Plote the Venn diagram
venn <- ggVennDiagram(asv_lista,label = "none", label_alpha = 0, show_label = FALSE) +
  scale_fill_gradient(low = "#FFFF00", high = "#C1FFC1") +
  theme_void() +
  theme(legend.position = "none") +
  ggtitle("")

venn_type_plot <-venn +
  annotate("text", x = -0.2, y = 8.5, label = "Apical", size = 8, fontface = "plain") +
  annotate("text", x = 0.15, y = -4.5, label = "Sporophylls", size = 8, fontface = "plain") +
  annotate("text", x = -1.4, y = 6.1, label = texto_apical, size = 8, hjust = 0, fontface = "italic") +
  annotate("text", x = -1.2, y = -1.7, label = texto_sporophylls, size = 8, hjust = 0, fontface = "italic") +
  annotate("text", x = -1.8, y = 2, label = texto_inter, size = 7, hjust = 0, fontface = "italic")


# CORE Date ----
#Subset data
August<- subset_samples(Bac_16S_filt, Date == "August")
March <- subset_samples(Bac_16S_filt, Date == "March")
#Make it relative
August_rel <- transform_sample_counts(August, function(x) x / sum(x))
March_rel <- transform_sample_counts(March, function(x) x / sum(x))

otu_table(August_rel) <- otu_table(August_rel, taxa_are_rows = TRUE)
otu_table(March_rel) <- otu_table(March_rel, taxa_are_rows = TRUE)
#Pick the core 
core_August <- core(August_rel, detection = 0.001, prevalence = 0.90)
core_March <- core(March_rel, detection = 0.001, prevalence = 0.90)

asv_August <- taxa_names(core_August)
asv_March <- taxa_names(core_March)

#Assign to list to manipulate data
asv_August <- c("xx,xx, ....." )
#Save list in .csv
write.csv(asv_August, "core_August.csv", row.names = TRUE)
asv_March <- c("xx, xx,....")
write.csv(asv_March, "core_March.csv", row.names = TRUE)

#Read the .csv modified with the names of bacterial genera
core_August <- read.csv("core_August.csv")
core_March <- read.csv("core_March.csv")

asv_August <- unique(core_Augustt$ASV)
asv_March <- unique(core_March$ASV)

asv_list <- list(
  August = asv_August,
  March = asv_march)

core_August$Origen <- "August"
core_March$Origen <- "March"
core_total <- bind_rows(core_August, core_March)

# Inter
asv_comun_date <- intersect(asv_August, asv_March)
asv_solo_August <- setdiff(asv_August, asv_March)
asv_solo_March <- setdiff(asv_March, asv_August)

# Comun ASVs
gen_inter_date <- core_total %>%
  filter(ASV %in% asv_comunes_date) %>%
  distinct(ASV, .keep_all = TRUE) %>%
  group_by(Genus) %>%
  summarise(n = n()) %>%
  mutate(texto = paste0(Genus, " (", n, ")"))

# only ASVs August
gen_August <- core_August %>%
  filter(ASV %in% asv_solo_August) %>%
  group_by(Genus) %>%
  summarise(n = n()) %>%
  mutate(texto = paste0(Genus, " (", n, ")"))

# only ASVs March
gen_March <- core_March%>%
  filter(ASV %in% asv_solo_March) %>%
  group_by(Genus) %>%
  summarise(n = n()) %>%
  mutate(texto = paste0(Genus, " (", n, ")"))


texto_August <- paste(gen_August$texto, collapse = "\n")
texto_March <- paste(gen_March$texto, collapse = "\n")
texto_inter_date <- paste(gen_inter_date$texto, collapse = "\n")

names(asv_lista) <- c("", "") 
names(asv_lista) <- c("", "") 

#PLot Venn diagram
venn <- ggVennDiagram(asv_lista, label = "none", label_alpha = 0, show_label = FALSE) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_void() +
  theme(legend.position = "none") +
  ggtitle("")

venn_date_plot <-venn +
  annotate("text", x = -0.2, y = 8.5, label = "August", size = 8, fontface = "plain") +
  annotate("text", x = 0.15, y = -4.5, label = "March", size = 8, fontface = "plain") +
  annotate("text", x = -1.2, y = 6, label = texto_august, size = 6, hjust = 0, fontface = "italic") +
  annotate("text", x = -1.5, y = 2.2, label = texto_inter_date, size = 6, hjust = 0, fontface = "italic") +
  annotate("text", x = -1.1, y = -1.7, label = texto_march, size = 6, hjust = 0, fontface = "italic")

# CORE Site ----
#Subset data
BB <- subset_samples(Bac_16S_filt, Site == "Bahia_Buzo")
SG <- subset_samples(Bac_16S_filt, Site == "San_Gregorio")
BQ <- subset_samples(Bac_16S_filt, Site == "Buque_Quemado")


#Make it relative
BB_rel <- transform_sample_counts(BB, function(x) x / sum(x))
SG_rel <- transform_sample_counts(SG, function(x) x / sum(x))
BQ_rel <- transform_sample_counts(BQ, function(x) x / sum(x))

otu_table(BB_rel) <- otu_table(BB_rel, taxa_are_rows = TRUE)
otu_table(SG_rel) <- otu_table(SG_rel, taxa_are_rows = TRUE)
otu_table(BQ_rel) <- otu_table(BQ_rel, taxa_are_rows = TRUE)

#Pick the core 
core_BB <- core(BB_rel, detection = 0.001, prevalence = 0.90)
core_SG <- core(SG_rel, detection = 0.001, prevalence = 0.90)
core_BQ <- core(BQ_rel, detection = 0.001, prevalence = 0.90)

asv_BB <- taxa_names(core_BB)
asv_SG <- taxa_names(core_SG)
asv_BQ <- taxa_names(core_BQ)

#Assign to list to manipulate data
asv_BB <- c("ASV_10","ASV_14","ASV_2","ASV_4","ASV_5","ASV_9" )
#Save list in .csv
#write.csv(asv_BB, "core_BB.csv", row.names = TRUE)
asv_SG <- c("ASV_2","ASV_3","ASV_4","ASV_9")
#write.csv(asv_SG, "core_SG.csv", row.names = TRUE)
asv_BQ <- c("ASV_10","ASV_13","ASV_17","ASV_18","ASV_19","ASV_2","ASV_21","ASV_22","ASV_27",
            "ASV_3","ASV_30","ASV_36","ASV_4","ASV_5","ASV_50","ASV_6","ASV_9","ASV_103","ASV_94")
#write.csv(asv_BQ, "core_BQ.csv", row.names = TRUE)

#Read the .csv modified with the names of bacterial genera
core_BB_df <- read.csv("core_BB.csv")
core_SG_df <- read.csv("core_SG.csv")
core_BQ_df <- read.csv("core_BQ.csv")


library(dplyr)

core_BB_df$ASV <- as.character(core_BB_df$ASV)
core_SG_df$ASV <- as.character(core_SG_df$ASV)
core_BQ_df$ASV <- as.character(core_BQ_df$ASV)


asv_list <- list(
  Bahia_Buzo = asv_BB,
  San_Gregorio = asv_SG,
  Buque_Quemado =  asv_BQ
  )

core_BB_df$Origen <- "Bahia_Buzo"
core_SG_df$Origen <- "San_Gregorio"
core_BQ_df$Origen <- "Buque_Quemado"

core_total <- bind_rows(core_BB_df, core_SG_df, core_BQ_df)


# Inter
asv_BB <- unique(core_BB_df$ASV)
asv_SG <- unique(core_SG_df$ASV)
asv_BQ <- unique(core_BQ_df$ASV)

asv_comun_site <- Reduce(intersect, list(asv_BB, asv_SG, asv_BQ))
asv_solo_BB <- setdiff(asv_BB, union(asv_SG, asv_BQ))
asv_solo_SG <- setdiff(asv_SG, union(asv_BB, asv_BQ))
asv_solo_BQ <- setdiff(asv_BQ, union(asv_BB, asv_SG))

asv_BB_SG <- setdiff(intersect(asv_BB, asv_SG), asv_BQ)
asv_BB_BQ <- setdiff(intersect(asv_BB, asv_BQ), asv_SG)
asv_SG_BQ <- setdiff(intersect(asv_SG, asv_BQ), asv_BB)

gen_por_region <- function(asvs, df) {
  df %>%
    filter(ASV %in% asvs) %>%
    distinct(ASV, .keep_all = TRUE) %>%
    group_by(Genus) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(texto = paste0(Genus, " (", n, ")"))
}
gen_BB      <- gen_por_region(asv_solo_BB, core_BB_df)
gen_SG      <- gen_por_region(asv_solo_SG, core_SG_df)
gen_BQ      <- gen_por_region(asv_solo_BQ, core_BQ_df)

gen_BB_SG   <- gen_por_region(asv_BB_SG, core_total)
gen_BB_BQ   <- gen_por_region(asv_BB_BQ, core_total)
gen_SG_BQ   <- gen_por_region(asv_SG_BQ, core_total)

gen_all     <- gen_por_region(asv_comun_site, core_total)
split_two_columns <- function(x) {
  n <- length(x)
  mid <- ceiling(n / 2)
  col1 <- x[1:mid]
  col2 <- x[(mid + 1):n]
  
  # rellenar con "" si están desbalanceadas
  if (length(col1) > length(col2)) {
    col2 <- c(col2, rep("", length(col1) - length(col2)))
  }
  
  paste(col1, col2, sep = "    ", collapse = "\n")
}

texto_BB    <- paste(gen_BB$texto, collapse = "\n")
texto_SG    <- paste(gen_SG$texto, collapse = "\n")
texto_BQ <- split_two_columns(gen_BQ$texto)

texto_BB_SG <- paste(gen_BB_SG$texto, collapse = "\n")
texto_BB_BQ <- paste(gen_BB_BQ$texto, collapse = "\n")
texto_SG_BQ <- paste(gen_SG_BQ$texto, collapse = "\n")

texto_all   <- paste(gen_all$texto, collapse = "\n")

library(ggVennDiagram)


venn <- ggVennDiagram(
  asv_list,
  label = "none",        # no números
  set_label = "none"     # <- ESTO apaga los nombres de los sitios
) +
  scale_fill_gradient(low = "#98F5FF", high = "#53868B") +
  theme_void() +
  theme(legend.position = "none") +
  ggtitle("")

venn_site_plot <- venn +
  # títulos
  annotate("text", x = -1, y = 4.5, label = "Bahía Buzo", size = 7) +
  annotate("text", x =  4.5, y = 4.5, label = "San Gregorio", size = 7) +
  annotate("text", x =  2.1, y = -8.5, label = "Buque Quemado", size = 7) +
  
  # regiones exclusivas
  annotate("text", x = -3.5, y =  0.5, label = texto_BB, size = 4.5, hjust = 0, fontface = "italic") +
  annotate("text", x =  6.0, y =  1.5, label = texto_SG, size = 4.5, hjust = 0, fontface = "italic") +
  annotate("text", x = -0.7, y = -5.5, label = texto_BQ, size = 4.5, hjust = 0, fontface = "italic") +
  
  # intersecciones dobles
  annotate("text", x =  1.0, y =  2.5, label = texto_BB_SG, size = 4.5, hjust = 0, fontface = "italic") +
  annotate("text", x = -1.5, y = -3, label = texto_BB_BQ, size = 4.5, hjust = 0, fontface = "italic") +
  annotate("text", x =  3, y = -3, label = texto_SG_BQ, size = 4.5, hjust = 0, fontface = "italic") +
  
  # intersección triple
  annotate("text", x =  0.35, y =  -1, label = texto_all, size = 5, hjust = 0, fontface = "italic")

venn_site_plot


venn_type_plot <- venn_type_plot +
  annotate("text", x = -4, y = 7, label = "A", size = 10, fontface = "bold")

venn_date_plot <- venn_date_plot +
  annotate("text", x = -4, y = 7, label = "C", size = 10, fontface = "bold")

venn_site_plot <- venn_site_plot +
  annotate("text", x = -4, y = 7, label = "B", size = 10, fontface = "bold")

combined_plot <- venn_type_plot + venn_site_plot+ venn_date_plot  + 
  plot_layout(ncol = 3)

combined_plot 


# CORE global ----
Bac_rel <- transform_sample_counts(Bac_16S_filt, function(x) x / sum(x))
otu_table(Bac_rel) <- otu_table(Bac_rel, taxa_are_rows = TRUE)

core_global <- core(Bac_rel, detection = 0.001, prevalence = 0.90)
asv_global <- taxa_names(core_global)

length(asv_global)

core_global_df <- as.data.frame(tax_table(core_global))
core_global_df$ASV <- rownames(core_global_df)

core_rel <- prune_taxa(asv_global, Bac_rel)

library(phyloseq)
library(dplyr)

core_long <- psmelt(core_rel)


library(dplyr)

core_long <- core_long %>%
  mutate(
    OTU = as.factor(OTU),
    Type = as.factor(Type),
    Date = as.factor(Date),
    Site = as.factor(Site),
    Abundance = as.numeric(Abundance)
  )
core_type_sum <- core_long %>%
  group_by(Type, OTU) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop")

core_date_sum <- core_long %>%
  group_by(Date, OTU) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop")

core_site_sum <- core_long %>%
  group_by(Site, OTU) %>%
  summarise(mean_abund = mean(Abundance), .groups = "drop")

Ty <- ggplot(core_type_sum, aes(x = Type, y = mean_abund, fill = OTU)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "Mean relative abundance", x = "Blade type", fill = "Core ASVs")
Da <- ggplot(core_date_sum, aes(x = Date, y = mean_abund, fill = OTU)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "Mean relative abundance", x = "Sampling date", fill = "Core ASVs")
Si <- ggplot(core_site_sum, aes(x = Site, y = mean_abund, fill = OTU)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(y = "Mean relative abundance", x = "Site", fill = "Core ASVs")
library(ggplot2)
library(patchwork)

Cor <- Ty + Si + Da
Cor



