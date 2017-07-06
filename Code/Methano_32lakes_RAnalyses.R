###### Environmental drivers of methanogen community assembly in lake sediments
###### Brittni L. Bertolet, William E. West, David W. Armitage, and Stuart E. Jones
###### Statistical Analyses

rm(list = ls())

# Load packages 
library(ggplot2)
library(cowplot)
library(phyloseq)
library(vegan)

# Set working directory
setwd("/Users/brittnibertolet/Dropbox/methanogenData/Bertolet/Final Draft/Data")

####################################################################
# Read in data as phyloseq object
otufile = "pelagic_otus.txt"                  # OTU table
mapfile = "envData_32lakes_map.txt"           # Env data
trefile = "methanogen_phylo.tre"              # Phylogenetic tree
refseqfile= "methanogen_rep_set.fna"          # Representative sequences
methanos = import_qiime(otufile,mapfile,trefile,refseqfile) # Merge into phyloseq object

# Root phylogeny for UniFrac
phy_tree(methanos) <- root(phy_tree(methanos), sample(taxa_names(methanos), 1), resolve.root = TRUE)

# Recover OTU table if it is needed elsewhere
OTU <- as(otu_table(methanos), "matrix")


####################################################################
# Alpha diversity calculation
observed_richness = estimate_richness(methanos, measures ="Shannon")

diversity = observed_richness$Shannon
evenness = observed_richness$Shannon/max(observed_richness$Shannon)

divdat = data.frame(sample_data(methanos), diversity, evenness)
divdat$pH2 = divdat$pH^2

# Do linear regression to test for unimodal effect of pH on OTU alpha diversity
mod1 = glm(diversity ~ pH, data = divdat, family = gaussian)
mod2 = glm(diversity ~ pH + pH2, data = divdat, family = gaussian)
mod0 = glm(diversity ~ 1, data = divdat, family = gaussian)

summary(mod1)
summary(mod2)
summary(mod0)

# Plot results of regression
shannon = ggplot(divdat, aes(x = pH, y = diversity)) + geom_point(pch = 1, size = 4) +
  stat_smooth(method = "lm", formula = y~x+I(x^2), se = F, color = "black") +
  scale_y_continuous("Shannon diversity (H')") + scale_x_continuous("Lake pH")

# Save plot to pdf
save_plot("Figure1.pdf", shannon, base_width = 3, base_height = 3)

####################################################################
### Standardizations
# Relative abundance transformation 
methanos -> methanos_rel
otu_table(methanos_rel) <- otu_table(decostand(otu_table(methanos), MARGIN = 2, method = "total"), taxa_are_rows=TRUE)
# Rarification transformation
set.seed(666)
methanos_rare = rarefy_even_depth(methanos)

####################################################################
# Create dissimilarity matrix for PERMANOVA
dist_rel = distance(methanos_rel, "bray") # with relative abundance transformation
dist_rare = distance(methanos_rare, "bray") # with rarifaction

# Perform PERMANOVAs with both transformations
### pH
fit_rel <- adonis(dist_rel ~ methanos_rel@sam_data$pH)
fit_rare <- adonis(dist_rare ~ methanos_rare@sam_data$pH)
### DOC
fit_rel <- adonis(dist_rel ~ methanos_rel@sam_data$DOC)
fit_rare <- adonis(dist_rare ~ methanos_rare@sam_data$DOC)
### chla
fit_rel <- adonis(dist_rel ~ methanos_rel@sam_data$chlA)
fit_rare <- adonis(dist_rare ~ methanos_rare@sam_data$chlA)
### TP
fit_rel <- adonis(dist_rel ~ methanos_rel@sam_data$TP)
fit_rare <- adonis(dist_rare ~ methanos_rare@sam_data$TP)
### TN 
fit_rel <- adonis(dist_rel ~ methanos_rel@sam_data$TN)
fit_rare <- adonis(dist_rare ~ methanos_rare@sam_data$TN)

####################################################################
# Perform PCoAs with different transformations and ordination metrics 
ordu = ordinate(methanos_rel, "PCoA", "bray")
ordu_rare = ordinate(methanos_rare, "PCoA", "bray")
ordu_unifrac = ordinate(methanos_rel, "PCoA", "unifrac")
ordu_unifracW = ordinate(methanos_rel, "PCoA", "wunifrac")

# Plot Relative Abundance & Bray-Curtis PCoA
p1 = plot_ordination(methanos_rel, ordu, color="pH") + geom_point(size = 4) + 
  scale_colour_gradient(low = "black", high = "grey85") +
  scale_x_continuous(breaks = c(-0.2,0,0.2,0.4)) +  scale_y_continuous(breaks = c(-0.2,0,0.2,0.4))

# Pull PCoA axes 1
divdat$bray_rel1 = ordu$vectors[,1]
divdat$bray_rare1 = ordu_rare$vectors[,1]
divdat$uni1 = ordu_unifrac$vectors[,1]
divdat$wuni1 = ordu_unifracW$vectors[,1]
divdat$wuni2 = ordu_unifracW$vectors[,2]

# Do regressions on pH vs PCO1 or 2
summary(lm(divdat$bray_rel1~divdat$pH))
summary(lm(divdat$bray_rare1~divdat$pH))
summary(lm(divdat$uni1~divdat$pH))
summary(lm(divdat$wuni1~divdat$pH))
summary(lm(divdat$wuni2~divdat$pH))

# Make R-squared labels for figs
lb1 <- paste("R^2 == ", 0.50)

# Plot pH ~ PCO1 graph
q1 = ggplot(divdat, aes(y=bray_rel1, x = pH, color = pH)) + geom_point(size = 4) + stat_smooth(method = "lm", se = F, color = "black") +
  scale_colour_gradient(low = "black", high = "grey85") +   scale_y_continuous("PCoA 1") + xlab("Lake pH") +
  annotate("text", x=-Inf, y=-Inf, hjust=-0.25, vjust=-1, label=lb1, parse=TRUE)

# Arrange into panel
prow <- plot_grid( p1 + theme(legend.position="none"),
                   q1 + theme(legend.position="none"),
                   align = 'vh',
                   labels = c("A", "B"),
                   hjust = -1,
                   nrow = 1
)
legend <- get_legend(p1)
p <- plot_grid(prow, legend, nrow = 1, rel_widths = c(3, .3))

# Save plot
save_plot("Figure2.pdf",p, base_width = 6.5, base_height = 3)

####################################################################
# Perform order level analyses 
# Aggregate to order level
OTU_rel <- as.data.frame(methanos_rel@otu_table) # Get OTU table
OTU_rel <- OTU_rel[order(row.names(OTU_rel)),] # Order by row name
tax <- as.data.frame(methanos_rel@tax_table) # Get taxonomy table
tax <- tax[order(row.names(tax)),] # Order by row name
OTU_rel$Order <- tax$Order # Add order taxonomy to OTU table

ord_table <- matrix(NA, ncol = 5, nrow = 32) # Create output table
row.names(ord_table) <- colnames(OTU_rel[,1:32])
orders <- as.vector(unique(tax$Order))
colnames(ord_table) <- orders

for(i in 1:32){
  for(j in 1:5){
    ord_table[i,j] <- sum(OTU_rel[OTU_rel$Order == orders[j],i])
  }
} # for loop for summing relative abundance of each order for each lake 

# Create dissimilarity matrix for Order level
ord_dist <- vegdist(ord_table, "bray")
# Compare Order level dissimilarity matrix to OTU level 
fit <- mantel(dist_rel, ord_dist)
# Put relative abundances of Orders into divdat table
divdat <- cbind(divdat, ord_table)

# Do correlations on pH vs Order
summary(lm(divdat$Methanobacteriales~divdat$pH))
summary(lm(divdat$Methanosarcinales~divdat$pH))
summary(lm(divdat$Methanomicrobiales~divdat$pH))
summary(lm(divdat$E2~divdat$pH))
summary(lm(divdat$Methanocellales~divdat$pH))

# Create order composition graph 
library(reshape2)
ord_table2 <- melt(ord_table)
palette1 <- c('black','darkgrey','grey90','rosybrown2','red3')
ord_table2$Var1 <- as.character(ord_table2$Var1)
ord_table2$Var1 <- factor(ord_table2$Var1, 
                               levels = c("CB", "TR", "NG", "BE", "HB", "BO", "FO", "CR", "WL", "TU", "BA", "PA", "SN", "JU", "PE", "BC",
                               "DI", "DO", "MO", "CH", "CO", "GA", "MA", "RO", "LW", "PL", "ST", "UF", "NC", "BR", "GV", "HU"))

p2 <- ggplot(ord_table2, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity") + ylab("Relative Abundance") + xlab("Lake ID") +
  scale_fill_manual(values=palette1) + 
  theme(text = element_text(size=12), axis.text.x = element_text(size = 8), axis.text.y = element_text(angle=90, hjust=1))

# Create orderVpH graphs
quartz()
p3 <- ggplot(divdat, aes(x = pH, y = Methanocellales)) + geom_point() +
  stat_smooth(method = "lm", se = F, color = "black") + theme(text = element_text(size=10),axis.text.y = element_text(angle=90, hjust=1)) +
  theme(axis.text.y = element_text(size=10),axis.text.x = element_text(size=10))
p4 <- ggplot(divdat, aes(x = pH, y = E2)) + geom_point() + theme(text = element_text(size=10),axis.text.y = element_text(angle=90, hjust=1)) +
  theme(axis.text.y = element_text(size=10),axis.text.x = element_text(size=10))
p5 <- ggplot(divdat, aes(x = pH, y = Methanomicrobiales)) + geom_point() +
  stat_smooth(method = "lm", se = F, color = "black") + theme(text = element_text(size=10),axis.text.y = element_text(angle=90, hjust=1)) +
  theme(axis.text.y = element_text(size=10),axis.text.x = element_text(size=10))
p6 <- ggplot(divdat, aes(x = pH, y = Methanosarcinales)) + geom_point() + theme(text = element_text(size=10),axis.text.y = element_text(angle=90, hjust=1)) +
  theme(axis.text.y = element_text(size=10),axis.text.x = element_text(size=10))
p7 <- ggplot(divdat, aes(x = pH, y = Methanobacteriales)) + geom_point() + theme(text = element_text(size=10),axis.text.y = element_text(angle=90, hjust=1)) +
  theme(axis.text.y = element_text(size=10),axis.text.x = element_text(size=10))


p8 <- plot_grid(p3 + theme(legend.position="none"),
                   p4 + theme(legend.position="none"),
                   p5 + theme(legend.position="none"),
                   p6 + theme(legend.position="none"),
                   p7 + theme(legend.position="none"),
                   nrow = 1
)
legend <- get_legend(p2)

# Put together with composition graph
p9 <- plot_grid(legend,
                  p2 + theme(text = element_text(size=12), legend.position = "none"),
                  p8 + theme(text = element_text(size=10)),
                  labels = c("A", "B"),
                  nrow = 3, rel_heights = c(1, 1.5, 1)
)
save_plot("Figure3.pdf",prow2, base_width = 7, base_height = 6)

####################################################################
# Perform order level correlations after accounting for pH
summary(lm(divdat$E2~divdat$pH + divdat$Methanobacteriales)) 
summary(lm(divdat$E2~divdat$pH + divdat$Methanosarcinales))
summary(lm(divdat$E2~divdat$pH + divdat$Methanomicrobiales)) # this is the only one that is significant
summary(lm(divdat$E2~divdat$pH + divdat$Methanocellales))

summary(lm(divdat$Methanobacteriales~divdat$pH + divdat$Methanosarcinales))
summary(lm(divdat$Methanobacteriales~divdat$pH + divdat$Methanomicrobiales))
summary(lm(divdat$Methanobacteriales~divdat$pH + divdat$Methanocellales))

summary(lm(divdat$Methanosarcinales~divdat$pH + divdat$Methanomicrobiales))
summary(lm(divdat$Methanosarcinales~divdat$pH + divdat$Methanocellales))

summary(lm(divdat$Methanomicrobiales~divdat$pH + divdat$Methanosarcinales))


# Get the residuals
divdat$r.fit <- resid(lm(divdat$E2~divdat$pH))
summary(lm(divdat$r.fit~divdat$Methanomicrobiales))
R = paste("R^2 == ", 0.42)

q1 = ggplot(divdat, aes(y=r.fit, x = Methanomicrobiales, color = pH)) + geom_point(size = 4) +
  stat_smooth(method = "lm", se = F, color = "black") + scale_colour_gradient(low = "black", high = "grey85") + 
  scale_y_continuous("Residuals") + xlab("Methanomassiliicoccales") +
  annotate("text", x=-Inf, y=-Inf, hjust=-0.25, vjust=-1, label = R, parse=TRUE) 

save_plot("Figure4.pdf",q1, base_width = 5, base_height = 5)

####################################################################
# Perform UNDERC analysis
otufile = "UNDERC_otu.txt"             # OTU table
mapfile = "UNDERC_envData_map.txt"      # Env data
underc = import_qiime(otufile,mapfile) # Merge into phyloseq object

# Relative abundance transformation 
underc -> underc_rel
otu_table(underc_rel) <- otu_table(decostand(otu_table(underc), method = "total"), taxa_are_rows=TRUE)

# Create dissimilarity matrix for PERMANOVA
dist_rel = distance(underc_rel, "bray")

# Perform PCoAs with different transformations and ordination metrics 
ordu = ordinate(underc_rel, "PCoA", "bray")

# Determine correlation with PCoA 2
summary(lm(underc_rel@sam_data$Prod~underc_rel@sam_data$chlA + ordu$vectors[,2]))
fit.r <- resid(lm(underc_rel@sam_data$Prod~underc_rel@sam_data$chlA)) # residuals
cor(fit.r, ordu$vectors[,2]) # correlation, r = 0.51


