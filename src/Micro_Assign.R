
# Soil Microorganisms 2023 - Diversity Analysis ####

library(tidyr) # data management
library(dplyr) # data management
library(iNEXT) # Hill numbers
library(readr) # for write_csv (within extrar_divmetrics function)

# 1. Unassigned ASV exploration ####

# Micro_2023 <- read.csv("data/Micro_2023.csv", fileEncoding = "latin1", na.strings = c("", "NA"))
#   
# Micro_2023 <- Micro_2023 %>% # Transpose data frame
#         pivot_longer(cols = starts_with("M"), # Transposing will be applied only to sample name (M___) columns
#                      names_to = "Sample",      # Transposed columnes are named "Sample"
#                      values_to = "Abundance") # All values from transposed columns go to an "Abundance" column
# 
# Sample_info <- read.csv("data/Sample_info.csv", fileEncoding = "latin1", na.strings = c("", "NA")) # Importing Sample information
# 
# Micro_2023 <- Micro_2023 %>% 
#         left_join(Sample_info, by = "Sample") # Including sample information
#   
# Micro_2023 <- Micro_2023 %>%   
#         filter(!Sampling %in% c(5, 6)) %>%  # Samplings 5 and 6 dropped as they were collected during Fallow Season.
#         filter(Abundance != 0) # Remove all rows with 0 Abundance
# 
# # Checking weight of non assigned Taxonomy:
# filter(Micro_2023, Abundance > 0, grepl("^unclassified", Class)) %>%  nrow() # rows where Class = "unclassified..." = 5661 (1.63%)
# filter(Micro_2023, Abundance > 0, grepl("^uncultured", Class)) %>%  nrow() # rows where Class = "uncultured..." = 0
# filter(Micro_2023, Abundance > 0, grepl("^unclassified", Order)) %>%  nrow() # rows where Order = "unclassified..." = 22439 (6.46%)
# filter(Micro_2023, Abundance > 0, grepl("^uncultured", Order)) %>%  nrow() # rows where Order = "uncultured..." = 885 (0.25%)
# filter(Micro_2023, Abundance > 0, grepl("^unclassified", Family)) %>%  nrow() # rows where Family = "unclassified..." = 58131 (16.75%)
# filter(Micro_2023, Abundance > 0, grepl("^uncultured", Family)) %>%  nrow() # rows where Family = "uncultured..." = 26001 (7.49%) 
# filter(Micro_2023, Abundance > 0, grepl("^unclassified", Genus)) %>%  nrow() # rows where Genus = "unclassified..." = 177569 (51.19%) 
# filter(Micro_2023, Abundance > 0, grepl("^uncultured", Genus)) %>%  nrow() # rows where Genus = "uncultured..." = 33983 (9.79%) 
# 
# # Checking weight of non assigned Family per Treat:
# filter(Micro_2023, Treat == "CON", Abundance > 0, grepl("^unclassified", Family)) %>%  nrow() # rows where Family = "unclassified..." = 18667 (%)
# filter(Micro_2023, Treat == "MSD", Abundance > 0, grepl("^unclassified", Family)) %>%  nrow() # rows where Family = "unclassified..." = 19164 (%)
# filter(Micro_2023, Treat == "AWD", Abundance > 0, grepl("^unclassified", Family)) %>%  nrow() # rows where Family = "unclassified..." = 20300 (%)
# filter(Micro_2023, Treat == "CON", Abundance > 0, grepl("^uncultured", Family)) %>%  nrow() # rows where Family = "unclassified..." = 8528 (%)
# filter(Micro_2023, Treat == "MSD", Abundance > 0, grepl("^uncultured", Family)) %>%  nrow() # rows where Family = "unclassified..." = 8591 (%)
# filter(Micro_2023, Treat == "AWD", Abundance > 0, grepl("^uncultured", Family)) %>%  nrow() # rows where Family = "unclassified..." = 8882 (%)

# 2. Diversity at ASV level ####

## 2.1 Preparing data for Hill numbers calculation ####

Micro_2023_ASV <- read.csv("data/ASV_TABLE_PARCELA.csv", fileEncoding = "latin1", na.strings = c("", "NA"))

Micro_2023_ASV <- Micro_2023_ASV %>% # Transpose data frame
        pivot_longer(cols = starts_with("P"), # Transposing will be applied only to sample name (P___) columns
                     names_to = "Plot",      # Transposed columnes are named "Plot"
                     values_to = "Abundance") %>% # All values from transposed columns go to an "Abundance" colum
        filter(Abundance != 0) # Remove all rows with 0 Abundance

# Assigning Treats to plots:
Plot <- c("P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P12", "P13", "P14", "P15")
Treat <- c("AWD", "MSD", "CON", "MSD", "AWD", "CON", "MSD", "CON", "AWD", "AWD", "MSD", "CON", "MSD", "AWD", "CON")
Plot_Treat <- data.frame(Plot, Treat)

Micro_2023_ASV <- Micro_2023_ASV %>% left_join(Plot_Treat, by = "Plot")

Micro_2023_ASV <- Micro_2023_ASV %>%
        mutate(siteID = paste0(Plot, "_", Treat)) # Create siteID column to use as variable within extrar_divmetrics

## 2.2.  Calculate Hill numbers through extrar_divmetrics function ####

extrar_divmetrics <- function (datafile){ # Function that calculates Hill numbers using the iNEXT package.
  
        list.sites <- datafile %>%
          group_by(siteID) %>%
          summarise(list(Abundance))
        
        list.sites2 <- list.sites[[2]]
        
        names(list.sites2) <- list.sites[[1]]
        
        x <- iNEXT::iNEXT(list.sites2, q = 0, datatype = "abundance", size = NULL)
        #ggiNEXT(x, type=1, se=TRUE, facet.var="none", color.var="none", grey=FALSE)
        
        inextparams <- as.data.frame(x[3]) %>%
          dplyr::select(siteID = AsyEst.Assemblage,
                        divmetric = AsyEst.Diversity,
                        obsmetric = AsyEst.Observed,
                        estmetric = AsyEst.Estimator,
                        se.metric = AsyEst.s.e.)
        
        rich.df <- inextparams %>%
          filter(divmetric == "Species richness") %>%
          dplyr::select(siteID,
                        q0.obs = obsmetric,
                        q0.est = estmetric,
                        q0.se = se.metric)
        
        shannon.df <- inextparams %>%
          filter(divmetric == "Shannon diversity") %>%
          dplyr::select(siteID,
                        q1.obs = obsmetric,
                        q1.est = estmetric,
                        q1.se = se.metric)
        
        simpson.df <- inextparams %>%
          filter(divmetric == "Simpson diversity") %>%
          dplyr::select(siteID,
                        q2.obs = obsmetric,
                        q2.est = estmetric,
                        q2.se = se.metric)
        
        div.params <- rich.df %>%
          left_join(., shannon.df, by = "siteID") %>%
          left_join(., simpson.df, by = "siteID")
        
        dataname <-  deparse(substitute(datafile))
        write_csv(div.params,paste0("data/2023/", dataname,"_inext_params.csv"))
  
}

Micro_2023_ASV_inext_params <- extrar_divmetrics(Micro_2023_ASV)                                     ######### THIS STEP TAKES TIME TO RUN #########

Hills_Micro_2023_ASV <- Micro_2023_ASV %>%
        group_by(Plot, Treat, siteID) %>%
        summarise(siteID = first(siteID), Plot = first(Plot), Treat = first(Treat)) %>%
        left_join(Micro_2023_ASV_inext_params, by = "siteID")

Hills_Micro_2023_ASV$Treat <- factor(Hills_Micro_2023_ASV$Treat, levels = c("CON", "MSD", "AWD"))  # Reorder the Treat variable
Hills_Micro_2023_ASV <- Hills_Micro_2023_ASV[order(Hills_Micro_2023_ASV$Treat), ]

## 2.3.  Visualization ####

Hill_coef <- c("q0", "q1")
rm(Micro_plots_ASV) # emptying the list before re-running
Micro_plots_ASV <- list()  # creating an empty list to store plots

for (var in Hill_coef) {
  
        avg_var <- paste0(var, "_avg")
        se_var <- paste0("se_", var, ".obs")
        obs_var <- paste0(var, ".obs")    
  
Micro_summary <- Hills_Micro_2023_ASV %>% 
        group_by(Treat) %>% 
        dplyr::summarise(!!avg_var := mean(.data[[obs_var]]),
                     !!se_var := sd(.data[[obs_var]]) / sqrt(n()))
  
Micro_allPlot_ASV <- ggplot(Hills_Micro_2023_ASV, aes_string("Treat", obs_var, group = "Treat", colour = "Treat", fill = "Treat")) +
        geom_point(position = position_jitterdodge (0.80, jitter.width = 0.2, jitter.height = 0), alpha = 0.2,shape = 21,colour = "black",size = 10) +
        scale_colour_manual(name = "Treatment", values = c("#002B5B", "#03C988", "#C31616")) +
        scale_fill_manual(values = c("#002B5B", "#03C988", "#C31616"), guide = "none") +
        theme_bw() +
        ylab(if (var == "q0") {expression("Species richness (q"[0]*")")}
             else if (var == "q1") {expression("Shannon diversity (q"[1]*")")}) +
        ggtitle("") +
        # annotate(geom="text", x=2, y=, label= "Sampled Group/s: \nSoil Microorganisms", color="black", size = 6, family = "serif", fontface = "bold") +
        theme(plot.title = element_text(size=20, hjust=0.5)) +
        theme(axis.title = element_text(size = 20), axis.text = element_text(size = 14), strip.text = element_text(size = 14),
              axis.title.y = element_text(size = 20, margin = margin(r = 1)), axis.title.x = element_blank(), legend.position = "none", 
              axis.text.y = element_text(size = 20, margin = margin(r = 0)), axis.text.x = element_text(size = 20), panel.border = element_rect(size = 1)) +
        geom_point(data = Micro_summary, aes_string(x = "Treat", y = avg_var), shape = 19, colour = "black", size = 12) +
        geom_point(data = Micro_summary, aes_string(x = "Treat", y = avg_var), shape = 19, size = 10) +
        geom_errorbar(data = Micro_summary, aes_string(x = "Treat", y = avg_var, ymin = paste0(avg_var, " - ", se_var), ymax = paste0(avg_var, " + ", se_var)), width = 0.3, size = 1) 
  
        Micro_plots_ASV[[paste0("Micro_", var)]] <- Micro_allPlot_ASV # adding plots to the list
  
}

pdf("outputs/Hill_Micro_ASV_q0.pdf", width = 5, height = 10)  
print(Micro_plots_ASV$Micro_q0) 
dev.off()

pdf("outputs/Hill_Micro_ASV_q1.pdf", width = 5, height = 10)
print(Micro_plots_ASV$Micro_q1) 
dev.off()

