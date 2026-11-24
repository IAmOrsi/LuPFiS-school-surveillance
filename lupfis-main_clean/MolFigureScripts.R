#####################################################Generating figures for Mol Epi results
###################################SARI DASHBOARD RELATIVE ABUNDANCE PLOTS
####calling data from csv file -  START HERE when running it de novo
setwd('~/LuPFiS/PCRanalyses')
Plotdata <- read.csv("Plotdata.csv")

####calling hospitalization data from file
SC2SARIdata <- read.csv("SC2SARIHospdata.csv")
RSVSARIdata <- read.csv("RSVSARIHospdata.csv")
FluABSARIdata <- read.csv("FluABSARIHospdata.csv")

####unify databases
names(RSVSARIdata)[names(RSVSARIdata) == "frequency"] <- "Hospitalized cases"
names(RSVSARIdata)[names(RSVSARIdata) == "midpoints"] <- "date.order.mod"

names(FluABSARIdata)[names(FluABSARIdata) == "frequency"] <- "Hospitalized cases"
names(FluABSARIdata)[names(FluABSARIdata) == "midpoints"] <- "date.order.mod"

####adjust observations from read in errors
##truncating sampling week data
SC2SARIdata$date.order.mod<- floor(SC2SARIdata$date.order.mod)
RSVSARIdata$date.order.mod<- floor(RSVSARIdata$date.order.mod)
FluABSARIdata$date.order.mod<- floor(FluABSARIdata$date.order.mod)

##cleaning up reading errors for RSV
RSVSARIdata$date.order.mod[11]<-17

missingRSVdata<-c(12, 'HOSP', 18, 0, 'Sampling weeks')
RSVSARIdata<-rbind(RSVSARIdata, missingRSVdata)
RSVSARIdata$id[is.na(RSVSARIdata$id)] <- 'HOSP'

##cleaning up reading errors for FluAB
FluABSARIdata$date.order.mod[1:3]<-FluABSARIdata$date.order.mod[1:3]+1
FluABSARIdata$date.order.mod[5:12]<-FluABSARIdata$date.order.mod[5:12]+1

FluABSARIdata$id[is.na(FluABSARIdata$id)] <- 'HOSP'
SC2SARIdata$id[is.na(SC2SARIdata$id)] <- 'HOSP'

SC2SARIdata$virus<-'SC2'
RSVSARIdata$virus<-'RSV'
FluABSARIdata$virus<-'Flu'

str(SC2SARIdata)
SC2SARIdata$virus<-as.factor(SC2SARIdata$virus)
SC2SARIdata$date.order.mod<-as.factor(SC2SARIdata$date.order.mod)

str(RSVSARIdata)
RSVSARIdata$virus<-as.factor(RSVSARIdata$virus)
RSVSARIdata$date.order.mod<-as.factor(RSVSARIdata$date.order.mod)
RSVSARIdata$`Hospitalized cases`<-as.integer(RSVSARIdata$`Hospitalized cases`)
names(RSVSARIdata)[names(RSVSARIdata) == "Hospitalized cases"] <- "Hospitalized.cases"

str(FluABSARIdata)
FluABSARIdata$virus<-as.factor(FluABSARIdata$virus)
FluABSARIdata$date.order.mod<-as.factor(FluABSARIdata$date.order.mod)
names(FluABSARIdata)[names(FluABSARIdata) == "Hospitalized cases"] <- "Hospitalized.cases"

figvars<-c('date.order.mod', 'Hospitalized.cases', 'virus')
SC2plotdata<-SC2SARIdata[figvars]
RSVplotdata<-RSVSARIdata[figvars]
Fluplotdata<-FluABSARIdata[figvars]

SARIdata<-rbind(SC2plotdata,RSVplotdata,Fluplotdata)

target_colours <- c(
  "SC2" = "#d31f11",   
  "RSV" = "#62c8d3",          
  "Flu" = "#f47a00"   # placeholder for stripes
)

SARIdata %>%
  group_by(date.order.mod) %>%
  mutate(rel_abund = Hospitalized.cases / sum(Hospitalized.cases)) %>%
  ggplot(aes(x = factor(date.order.mod), y = rel_abund, pattern = virus, fill = virus)) +
  geom_col_pattern(
    pattern = ifelse(SARIdata$virus == "Flu", "stripe", "none"),  # stripes only for Flu
    pattern_angle = 45,
    pattern_density = 0.5,
    pattern_spacing = 0.05,
    color = "black"
  ) +
  scale_fill_manual(values = target_colours) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),)
labs(
  x = "Sampling week",
  y = "Relative abundance",
  title = "Relative abundance from SARI dashboard"
)


###################################PATHOGEN ABUNDANCE - PCR DATA
setwd('~/LuPFiS/PCRanalyses')
Plotdata <- read.csv("Plotdata.csv")

####Excluding bacterial targets as they only serve as positive controls in the setting and typically overload 26k plates
Plotdata<-Plotdata[!grepl('Staphylococcus aureus', Plotdata$Target.Name), ]
Plotdata<-Plotdata[!grepl('Streptococcus pneumoniae', Plotdata$Target.Name), ]
Plotdata<- na.omit(Plotdata)

########################Center figure on abundance against ww and cumulative abundance for schools
####If running it on only environmental samples (air + swab) then run this code beforehand
Plotdata<-Plotdata[!grepl('T', Plotdata$sample.type), ]

####Averaging pathogen concentrations by sampling event
resultSumAbund <- Plotdata %>%
  group_by(date.order) %>%
  summarise(Average_Concentration = mean(Concentration), School_SE = sd(Concentration) / sqrt(n()))
DataSumAbund<-resultSumAbund

###Read in WW data
WWdata<-read.csv('WWdataMod.csv')
WWdata<- na.omit(WWdata)

####Averaging pathogen concentrations by sampling event
resultWWSumAbund <- WWdata %>%
  group_by(date.order) %>%
  summarise(Average_WW_Concentration = mean(sum_Salmon_TPM))
DataWWSumAbund<-resultWWSumAbund

##Merge databases for figure
PlotSumAbund_df<-left_join(resultSumAbund, resultWWSumAbund, by = "date.order") %>%
  rename(WW_Concentration = Average_WW_Concentration)
str(PlotSumAbund_df)
names(PlotSumAbund_df)[names(PlotSumAbund_df) == "Average_Concentration"] <- "School_Concentration"

##Generating the plot
scale_factor <- 1.5*(max(PlotSumAbund_df$School_Concentration, na.rm = TRUE) / max(PlotSumAbund_df$WW_Concentration, na.rm = TRUE))

bar_scale_factor <- 1 

max_val<-max(PlotSumAbund_df$School_Concentration+PlotSumAbund_df$School_SE, na.rm = TRUE)

SumAbund<-ggplot(PlotSumAbund_df, aes(x = date.order)) +
  geom_ribbon(aes(ymin = 0, ymax = WW_Concentration * scale_factor), 
              fill = "lightgrey", alpha = 0.3) +
  geom_line(aes(y = WW_Concentration * scale_factor), 
            color = "darkgrey", linewidth = 1) +
  geom_point(aes(y = WW_Concentration * scale_factor), 
             color = "darkgrey", size = 2) +
  geom_bar(aes(y = School_Concentration * bar_scale_factor), 
           stat = "identity", fill = "#4B4E6D", width = 0.6) +
  geom_errorbar(
    aes(ymin = (School_Concentration - School_SE) * bar_scale_factor,
      ymax = (School_Concentration + School_SE) * bar_scale_factor),
    width = 0.2, color = "black") +
  scale_y_continuous(
    name = "Pathogen abundance (school)",
    limits = c(0, max_val * 1.1),
    labels =  number_format(accuracy = 0.1),
    sec.axis = sec_axis(~ . / scale_factor, name = "Pathogen abundance (wastewater)")) +
  labs(x = "Sampling week", title = "Cumulative pathogen abundance measured in schools and wastewater") +
  theme_minimal(base_size = 14)



########################Figures for abundance targets in schools - RERUN DATABASE IMPORT WITHOUT !GREPL('T')
####Colour set for targets
target_colours <- c(
  "SARS-CoV-2" = "#d31f11",   
  "RSV" = "#62c8d3",          
  "Influenza A" = "#f47a00",  
  "Influenza B" = "#007191"   
)


####Subsetting data by school
##Groehrmuehlgasse
SchoolA_df<-Plotdata[ which(Plotdata$school %in% c('G')),]

max_total_conc_A <- SchoolA_df %>%
  filter(Target.Name %in% c("SARS-CoV-2", "RSV", "Influenza A", "Influenza B")) %>%
  group_by(date.order) %>%
  summarise(total_conc = sum(Concentration)) %>%
  pull(total_conc) %>%
  max(na.rm = TRUE)

scaling_factor_A <- max_total_conc_A  # scale total_conc down to max ~1 for plotting with relative_abundance

SchoolA<-SchoolA_df %>%
  filter(Target.Name %in% c("SARS-CoV-2", "RSV", "Influenza A", "Influenza B")) %>%
  group_by(date.order) %>%
  mutate(total_conc = sum(Concentration),
         relative_abundance = Concentration / total_conc) %>%
  ungroup() %>%
  ggplot(aes(x = date.order, y = relative_abundance, fill = Target.Name)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_line(aes(y = total_conc / scaling_factor_A, group = 1), color = "black", linewidth = 1) +
  geom_point(aes(y = total_conc / scaling_factor_A), color = "black", size = 2) +
  scale_fill_manual(values = target_colours) +
  scale_y_continuous(
    name = "Relative pathogen abundance",
    limits = c(0, 1.1),  # relative abundance max 1, give some headroom
    sec.axis = sec_axis(~ . * scaling_factor_A, name = "Cumulative pathogen abundance")
  ) +
  labs(
    x = "Sampling week",
    fill = "Target",
    title = "Stacked relative abundance of respiratory viruses with total concentration over time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 22, face = "bold"),
    axis.text.y = element_text(size = 22, face = "bold"),
  )

##Karl Popper Schule
SchoolB_df<-Plotdata[ which(Plotdata$school %in% c('K')),]

max_total_conc_B <- SchoolB_df %>%
  filter(Target.Name %in% c("SARS-CoV-2", "RSV", "Influenza A", "Influenza B")) %>%
  group_by(date.order) %>%
  summarise(total_conc = sum(Concentration)) %>%
  pull(total_conc) %>%
  max(na.rm = TRUE)

scaling_factor_B <- max_total_conc_B  # scale total_conc down to max ~1 for plotting with relative_abundance


SchoolB<-SchoolB_df %>%
  filter(Target.Name %in% c("SARS-CoV-2", "RSV", "Influenza A", "Influenza B")) %>%
  group_by(date.order) %>%
  mutate(total_conc = sum(Concentration),
         relative_abundance = Concentration / total_conc) %>%
  ungroup() %>%
  ggplot(aes(x = date.order, y = relative_abundance, fill = Target.Name)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_line(aes(y = total_conc / scaling_factor_B, group = 1), color = "black", linewidth = 1) +
  geom_point(aes(y = total_conc / scaling_factor_B), color = "black", size = 2) +
  scale_fill_manual(values = target_colours) +
  scale_y_continuous(
    name = "Relative pathogen abundance",
    limits = c(0, 1.1),  # relative abundance max 1, give some headroom
    sec.axis = sec_axis(~ . * scaling_factor_B, name = "Cumulative pathogen abundance")
  ) +
  labs(
    x = "Sampling week",
    fill = "Target",
    title = "Stacked relative abundance of respiratory viruses with total concentration over time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 22, face = "bold"),
    axis.text.y = element_text(size = 22, face = "bold"),
  )

##Montessori Campus
SchoolC_df<-Plotdata[ which(Plotdata$school %in% c('M')),]

max_total_conc_C <- SchoolC_df %>%
  filter(Target.Name %in% c("SARS-CoV-2", "RSV", "Influenza A", "Influenza B")) %>%
  group_by(date.order) %>%
  summarise(total_conc = sum(Concentration)) %>%
  pull(total_conc) %>%
  max(na.rm = TRUE)

scaling_factor_C <- max_total_conc_C  # scale total_conc down to max ~1 for plotting with relative_abundance


SchoolC<-SchoolC_df %>%
  filter(Target.Name %in% c("SARS-CoV-2", "RSV", "Influenza A", "Influenza B")) %>%
  group_by(date.order) %>%
  mutate(total_conc = sum(Concentration),
         relative_abundance = Concentration / total_conc) %>%
  ungroup() %>%
  ggplot(aes(x = date.order, y = relative_abundance, fill = Target.Name)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_line(aes(y = total_conc / scaling_factor_C, group = 1), color = "black", linewidth = 1) +
  geom_point(aes(y = total_conc / scaling_factor_C), color = "black", size = 2) +
  scale_fill_manual(values = target_colours) +
  scale_y_continuous(
    name = "Relative pathogen abundance",
    limits = c(0, 1.1),  # relative abundance max 1, give some headroom
    sec.axis = sec_axis(~ . * scaling_factor_C, name = "Cumulative pathogen abundance")
  ) +
  labs(
    x = "Sampling week",
    fill = "Target",
    title = "Stacked relative abundance of respiratory viruses with total concentration over time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 22, face = "bold"),
    axis.text.y = element_text(size = 22, face = "bold"),
  )


##Hagenmuellergasse
SchoolD_df<-Plotdata[ which(Plotdata$school %in% c('H')),]

max_total_conc_D <- SchoolD_df %>%
  filter(Target.Name %in% c("SARS-CoV-2", "RSV", "Influenza A", "Influenza B")) %>%
  group_by(date.order) %>%
  summarise(total_conc = sum(Concentration)) %>%
  pull(total_conc) %>%
  max(na.rm = TRUE)

scaling_factor_D <- max_total_conc_D  # scale total_conc down to max ~1 for plotting with relative_abundance


SchoolD<-SchoolD_df %>%
  filter(Target.Name %in% c("SARS-CoV-2", "RSV", "Influenza A", "Influenza B")) %>%
  group_by(date.order) %>%
  mutate(total_conc = sum(Concentration),
         relative_abundance = Concentration / total_conc) %>%
  ungroup() %>%
  ggplot(aes(x = date.order, y = relative_abundance, fill = Target.Name)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_line(aes(y = total_conc / scaling_factor_D, group = 1), color = "black", linewidth = 1) +
  geom_point(aes(y = total_conc / scaling_factor_D), color = "black", size = 2) +
  scale_fill_manual(values = target_colours) +
  scale_y_continuous(
    name = "Relative pathogen abundance",
    limits = c(0, 1.1),  # relative abundance max 1, give some headroom
    sec.axis = sec_axis(~ . * scaling_factor_D, name = "Cumulative pathogen abundance")
  ) +
  labs(
    x = "Sampling week",
    fill = "Target",
    title = "Stacked relative abundance of respiratory viruses with total concentration over time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 22, face = "bold"),
    axis.text.y = element_text(size = 22, face = "bold"),
  )

##Marchettigasse
SchoolE_df<-Plotdata[ which(Plotdata$school %in% c('R')),]

max_total_conc_E <- SchoolE_df %>%
  filter(Target.Name %in% c("SARS-CoV-2", "RSV", "Influenza A", "Influenza B")) %>%
  group_by(date.order) %>%
  summarise(total_conc = sum(Concentration)) %>%
  pull(total_conc) %>%
  max(na.rm = TRUE)

scaling_factor_E <- max_total_conc_E  # scale total_conc down to max ~1 for plotting with relative_abundance


SchoolE<-SchoolE_df %>%
  filter(Target.Name %in% c("SARS-CoV-2", "RSV", "Influenza A", "Influenza B")) %>%
  group_by(date.order) %>%
  mutate(total_conc = sum(Concentration),
         relative_abundance = Concentration / total_conc) %>%
  ungroup() %>%
  ggplot(aes(x = date.order, y = relative_abundance, fill = Target.Name)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_line(aes(y = total_conc / scaling_factor_E, group = 1), color = "black", linewidth = 1) +
  geom_point(aes(y = total_conc / scaling_factor_E), color = "black", size = 2) +
  scale_fill_manual(values = target_colours) +
  scale_y_continuous(
    name = "Relative pathogen abundance",
    limits = c(0, 1.1),  # relative abundance max 1, give some headroom
    sec.axis = sec_axis(~ . * scaling_factor_E, name = "Cumulative pathogen abundance")
  ) +
  labs(
    x = "Sampling week",
    fill = "Target",
    title = "Stacked relative abundance of respiratory viruses with total concentration over time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 22, face = "bold"),
    axis.text.y = element_text(size = 22, face = "bold"),
  )


##Wassermanngasse
SchoolF_df<-Plotdata[ which(Plotdata$school %in% c('W')),]

max_total_conc_F <- SchoolF_df %>%
  filter(Target.Name %in% c("SARS-CoV-2", "RSV", "Influenza A", "Influenza B")) %>%
  group_by(date.order) %>%
  summarise(total_conc = sum(Concentration)) %>%
  pull(total_conc) %>%
  max(na.rm = TRUE)

scaling_factor_F <- max_total_conc_F  # scale total_conc down to max ~1 for plotting with relative_abundance


SchoolF<-SchoolF_df %>%
  filter(Target.Name %in% c("SARS-CoV-2", "RSV", "Influenza A", "Influenza B")) %>%
  group_by(date.order) %>%
  mutate(total_conc = sum(Concentration),
         relative_abundance = Concentration / total_conc) %>%
  ungroup() %>%
  ggplot(aes(x = date.order, y = relative_abundance, fill = Target.Name)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_line(aes(y = total_conc / scaling_factor_F, group = 1), color = "black", linewidth = 1) +
  geom_point(aes(y = total_conc / scaling_factor_F), color = "black", size = 2) +
  scale_fill_manual(values = target_colours) +
  scale_y_continuous(
    name = "Relative pathogen abundance",
    limits = c(0, 1.1),  # relative abundance max 1, give some headroom
    sec.axis = sec_axis(~ . * scaling_factor_F, name = "Cumulative pathogen abundance")
  ) +
  labs(
    x = "Sampling week",
    fill = "Target",
    title = "Stacked relative abundance of respiratory viruses with total concentration over time"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 22, face = "bold"),
    axis.text.y = element_text(size = 22, face = "bold"),
  )

####Arrange figures in order

Phase1<- SchoolA / SchoolB / SchoolC
Phase2<- SchoolD / SchoolE / SchoolF

Phases<-Phase1 | Phase2

Final<- SumAbund / Phases

Final <- SumAbund / (Phase1 | Phase2) +
  plot_layout(heights = c(2, 3), widths = c(1, 1))
#-> not getting the hang of the patchwork to this detail, so final arrangement will be made in GIMP


###################################PATHOGEN DIVERSITY - SEQUENCE DATA
setwd('~/LuPFiS/Sequencing_analyses')
Data<-read_csv('SEQdata+sampledata.csv')
Data$species<-as.factor(Data$species)

#colour palette for targets
div_colours <- c(
  "Betacoronavirus 1" = "#913371",   
  "Human coronavirus 229E" = "#59c981",          
  "Human coronavirus NL63" = "#c55bb6",  
  "Human mastadenovirus B" = "#69b151",
  "Human metapneumovirus" = "#6973d8",
  "Rhinovirus C" = "#9dbb45",
  "Human orthorubulavirus 2" = "#573788",
  "Human respirovirus 3" = "#caa432",
  "Influenza B virus" = "#628dd4",
  "Primate bocaparvovirus 1" = "#c67d2d",
  "Respiratory syncytial virus" = "#c387d5",
  "Rhinovirus B" = "#607923",
  "SARS-related CoV" = "#dc6598",
  "Alphapolyomavirus quintihominis" = "#43c8ac",
  "Betapapillomavirus" = "#d85750",
  "Human coronavirus HKU1" = "#468646",
  "Human endogenous retrovirus K" = "#cd556e",
  "Human mastadenovirus C" = "#c0b15a",
  "Human orthopneumovirus" = "#8b264a",
  "Human respirovirus 1" = "#dd9360",
  "Influenza A virus" = "#893022",
  "Metapneumovirus sp." = "#8b6122",
  "Rabbit coronavirus HKU14" = "#cd535c",
  "Rhinovirus A" = "#c55530"
)


###########Number of counts per target overall
CountdataSP<-Data%>%count(species)
CountdataSP$species<-as.character(CountdataSP$species)
CountdataSP[CountdataSP=="Severe acute respiratory syndrome-related coronavirus"] <- "SARS-related CoV"

CountdataSP$species<-as.factor(CountdataSP$species)
CountdataSP$n<-as.numeric(CountdataSP$n)
str(CountdataSP)


ggplot(CountdataSP, aes(x = reorder(species, -n), y = n, fill=species)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = div_colours) +
  geom_text(aes(label = n), size = 3)+
  labs(
    title = "Cumulative species richness measured in schools",
    x = "Species",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = 'none')

###########Heatmap for breakdown by school
###Species level
SchoolheatmapdataSP<-Data%>%
  count(species, schools)

#Adding totals for targets and schools
target_totals <- SchoolheatmapdataSP %>%
  group_by(species) %>%
  summarise(n = sum(n)) %>%
  mutate(schools = "Total")

school_totals <- SchoolheatmapdataSP %>%
  group_by(schools) %>%
  summarise(n = sum(n)) %>%
  mutate(species = "Total")

wide_data <- SchoolheatmapdataSP %>%
  pivot_wider(names_from = schools, values_from = n, values_fill = 0)

wide_data$RowTotal <- rowSums(wide_data[,-1])
total_row <- c("Total", colSums(wide_data[,-1]))

# Convert total_row to dataframe and combine
total_row_df <- as.data.frame(t(total_row), stringsAsFactors = FALSE)
names(total_row_df) <- names(wide_data)
total_row_df[-1] <- lapply(total_row_df[-1], as.numeric)

combined <- bind_rows(wide_data, total_row_df)

plot_data <- combined %>%
  pivot_longer(-species, names_to = "schools", values_to = "n")

# Flag which cells are totals
plot_data <- plot_data %>%
  mutate(
    is_total = species == "Total" | schools == "Total",
    n = as.numeric(n)
  )

plot_data$schools <- factor(plot_data$schools, levels = c(setdiff(unique(plot_data$schools), "RowTotal"), "RowTotal"))

ggplot(plot_data, aes(x = schools, y = species)) +
  geom_tile(data = filter(plot_data, schools != "RowTotal" & species != "Total"),
            aes(fill = n),
            color = "white") + #Gradient only valid for non-total cells
  geom_tile(data = filter(plot_data, schools == "RowTotal" | species == "Total"),
            fill = "grey90",
            color = "white") + #Fixing the grey colour for all total cells
  geom_text(aes(label = n), size = 3) +  # <- Adds count to each cell
  scale_fill_gradientn(
    colors = c("white", "gold", "#A50021"),  # Low → Mid → High
    values = scales::rescale(c(0, median(SchoolheatmapdataSP$n), max(SchoolheatmapdataSP$n)))
  ) +
  labs(title = "Species frequency per school",
       x = "School", y = "Strain", fill = "Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold")
  )

###########Stacked barplot for richness for each sampling event in each school
###Species level

####Colour set for targets
str(Data$species)
summary(Data$species)

##Subsetting data for Groehrmuehlgasse
SchoolAdiv_df<-Data[ which(Data$schools %in% c('G')),]

SchoolAdiv_df$count<-1
result_SchoolAdiv_plot_df<-SchoolAdiv_df%>%
  group_by(Sampling_event, species) %>%
  summarise(count = count)

ggplot(result_SchoolAdiv_plot_df, aes(x = Sampling_event, y = count, fill = species)) +
  geom_bar(stat = "identity", width = 0.7) +
  ylim(0,24)+
  scale_fill_manual(values = div_colours) +
  labs(x = "Sampling event",
    fill = "Target") +
  theme_minimal(base_size = 14) +
  theme(legend.position = 'none', axis.text.y = element_text(face ='bold', size = 14), axis.text.x = element_text(face ='bold', size = 14))

##Subsetting data for Karl Popper Schule
SchoolBdiv_df<-Data[ which(Data$schools %in% c('K')),]

SchoolBdiv_df$count<-1
result_SchoolBdiv_plot_df<-SchoolBdiv_df%>%
  group_by(Sampling_event, species) %>%
  summarise(count = count)

ggplot(result_SchoolBdiv_plot_df, aes(x = Sampling_event, y = count, fill = species)) +
  geom_bar(stat = "identity", width = 0.7) +
  ylim(0,24)+
  scale_fill_manual(values = div_colours) +
  labs(x = "Sampling event",
       fill = "Target") +
  theme_minimal(base_size = 14) +
  theme(legend.position = 'none', axis.text.y = element_text(face ='bold', size = 14), axis.text.x = element_text(face ='bold', size = 14))

##Subsetting data for Montessori Campus
SchoolCdiv_df<-Data[ which(Data$schools %in% c('M')),]

SchoolCdiv_df$count<-1
result_SchoolCdiv_plot_df<-SchoolCdiv_df%>%
  group_by(Sampling_event, species) %>%
  summarise(count = count)

ggplot(result_SchoolCdiv_plot_df, aes(x = Sampling_event, y = count, fill = species)) +
  geom_bar(stat = "identity", width = 0.7) +
  xlim(0.6,8)+
  ylim(0,24)+
  scale_fill_manual(values = div_colours) +
  labs(x = "Sampling event",
       fill = "Target") +
  theme_minimal(base_size = 14) +
  theme(legend.position = 'none', axis.text.y = element_text(face ='bold', size = 14), axis.text.x = element_text(face ='bold', size = 14))


##Subsetting data for Hagenmuellergasse
SchoolDdiv_df<-Data[ which(Data$schools %in% c('H')),]

SchoolDdiv_df$count<-1
result_SchoolDdiv_plot_df<-SchoolDdiv_df%>%
  group_by(Sampling_event, species) %>%
  summarise(count = count)

ggplot(result_SchoolDdiv_plot_df, aes(x = Sampling_event, y = count, fill = species)) +
  geom_bar(stat = "identity", width = 0.7) +
  xlim(0.6,8)+
  ylim(0,24)+
  scale_fill_manual(values = div_colours) +
  labs(x = "Sampling event",
       fill = "Target") +
  theme_minimal(base_size = 14) +
  theme(legend.position = 'none', axis.text.y = element_text(face ='bold', size = 14), axis.text.x = element_text(face ='bold', size = 14))

##Subsetting data for Marchettigasse
SchoolEdiv_df<-Data[ which(Data$schools %in% c('R')),]

SchoolEdiv_df$count<-1
result_SchoolEdiv_plot_df<-SchoolEdiv_df%>%
  group_by(Sampling_event, species) %>%
  summarise(count = count)

ggplot(result_SchoolEdiv_plot_df, aes(x = Sampling_event, y = count, fill = species)) +
  geom_bar(stat = "identity", width = 0.7) +
  ylim(0,24)+
  scale_fill_manual(values = div_colours) +
  labs(x = "Sampling event",
       fill = "Target") +
  theme_minimal(base_size = 14) +
  theme(legend.position = 'none', axis.text.y = element_text(face ='bold', size = 14), axis.text.x = element_text(face ='bold', size = 14))

##Subsetting data for Wassermanngasse
SchoolFdiv_df<-Data[ which(Data$schools %in% c('W')),]

SchoolFdiv_df$count<-1

#to unify x axis values, we change this here
result_SchoolFdiv_plot_df<-SchoolFdiv_df%>%
  group_by(Sampling_event, species) %>%
  summarise(count = count)

#to this
str(result_SchoolFdiv_plot_df)
result_SchoolFdiv_plot_df <- result_SchoolFdiv_plot_df %>%
  ungroup()%>%
  mutate(Sampling_event = factor(Sampling_event, levels = as.character(1:8))) %>%
  tidyr::complete(Sampling_event, species, fill = list(count = 0))

ggplot(result_SchoolFdiv_plot_df, aes(x = Sampling_event, y = count, fill = species)) +
  geom_bar(stat = "identity", width = 0.7) +
  ylim(0,24)+
  scale_fill_manual(values = div_colours) +
  labs(x = "Sampling event",
       fill = "Target") +
  theme_minimal(base_size = 14) +
  theme(legend.position = 'none', axis.text.y = element_text(face ='bold', size = 14), axis.text.x = element_text(face ='bold', size = 14))


###################################SAMPLE COMPARISON - ABUNDANCE
####setting data
setwd('~/LuPFiS/PCRanalyses')
Plotdata <- read_csv("Plotdata.csv")
Plotdata$date.order<-as.numeric(Plotdata$date.order)
Plotdata$Concentration<-as.numeric(Plotdata$Concentration)

Plotdata<-na.omit(Plotdata)

Plotdata<-Plotdata[!grepl('Staphylococcus aureus', Plotdata$Target.Name), ]
Plotdata<-Plotdata[!grepl('Streptococcus pneumoniae', Plotdata$Target.Name), ]

###Creating df for plot, log transform not used in plot but i didn't want to rewrite the code
Plotdata_conc <- Plotdata %>%
  filter(Concentration > 0)

#running twwedie glm to produce residulas
SampleComp_Abund_Conc_mod <- glm(Concentration ~ sample.type,
                                 family = tweedie(var.power = 1.5, link.power = 0),
                                 data = Plotdata_conc)
summary(SampleComp_Abund_Conc_mod)

##extracting fit and upper/lower limits of each groupd of sample type into separate df
Tweedie_plot_df<-data.frame(sample.type = factor(unique(Plotdata_conc$sample.type),
                                                 levels = unique(Plotdata_conc$sample.type)))

##calculating predicted means for each group
predicted<-predict(SampleComp_Abund_Conc_mod, newdata = Tweedie_plot_df, type = 'link', se.fit = TRUE)

##adding fit, upper and lower limit to predicted means
Tweedie_plot_df$fit <- exp(predicted$fit)  # inverse of log link
Tweedie_plot_df$lower <- exp(predicted$fit - 1.96 * predicted$se.fit)
Tweedie_plot_df$upper <- exp(predicted$fit + 1.96 * predicted$se.fit)

#calculating estimated means for all sample types 
emm_abund<-emmeans(SampleComp_Abund_Conc_mod, specs = 'sample.type')

##pairwise comparison of groups and defining significance levels for later plotting
pairs<- pairs(emm_abund)

pairs_df<-as.data.frame(summary(pairs))

pairs_df$signif<-cut(pairs_df$p.value,
                     breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                     labels = c("***", "**", "*", "ns"))

sig_df <- data.frame(
  group1 = pairs_df$contrast %>% sapply(function(x) strsplit(x, " - ")[[1]][1]),
  group2 = pairs_df$contrast %>% sapply(function(x) strsplit(x, " - ")[[1]][2]),
  p.value = pairs_df$p.value,
  y.position = c(1600, 1800, 2000),
  labels = pairs_df$signif
)

#plotting
ggplot(Tweedie_plot_df, aes(x = sample.type, y = fit)) +
  geom_col(fill = "#4B4E6D") +
  geom_errorbar(aes(ymin = lower , ymax = upper), width = 0.2) +
  stat_pvalue_manual(sig_df, label = "labels", size = 8, bracket.size = 1 , tip.length = 0.01)+
  labs(y = "Predicted pathogen abundance", x = "Sample Type",
       title = "Predicted Pathogen Abundance by Sample Type") +
  theme(axis.text.y = element_text(face="bold", size=14), axis.text.x = element_text(face='bold', size = 14))

###################################SAMPLE COMPARISON - DIVERSITY
####setting data
setwd("~/LuPFiS/Sequencing_analyses")
Data<-read_csv('SEQdata+sampledata.csv')
divvarsSP<-c('schools', 'Sample_type', 'Sampling_event', 'species')
DivDataSP<-Data[divvarsSP]

###Create species richness for each sampling event in each school in a new df
DivDataSP<-DivDataSP%>%
  mutate(count = 1)
richness_df<- DivDataSP %>%
  group_by(schools, Sampling_event, Sample_type) %>%
  summarise(richness = sum(count), .groups = "drop")

#Remove NA cases
richness_df_NAOMIT<-na.omit(richness_df)

#Add order to the df
richness_df_NAOMIT$order<-with(richness_df_NAOMIT, ifelse(schools %in% c("G", "M", "K"),
                                                          Sampling_event, Sampling_event+8))

#running model to test for effects
richness_sampletype_Poisson_mod_full<-glm(richness~Sample_type+schools+order, data= richness_df_NAOMIT, family = poisson)
summary(richness_sampletype_Poisson_mod_full)

##running posthoc test for sample type effect
emm_sample_full<-emmeans(richness_sampletype_Poisson_mod_full,~Sample_type)
plot(emm_sample_full)
pairs(emm_sample_full, adjust='tukey')

emm_sample_full_df<-as.data.frame(emm_sample_full)
ggplot(emm_sample_full_df, aes(x = Sample_type, y = emmean)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(y = "Predicted Richness", x = 'Sample type', title = "Predicted Pathogen Richness by Sample Type") +
  theme(axis.text.y = element_text(face="bold", size=14), axis.text.x = element_text(face='bold', size = 14))

###################################CO2 EFFECTS- ABUNDANCE
####setting CO2 data
setwd('~/LuPFiS/CO2')
CO2data <- read_csv("CO2_summary_data.csv")

names(CO2data)[names(CO2data) == "start_timestamp"] <- "starttime"
names(CO2data)[names(CO2data) == "end_timestamp"] <- "endtime"
CO2data<-CO2data%>%
  mutate(starttime=gsub("T|Z", " ", starttime))
CO2data<-CO2data%>%
  mutate(endtime=gsub("T|Z", " ", endtime))

#only overall CO2 ranges for 24h time period
CO2_overall_df<-CO2data[ which(CO2data$timeslot %in% c('overall')),]

CO2_overall_df<-CO2_overall_df%>%
  mutate(starttime = ymd_hms(starttime)) %>%
  group_by(school) %>%
  arrange(starttime, .by_group = TRUE)%>%
  mutate(date.order = row_number())

##school data
setwd('~/LuPFiS/PCRanalyses')
Plotdata <- readr::read_csv("Plotdata.csv")

#excluding tissues due to significant effect + less likely contact with CO2 levels
Plotdata<-Plotdata[!grepl('T', Plotdata$sample.type), ]

Plotdata<-Plotdata[!grepl('Staphylococcus aureus', Plotdata$Target.Name), ]
Plotdata<-Plotdata[!grepl('Streptococcus pneumoniae', Plotdata$Target.Name), ]

Plotdata$school<-as.factor(Plotdata$school)

#calculating average concentration from environmental samples
resultPCR <- Plotdata %>%
  group_by(school, date.order, Target.Name) %>%
  summarise(Average_Concentration = mean(Concentration))
Plotdata_df<-resultPCR

Plotdata_transposed_df <- Plotdata_df %>%
  arrange(school, date.order) %>%
  pivot_wider(
    names_from = Target.Name,
    values_from = Average_Concentration
  )

#merging CO2 with abundance data
merged_CO2_df <- full_join(Plotdata_transposed_df, CO2_overall_df, by = c("school", "date.order")) %>%
  arrange(school, date.order)

merged_CO2_df <- merged_CO2_df %>%
  mutate(PathAbundPH = rowSums(across(c("Influenza A", "Influenza B",
                                        "RSV", "SARS-CoV-2")), na.rm = TRUE))
merged_CO2_df_naomit<-na.omit(merged_CO2_df)
str(merged_CO2_df_naomit)
merged_CO2_df_naomit$school<-as.factor(merged_CO2_df_naomit$school)

merged_CO2_df_naomit <- merged_CO2_df_naomit %>%
  mutate(
    School_code = case_when(
      school == "G" ~ "School A",  
      school == "M" ~ "School C", 
      school == "K" ~ "School B", 
      school == "R" ~ "School E",
      school == "H" ~ "School D",
      school == "W" ~ "School F",
      TRUE ~ NA_character_          
    )
  )

###plotting median
PathAbundPH_CO2_mod<-glm(I(PathAbundPH + 1e-3)~co2.ppm.median+School_code, data= merged_CO2_df_naomit, family = Gamma(link='log'))
summary(PathAbundPH_CO2_mod)

PathAbundPH_CO2_PLOT_mod<-glm(I(PathAbundPH + 1e-3)~School_code, data= merged_CO2_df_naomit, family = Gamma(link='log'))
PlotPathAbundPH_medianCO2 <- data.frame(
  residuals = residuals(PathAbundPH_CO2_PLOT_mod, type = "deviance")
)
PlotPathAbundPH_medianCO2$CO2median<-merged_CO2_df_naomit$co2.ppm.median
PlotPathAbundPH_medianCO2$School_code<-merged_CO2_df_naomit$School_code

ggplot(PlotPathAbundPH_medianCO2, aes(x = CO2median, y = residuals, color=School_code)) +
  geom_point(alpha = 0.6, size=2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 35), se = TRUE, color = "blue") +
  annotate("text", x = 600, y = 2, label = paste("p =", signif(summary(PathAbundPH_CO2_mod)$coefficients[2,4], 3)))+
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5))+
  labs(x = "Median of daily CO2 levels (ppm)",
    y = "Pathogen abundance (model residuals)",
    title = "Pathogen abundance by median CO2 levels")

###plotting high risk
PathAbundPH_CO2_high_mod<-glm(I(PathAbundPH + 1e-3)~threshold.high.duration.MINUTES+School_code, data= merged_CO2_df_naomit, family = Gamma(link='log'))
summary(PathAbundPH_CO2_high_mod)

PlotPathAbundPH_medianCO2$highdur<-merged_CO2_df_naomit$threshold.high.duration.MINUTES

ggplot(PlotPathAbundPH_medianCO2, aes(x = highdur, y = residuals, color=School_code)) +
  geom_point(alpha = 0.6, size=4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 20), se = TRUE, color = "black") +
  annotate("text", x = 500, y = 2, label = paste("p =", signif(summary(PathAbundPH_CO2_high_mod)$coefficients[2,4], 3)))+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),)+
  labs(x = "Duration of daily high-risk CO2 levels (min)",
    y = "Pathogen abundance (model residuals)",
    title = "Pathogen abundance by duration of high-risk CO2 levels")

###plotting area under the curve
#calling data from different analyses
setwd('~/LuPFiS/CO2')
DIV_CO2auc_df<-read_csv('Richness_CO2auc.csv')
DIV_CO2auc_df_naomit<-na.omit(DIV_CO2auc_df)

setwd('~/LuPFiS/PCRanalyses')
Plotdata <- read.csv("Plotdata.csv")

Plotdata<-Plotdata[!grepl('T', Plotdata$sample.type), ]

Plotdata<-Plotdata[!grepl('Staphylococcus aureus', Plotdata$Target.Name), ]
Plotdata<-Plotdata[!grepl('Streptococcus pneumoniae', Plotdata$Target.Name), ]

resultPCR <- Plotdata %>%
  group_by(school, date.order, Target.Name) %>%
  summarise(Average_Concentration = mean(Concentration))
Plotdata_df<-resultPCR

Plotdata_df$order<-with(Plotdata_df, ifelse(school %in% c("G", "M", "K"),
                                            date.order, date.order+8))

Plotdata_transposed_df <- Plotdata_df %>%
  arrange(school, date.order) %>%
  pivot_wider(
    names_from = Target.Name,
    values_from = Average_Concentration)

#Unifying columns between auc df and abund df for merging
Abund_df <- Plotdata_transposed_df %>%
  mutate(PathAbund = rowSums(across(c("Adenovirus", "Enterovirus", "Influenza A", "Influenza B",
                                      "Norovirus", "RSV", "Rhinovirus", "SARS-CoV-2")), na.rm = TRUE))
Abund_df <- Abund_df %>%
  select(-order) %>%              
  rename(order = date.order)      

Abund_df_merge<-Abund_df%>%
  select(school, order, PathAbund)

DIV_CO2auc_df<-DIV_CO2auc_df%>%
  left_join(Abund_df_merge, by = c('school', 'order'))

DIV_CO2auc_df_naomit<-na.omit(DIV_CO2auc_df)

DIV_CO2auc_df_naomit <- DIV_CO2auc_df_naomit %>%
  mutate(
    School_code = case_when(
      school == "G" ~ "School A",  
      school == "M" ~ "School C", 
      school == "K" ~ "School B", 
      school == "R" ~ "School E",
      school == "H" ~ "School D",
      school == "W" ~ "School F",
      TRUE ~ NA_character_))

###PathAbund is still Poisson from the previous analyses, so model uses poisson family
##Viral diversity analyses revealed sample type has no effect on richness, only school, so that is the only co-variant to correct model
CO2_abund_auc_mod<-glm(I(PathAbund + 1e-3)~auc+School_code, data= DIV_CO2auc_df_naomit, family = poisson(link='log'))
summary(CO2_abund_auc_mod)

##Plot significant result for auc
CO2_abund_auc_PLOT_mod<-glm(I(PathAbund + 1e-3)~School_code, data= DIV_CO2auc_df_naomit, family = poisson(link='log'))
PlotCO2_Abund_auc <- data.frame(
  residuals = residuals(CO2_abund_auc_PLOT_mod, type = "deviance"))

PlotCO2_Abund_auc$auc<-DIV_CO2auc_df_naomit$auc
PlotCO2_Abund_auc$School_code<-DIV_CO2auc_df_naomit$School_code
ggplot(PlotCO2_Abund_auc, aes(x = auc, y = residuals, color=School_code)) +
  geom_point(alpha = 0.6, size=4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 24), se = TRUE, color = "black") +
  annotate("text", x = 500, y = 2, label = paste("p =", signif(summary(CO2_abund_auc_mod)$coefficients[2,4], 3)))+
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),)+
  labs(x = "Area under the CO2 curve",
    y = "Pathogen abundance (model residuals)",
    title = "Pathogen abundance by exposure to CO2")

###For graphical abstract
ggplot(PlotCO2_Abund_auc, aes(x = auc, y = residuals)) +
  geom_point(alpha = 0.6, size=4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 24), se = TRUE, color = "black") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),)+
  labs(x = "CO2 level",
       y = "Pathogen abundance")

###################################CO2 EFFECTS- DIVERSITY
####setting CO2 data
setwd('~/LuPFiS/Sequencing_analyses')
DataDIVCO2<-read_csv('Viral_diversity_species_richness.csv')
str(DataDIVCO2)
DataDIVCO2$school<-as.factor(DataDIVCO2$school)

DataDIVCO2 <- DataDIVCO2 %>%
  mutate(
    School_code = case_when(
      school == "G" ~ "School A",  
      school == "M" ~ "School C", 
      school == "K" ~ "School B", 
      school == "R" ~ "School E",
      school == "H" ~ "School D",
      school == "W" ~ "School F",
      TRUE ~ NA_character_))


###plotting median
richness_df_naomit<-na.omit(DataDIVCO2)
richness_median_mod<-glm(richness~Median_CO2+School_code, data= richness_df_naomit, family = Gamma(link='log'))
summary(richness_median_mod)

##Plot non-significant effects for median
CO2_div_med_PLOT_mod<-glm(richness~School_code, data= richness_df_naomit, family = Gamma(link='log'))
PlotCO2_div_med <- data.frame(
  residuals = residuals(CO2_div_med_PLOT_mod, type = "deviance"))

PlotCO2_div_med$Median_CO2<-richness_df_naomit$Median_CO2
PlotCO2_div_med$School_code<-richness_df_naomit$School_code
ggplot(PlotCO2_div_med, aes(x = Median_CO2, y = residuals, color=School_code)) +
  geom_point(alpha = 0.6, size=2) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 24), se = TRUE, color = "blue") +
  theme_minimal() +
  labs(x = "Median CO2 level",
       y = "Pathogen richness (model residuals)",
       title = "Pathogen diversity by median CO2 levels")

###plotting high risk
richness_high_mod<-glm(richness~Duration_high+School_code, data= richness_df_naomit, family = Gamma(link='log'))
summary(richness_high_mod)

##Plot non-significant effects for high-risk
CO2_div_high_PLOT_mod<-glm(richness~School_code, data= richness_df_naomit, family = Gamma(link='log'))
PlotCO2_div_high <- data.frame(
  residuals = residuals(CO2_div_high_PLOT_mod, type = "deviance"))

PlotCO2_div_high$Duration_high<-richness_df_naomit$Duration_high
PlotCO2_div_high$School_code<-richness_df_naomit$School_code
ggplot(PlotCO2_div_high, aes(x = Duration_high, y = residuals, color=School_code)) +
  geom_point(alpha = 0.6, size=4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10), se = TRUE, color = "blue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),)
  labs(x = "Duration of high-risk CO2 levels",
       y = "Pathogen richness (model residuals)",
       title = "Pathogen diversity by duration of high-risk CO2 levels")

###plotting auc
setwd('~/LuPFiS/CO2')
DIV_CO2auc_df<-read_csv('Richness_CO2auc.csv')
DIV_CO2auc_df_naomit<-na.omit(DIV_CO2auc_df)

DIV_CO2auc_df_naomit <- DIV_CO2auc_df_naomit %>%
  mutate(
    School_code = case_when(
      school == "G" ~ "School A",  
      school == "M" ~ "School C", 
      school == "K" ~ "School B", 
      school == "R" ~ "School E",
      school == "H" ~ "School D",
      school == "W" ~ "School F",
      TRUE ~ NA_character_))

CO2_auc_mod<-glm(richness~auc+School_code, data= DIV_CO2auc_df_naomit, family = poisson(link='log'))
summary(CO2_auc_mod)

##Plot non-significant effects for high-risk
CO2_div_auc_PLOT_mod<-glm(richness~School_code, data= DIV_CO2auc_df_naomit, family = poisson(link='log'))
PlotCO2_div_auc <- data.frame(
  residuals = residuals(CO2_div_auc_PLOT_mod, type = "deviance"))

PlotCO2_div_auc$auc<-DIV_CO2auc_df_naomit$auc
PlotCO2_div_auc$School_code<-DIV_CO2auc_df_naomit$School_code
ggplot(PlotCO2_div_auc, aes(x = auc, y = residuals, color=School_code)) +
  geom_point(alpha = 0.6, size=4) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10), se = TRUE, color = "blue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),)
  labs(x = "Area under the CO2 curve",
       y = "Pathogen richness (model residuals)",
       title = "Pathogen diversity by exposure to CO2")

###################Recreating figure from Nature article for Fig 5.
setwd('~/LuPFiS/Figures')
Fig5data<-read.csv("Fig5Data.csv", sep=";")

quad_colors <- c("#fde0dd", "#e0ecf4", "#e6f5d0", "#fff7bc")

cat_colors <- c("#8B0000", "#00008B", "#006400", "#FF8C00")


ggplot(Fig5data, aes(x = Accessibility, y = Coverage, color = Category, label = Data)) +
  annotate("rect", xmin = -Inf, xmax = 5, ymin = -Inf, ymax = 10, fill = quad_colors[1], alpha = 0.3) +
  annotate("rect", xmin = 5, xmax = Inf, ymin = -Inf, ymax = 10, fill = quad_colors[2], alpha = 0.3) +
  annotate("rect", xmin = -Inf, xmax = 5, ymin = 10, ymax = Inf, fill = quad_colors[3], alpha = 0.3) +
  annotate("rect", xmin = 5, xmax = Inf, ymin = 10, ymax = Inf, fill = quad_colors[4], alpha = 0.3) +
  geom_point(size = 3) +
  geom_text_repel(show.legend = FALSE) +  # avoids overlapping text
  scale_color_manual(values = cat_colors) +
  theme_minimal() +
  theme(axis.text = element_blank(),axis.ticks = element_blank()) +
  labs(x = "Accessibility",
    y = "Coverage",
    color = "Data types")

#####Figure 5 as discussed on lab retreat with Greg
Fig5data<-read.csv("Fig5dataV2.csv", sep=";")

quad_colors <- c("#fde0dd", "#e0ecf4", "#e6f5d0", "#fff7bc")

cat_colorsV2 <- c("#005f56", "#8c3600", "#3f3a7a", "#8c1c52","#e31a1c", "#3d5c0d", "#8c6d00")

# Reorder factor levels
Fig5data$Surveillance.source <- factor(Fig5data$Surveillance.source,
                            levels = c("Transdisciplinary epidemiology",
                                       "Census",
                                       "Online (searches + mobility)",
                                       "Patient data",
                                       "Prevalence surveys",
                                       "Vaccination registry",
                                       "Wastewater surveillance"))

# Keep the same color mapping (red stays with Transdisciplinary epidemiology)
cat_colorsV2_assign <- c("Transdisciplinary epidemiology" = "#e31a1c",
                "Census" = "#005f56",
                "Online (searches + mobility)" = "#8c3600",
                "Patient data" = "#3f3a7a",
                "Prevalence surveys" = "#8c1c52",
                "Vaccination registry" = "#3d5c0d",
                "Wastewater surveillance" = "#8c6d00")

x_mid<-6
y_mid<-11

ggplot(Fig5data, aes(x = Data.quality, y = Coverage.for.high.risk, color = Surveillance.source, shape = Type)) +
  annotate("rect", xmin = -Inf, xmax = x_mid, ymin = -Inf, ymax = y_mid, fill = quad_colors[1], alpha = 0.3) +
  annotate("rect", xmin = x_mid, xmax = Inf, ymin = -Inf, ymax = y_mid, fill = quad_colors[2], alpha = 0.3) +
  annotate("rect", xmin = -Inf, xmax = x_mid, ymin = y_mid, ymax = Inf, fill = quad_colors[3], alpha = 0.3) +
  annotate("rect", xmin = x_mid, xmax = Inf, ymin = y_mid, ymax = Inf, fill = quad_colors[4], alpha = 0.3) +
  geom_point(size = 4) +
  guides(color = guide_legend(override.aes = list(shape = 15))) +
  scale_color_manual(values = cat_colorsV2_assign) +
  scale_shape_manual(values = c(16, 4, 17)) +
  theme_minimal() +
  theme(axis.text = element_blank(),axis.ticks = element_blank(), legend.position = "right", 
        axis.title.x = element_text(size = 18, face = "bold"),  
        axis.title.y = element_text(size = 18, face = "bold")) +
  labs(x = "Data quality",
       y = "Coverage for high risk",
       color = "Surveillance source", shape = "Data types")

str(Fig5data)
Fig5data$Data.quality<-as.numeric(Fig5data$Data.quality)
Fig5data$Data.quality[2] <- "8.5"
Fig5data$Data.quality[5] <- "8.5"
Fig5data$Data.quality[12] <- "10"

Fig5data$Coverage.for.high.risk<-as.numeric(Fig5data$Coverage.for.high.risk)
Fig5data$Coverage.for.high.risk[8]<-"18.5"

Fig5data$Surveillance.source<-as.factor(Fig5data$Surveillance.source)
Fig5data$Type<-as.factor(Fig5data$Type)
