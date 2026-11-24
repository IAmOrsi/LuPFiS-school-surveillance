# Data and code from the LuPFiS transdiscplinary epidemiological study

Data and code are organized into subfolders named after the analyses they hold information on. All subfolders are in the lupfis-main_clean folder, representing the source branch the information is shared from. This is a clean, anonymized version of the code, meaning school identifiers have been removed from the code.
Code for reproducing molecular analyses and Fig 2DEFG and Fig 5 is stored in Analyses_Sample_CO2_Figure2DEFG_Fig5.R found in the main folder. Code for reproducing Fig 2ABC is stored in Fig2ABC.md

## Analysis Reports
This contains the rendered html files explaining all steps of statistical analyses condcted with molecular data, with code made visible. Files are separated based on steps (Descriptive, Exploratory, Sample effect testing, CO2 effect testing) and whetehr the dependent variable was Pathoen Abundance (Abundance) or Richness (Diversity)

## CO2
This folder holds data collected from CO2 measruing device and is required for all analyses conducted with CO2 levels
Code for calculating CO2 metrics and analyses is stored Analysis Reports

## PCR Analyses
This holds data for collected from PCR analyses and is required for all analyses conducted with Pathogen abundance.
Code using this data source is stored in Analysis Reports and Analyses_Sample_CO2_Figure2DEFG_Fig5.R 

## Sequencing analyses
This holds data collected from Sequencing pipeline, raw sequence data is stored on European Nucleotide Archive (ENA) database under accession PRJEB130317, code for processing sequence data is currently under publication and will be shared upon request.
