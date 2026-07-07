## Octocoral acute warming
This repository contains data and code used to compile the manuscript: Rakka M, Metaxas A, Bilan M. Rapid response of two deep-water octocorals to acute warming reveals their vulnerability to extreme warming events.
The repository contains three folders with code, data and outputs

## Code

### Extr_warming_part1_survey.R
Here we use data collected during a field survey aiming at describing the health status of the two target species. The code creates figures and uses ordinal regression to analyze the results.

### Extr_warming_part2_experiment.R
Here we analyze data collected during a short experiment aiming at determining the thermal tolerance of the two target species. The code starts by creating tables of temperature and salinity recorded during the
experiment, and then analyzes collected data on polyp activity, tissue loss, oxygen consumption and Structural Area Index which is a proxy of branch area.

## Data

### DailyMonitoring.csv
The file contains observations of temperature, salinity and the state of coral fragments, collected throughout the experiment.

-Species (character): The coral species, Primnoa for Primnoa resedaeformis, Paragorgia for Paragorgia arborea
-mydate (date): Includes information on the date and time of the observation
-expTime (numeric): Time since the beginning of the experiment in hours
-Treatment (numeric): Target temperature of the experimental treatment, one of 3,5,8,12
-Flask (character): The code of the flask where the coral fragment was maintained
-colony (character): The code of the colony that the fragment belonged to
-Temperature (numeric): Measured temperature in C
-Salinity (numeric): Measured salinity
-Polyps (character): Whether the polyps of the colony were open or closed
-opencols (binary): binary version of the Polyps variable (o=closed, 1=open)

### alldat.csv
The file contains all other measured variables collected at the end of the experiment, including mortality, Tissue loss, oxygen consumption, as well as measurements of the coral fragments.
-Species (character): The coral species, Primnoa for Primnoa resedaeformis, Paragorgia for Paragorgia arborea
-Treatment (numeric): Target temperature of the experimental treatment, one of 3,5,8,12
-Flask (character): The code of the flask where the coral fragment was maintained
-colony (character): The code of the colony that the fragment belonged to
-frag (character): The code of the coral fragment
-coral code (character): The complete code of the coral fragment
-Mortality_48h (binary): whether the fragment was alive (0) or dead (1)
-Tissue_loss_48h (ordered factor): Describes the presence of tissue loss, none (healthy), mild tissue loss, severe tissue loss     
-Ox_consumption_volh_dw (numeric): Oxygen consumption in mg/L per coral dry weight (g)
-Surface (numeric): The total surface of the coral fragment in squared meters          
-Volume (numeric): The total volume of the coral fragment in cubic meters
-Branch_length (numeric): The total branch length of the coral fragment in m
-polyp_nb (numeric): The total number of polyps in each fragment               
-tis_index (numeric): structural area index (SAI), estimated as the total area divided by the branch length

### CoralHealthField.csv
This file contains data collected during video surveys. Each row represents a coral colony. The file includes the following variables:
-Fishing_lines (binary): Describes the presence (1) or absence (0) of fishing lines in the video frame that contained the colony
-Sediment (binary): Describes the presence (1) or absence (0) of sediment on the colony
-Hydroids (binary): Describes the presence (1) or absence (0) of epibiotic hydroids on the colony
-dive (character): the code of the dive        
-Species (character): The coral species, Primnoa for Primnoa resedaeformis, Paragorgia for Paragorgia arborea
-Condition (ordinal factor): Condition of the coral colony, healthy, mild tissue loss, severe tissue loss
-Year (numeric): The year of the video survey        
-location (character): The location of the video survey
