#clear global environment

rm(list = ls())

#setwd("Please insert your working directory filepath here")  # (Example working directory below)

setwd()


# Install and load Boruta with dependency 'Ranger' and tidyverse package
install.packages('Boruta')
install.packages('ranger')
install.packages('tidyverse')

library(Boruta)
library(ranger)
library(tidyverse)

# Upload cell line data from your working directory
AB.data <- read.csv("CNV loss 1.csv")

# Execute Boruta algorithm, focusing on response (to rucaparib)
AB.RucClono.Boruta <- Boruta(response~.,AB.data, maxRun =5000,pValue=0.05)

# View basic results
print(AB.RucClono.Boruta)

#to get full list of important attributes

getSelectedAttributes(AB.MiST1.Boruta)

# Visualise results with a boxplot (significant interactions coloured green)
plot(AB.RucClono.Boruta)


#### Repeat Boruta execution up to 5,000 iterations with aooropriate p value e.g.p=0.005
successes <- 0
attempts <- 0
Boruta5k <- repeat{
  repeat{
    attempts <- attempts + 1
    results <- Boruta(response ~., AB.data, maxRuns=5000, pvalue=0.0001)
    duration <- as.numeric(results$timeTaken, units='secs')
    print(paste0('Attempt ', attempts, ' lasted ', duration, ' seconds.'))
    if (duration > 70){
      print("Sufficient iterations reached");
      successes <- successes + 1
      imps.stats <- attStats(results) %>% filter(decision=='Confirmed') %>% 
        mutate('Execution'=successes) 
      imps.stats <- imps.stats %>% mutate('Variable'=rownames(imps.stats))
      print(getSelectedAttributes(results))
      print(results)
      break
    }
  }
  ifelse(successes==1, collated <- imps.stats, collated <- bind_rows(collated, imps.stats))
  if (successes>=10) {;
    break
  }
}


# Add a column containing consistent variable names
Boruta5k.results <- collated


### Produce graphs
## Condense by variable counts
imp.counts <- summarise(Boruta5k.results %>% group_by(Variable), 'Count' = n())
# Plot variable counts
imp.count.plot <- ggplot(imp.counts, aes(Variable, Count)) + 
  geom_col(fill='darkgreen') + theme_classic() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  scale_y_continuous(breaks=seq(0,150,30)) + labs(title='Mesothelioma primary cell lines CNV  - Boruta Frequency of Importance')


## Visualise median normalised hits
normHits.plot <- ggplot(Boruta5k.results, aes(Variable, normHits)) + 
  geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  labs(title='Mesothelioma primary cell lines CNV - Boruta Normalised Hits')

#plots

plot(imp.count.plot)

plot(normHits.plot)
