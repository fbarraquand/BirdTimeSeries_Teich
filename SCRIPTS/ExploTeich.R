###############################################################################################################################
### Frédéric Barraquand 16/05/2016 - Exploratory analysis Teich wetland bird (Le Teich Reserve, PNR Landes Gascogne) dataset
### Edited 01/07/2016 - Production of the "monthly photograph" = max of summed abondance over all LDs per month
#################################################################################################################################
### EDITED 27/04/2018 - Production of the "monthly photograph" with added of new data (2007-2008)
### EDITED 05/09/2019 - Removed hard-coded paths
#################################################################################################################################

rm(list=ls())
graphics.off()
library("lattice")

DIRECTORY_ORIGIN = "."
setwd(DIRECTORY_ORIGIN)
#setwd(paste(DIRECTORY_ORIGIN,"IN/",sep=""))

#DB<-read.delim(file="data_ROT20160324.txt",header=TRUE)#"Latin-1"
#DB<-read.csv(file="/home/caluome/Documents/DATA/DATA/le_Teich/data_ROT20160324.csv",header=TRUE,sep="\t")
DB<-read.csv(file="./IN/Initial_files/data_ROT20160324.csv",header=TRUE,sep="\t")
DBori=DB
#setwd(paste(DIRECTORY_ORIGIN,"OUT/photo/",sep=""))
##############################################
### Exploratory analyses and dataset cleaning
##############################################

#### How is the sampling effort distributed?
hist(DB$Annee, main = "Distribution of records, 1 per species") ### Toutes espèces confondues

unique(DB$Nom_latin) ### 300 espèces
## Quelques espèces faciles - "golden standard"

#hist(DB$Annee[DB$Nom_latin=="Anas strepera"]) 
#hist(DB$Annee[DB$Nom_latin=="Anas clypeata"]) 
### Around 2010 - large increase in sampling effort?

### Spatial variation in sampling? 
unique(DB$ID.Lieu.dit)
unique(DB$Lieu_dit)
### Many things...

### Let's see how many Lieu dits truly contain many observations
vec_lieu_dit=as.character(unique(DB$Lieu_dit))
vec_lieu_dit
n_obs=rep(0,length(vec_lieu_dit))
for (i in 1:length(vec_lieu_dit)){
  n_obs[i]=sum(DB$Lieu_dit==vec_lieu_dit[i])
}

#barplot(log10(n_obs),ylab="log10(n_obs at Lieu dit)",xlab="Lieu dit")
### histogram(log(n_obs)) should work as well but somehow doesn't

############################################################################################################
### Restrict the dataset based on Lieu-dits with many obs (to do spatially after, i.e. within Le Teich)
############### We don't do this anymore - see below
#vec_lieu_dit_frequent=unique(DB$Lieu_dit)[n_obs>100] #20 Lieu dits avec >100 obs
#n_obs[n_obs>100]
### Restricting the dataset to these 
#DB=subset(DB,DB$Lieu_dit %in% vec_lieu_dit_frequent)
#unique(DB$Lieu_dit)
#hist(DB$Annee) ### No change
### Have we suppressed stuff from the 1970s? Perhaps to check later. Now everything is in Le Teich. 
#unique(DB$Nom_latin) #285 species now
####################################################################################################
### Lieu dit vs year
#table_sampling_lieux_dits_Teich=table(as.character(DB$Lieu_dit),DB$Annee)
#write.csv2(table_sampling_lieux_dits_Teich,file = "Table_sampling_lieux_dits_Teich.csv")
#####################################################################################################

########################################### Info about the dataset temporal/spatial structure ########################
# Before 2006, everything in "USN00-Réserve ornithologique (générique)"
# Except in 2003 where most obs. are in "USN01-Artigues-Réserve ornithologique" (error?) // Ask Claude. 
# >=2007 everything is "regionalized". So we have a non-spatial dataset til 2007 and spatial after. 

### Question : do the "USN00-Réserve ornithologique (générique)" sum up all the other USN observations?
### Answer : yes (well, it should - it sums all the "protocoled" obs)
########################################################################################################################

###########################################################################################
### Complemetely preliminary exploratory analyses below. 
### Kept as comment in case we need to borrow code from there. 

## Let's do something crazy - brute force trends and community ecology descriptors
## The data is summed before 2007 - so any comparison post 2007 is crap for most species

# ### Number of species over time (we can see later how that evolves with the sampling effort)
# ### Combine this with a plot of the effort
# pdf("Sampling_and_numberOfSpeciesTeich.pdf",width=8,height=4)
# par(mfrow=c(1,2))
# hist(DB$Annee, main = "Distribution of records, 1 record/species",xlab="Year") ### Toutes espèces confondues
# abline(v=2007,lwd=2,col="red")
# 
# total_species_list=unique(DB$Nom_latin)
# n_species=rep(0,max(DB$Annee)-min(DB$Annee)+1)
# for (year in min(DB$Annee):max(DB$Annee)){
#   k=year-min(DB$Annee)+1
#   n_species[k]=length(unique( as.character(DB$Nom_latin[DB$Annee==year]) ))
# }
# plot(min(DB$Annee):max(DB$Annee),n_species,ylab="Number of species",xlab="Year", main = "~ 300 species total")
# abline(v=2007,lwd=2,col="red")
# ### Probably strongly correlated to sampling effort // increase in 2010 (or 2007?)
# dev.off()
# ### The abondances, however, averaged over individuals, should not be affected by sampling effort (in theory)
# 
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas strepera"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas strepera"]) # perhaps a real trend?
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas clypeata"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas clypeata"])
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas platyrhynchos"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas platyrhynchos"]) #colvert
# ### I hope they weren't summing observations in the 1990s and we just don't know about it...
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas platyrhynchos"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas platyrhynchos"])
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas penelope"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas penelope"])
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas crecca"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas crecca"])
# ### Perhaps there are just some changes in the duck community - or the duck numbers are not summed since 2007 over the various areas...
# 
# ### If this is problematic, we can sort the dataset in two parts, pre- and post-2007 for now. 
# 
# DB_lieux-dits_frequent=DB ### store
# 
# ############## Zoom on some years
# 
# ### Let's zoom on 2000-2005
# DB=subset(DB,DB$Annee>=2000 & DB$Annee<=2005)
# 
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas strepera"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas strepera"]) # perhaps a real trend?
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas clypeata"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas clypeata"])
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas platyrhynchos"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas platyrhynchos"]) #colvert
# ### I hope they weren't summing observations in the 1990s and we just don't know about it...
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas penelope"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas penelope"])
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas crecca"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas crecca"])
# 
# ### Let's zoom on 2010-2015
# DB=subset(DBori,DBori$Annee>=2010 & DBori$Annee<=2015)
# 
# 
# 
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas strepera"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas strepera"]) # perhaps a real trend?
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas clypeata"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas clypeata"])
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas platyrhynchos"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas platyrhynchos"]) #colvert
# ### I hope they weren't summing observations in the 1990s and we just don't know about it...
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas penelope"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas penelope"])
# 
# plot(as.Date(DB$Date[DB$Nom_latin=="Anas crecca"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas crecca"])
# #lines(as.Date(DB$Date[DB$Nom_latin=="Anas crecca"],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin=="Anas crecca"])
# ########################################################################################################################################
# #library(Hmisc)
# 
# ### Looking at all species and all years
# 
# DB=DB_lieux-dits_frequent
# vec_species_names=as.character(unique(DB$Nom_latin))
# vec_species_french=as.character(unique(DB$Nom_espece))
# 
# ### First cleaning all non-abundant species 
# n_obs=rep(0,length(vec_species_names))
# for (i in 1:length(vec_species_names)){
#   n_obs[i]=sum(DB$Nom_latin==vec_species_names[i])
# }
# n_obs
# n_obs[n_obs>75]
# vec_species_names_frequent=vec_species_names[n_obs>75]
# vec_species_french_frequent=vec_species_french[n_obs>75]
# DB=subset(DB,DB$Nom_latin %in% vec_species_names_frequent)
# 
# ### When are other locations counted?
# min_date_lieu_dit=as.Date(rep(NA,length(vec_lieu_dit_frequent)))
# for (i in 1:length(vec_lieu_dit_frequent)){
#   min_date_lieu_dit[i]=(as.Date(DB$Date[DB$Lieu_dit==vec_lieu_dit_frequent[i]],format = "%d.%m.%Y"))[1]
#   }
# as.Date(min_date_lieu_dit)
# 
# pdf("TrendsAbundance_Teich_raw.pdf",width=8,height=12)
# par(mfrow=c(5,1))#cex=1.5
# for (i in 1:length(vec_species_names_frequent)){
#   plot(as.Date(DB$Date[DB$Nom_latin==vec_species_names_frequent[i]],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin==vec_species_names_frequent[i]],ylab = c(vec_species_french_frequent[i],vec_species_names_frequent[i]),xlab="Date")
#   #minor.tick(nx=10, tick.ratio=0.5)
#   #axis(tck=-0.015,at=c(seq(from=min(as.Date(DB$Date),format = "%d.%m.%Y"),to=,by=100)))
#   abline(v = as.Date("01.01.2007",format = "%d.%m.%Y"),col="red",lwd=2)
#   lines(as.Date(DB$Date[DB$Nom_latin==vec_species_names_frequent[i]],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin==vec_species_names_frequent[i]])
#   }
# dev.off()
# 
# ### Are species that increase in the end "localised" - it would explain, if they are found in only one place, why they cannot decrease through subcounting?
# i=50
# plot(as.Date(DB$Date[DB$Nom_latin==vec_species_names_frequent[i]],format = "%d.%m.%Y"),DB$Nombre[DB$Nom_latin==vec_species_names_frequent[i]],ylab = c(vec_species_french_frequent[i],vec_species_names_frequent[i]),xlab="Date")
# 
# table(as.character(DB$Lieu_dit[DB$Nom_latin==vec_species_names_frequent[i]]),DB$Annee[DB$Nom_latin==vec_species_names_frequent[i]])
# 
# table(as.Date(DB$Date[DB$Nom_latin==vec_species_names_frequent[i]&DB$Annee==2010],format = "%d.%m.%Y"),as.character(DB$ID.Lieu.dit[DB$Nom_latin==vec_species_names_frequent[i]&DB$Annee==2010]))
# unique(DB$Lieu_dit[DB$ID.Lieu.dit==6663])
# 
# DB_Recur2010=subset(DB,DB$Nom_latin==vec_species_names_frequent[i]&DB$Annee==2010)
# 
# DB_Recur2010$Lieu_dit[DB_Recur2010$Nombre>100] # quasi-tous sur vasière spatule. Probablement nous aurons à faire beaucoup de cas par cas comme celui-là. 
# 
# 
# DB_Strepera2010=subset(DB,DB$Nom_latin==vec_species_names_frequent[3]&DB$Annee==2010)
# DB_Strepera2010$Lieu_dit[DB_Strepera2010$Nombre>10] 
# 
# ### Observer effect?
# unique(DB$Nom) 
# unique(DB$Personne.morale) 
# unique(DB$Prénom.transmetteur) 
# unique(DB$Nom.transmetteur) 
# ### Quite obviously we won't find any. 
# #################################################### End of crazy exploratory stuff ########################################################

########################################################################################################################################
########### 17/06/2015 - Re-analysis after remarks by Claude. 
########################################################################################################################################

# Restrict database to USN. (remove Fleury and Pointe de l'Eyre, add USN11)
DB=DBori
head(DB)

vec_lieu_dit
DB=DB[grep("USN",DB$Lieu_dit),] ## select all lieux-dits within the reserve
table_sampling_lieux_dits_Teich=table(as.character(DB$Lieu_dit),DB$Annee)
table_sampling_lieux_dits_Teich
write.csv2(table_sampling_lieux_dits_Teich,file = "./Preliminary_analyses/OUT/photo/Table_sampling_lieux_dits_withinTeich.csv")

# Clean the date 
DB$Date=as.Date(DB$Date,format = "%d.%m.%Y")

# Sub-select fields to avoid crazy big database 
DB=subset(DB, select=c("Ref", "Nom_espece","Nom_latin","Date","Jour","Mois","Annee","Jour_de_l_année","ID.Lieu.dit","Lieu_dit","Nombre","Commentaire","Remarque_privee"))

## See how the data post-2007 is structured between dates and lieu-dit
#table_sampling2010_lieux_dits_Teich=table(as.character(DB$Lieu_dit[DB$Annee==2010]),as.Date(DB$Date[DB$Annee==2010],format = "%d.%m.%Y"))
#write.csv2(table_sampling2010_lieux_dits_Teich,file = "Table_sampling2010_lieux_dits_byDate.csv")

table_sampling2003_lieux_dits_Teich=table(as.character(DB$Lieu_dit[DB$Annee==2003]),DB$Date[DB$Annee==2003])
#write.csv2(table_sampling2003_lieux_dits_Teich,file = "Table_sampling2003_lieux_dits_byDate.csv")

## Look at 2002 to compare
table_sampling2002_lieux_dits_Teich=table(as.character(DB$Lieu_dit[DB$Annee==2002]),DB$Date[DB$Annee==2002])
#write.csv2(table_sampling2002_lieux_dits_Teich,file = "Table_sampling2002_lieux_dits_byDate.csv")


##################################
# Create "monthly photograph"
# Before 2007, use the 15 of the Month
CreateMonthly_photograph=function(DB){ 
  DB$Protocol=0
  DB$number_of_LDs_surveyed=0
  ### Show data of the fifth of the month. 
  head(DB[(DB$Annee<2007)&(as.POSIXlt(DB$Date)$mday==15),])
  ### Pas besoin de truc si compliqué...
  DB[(DB$Annee<2007)&(DB$Jour==15),]$Protocol=1
  ### All data entered the 15th of the Month is "protocoled"
  
  
  # ############################# New code to define "montly photograph after 2007" ###################
  # ### After 2006, we need to find the dates when at least 10 to 12 LD or counted
  ##### Then we sum all the LDs, per species, per year and per month 
  # before 2007, the insert date is always the 15th of each month (see above prtocol 1), we don't considere the other dates.
  
  vec_species_names=as.character(unique(DB$Nom_latin))
  vec_species_french=as.character(unique(DB$Nom_espece))
  #setwd(paste(DIRECTORY_ORIGIN,"IN/",sep=""))
  #Now subset the max in the month
  minyear=2007
  maxyear=max(unique(DB$Annee))
  for (year in minyear:maxyear) # Loop on years
    {
      for (month in 1:12) # Loop on month
        {
          DB2=subset(DB,(DB$Annee==year)&(DB$Mois==month)) # Select a dataset for the month and year, for all species
          head(DB2)
          ### dates_vector=unique(DB2$Date) # This has drawbacks...
          dates_vector=DB2$Date
          length(dates_vector) #length of DB2
          if (nrow(DB2)>0) # If there is data collected this month
          {
          number_of_LDs_surveyed=rep(0,nrow(DB2)) # initialize vector
          for (i in 1:length(dates_vector)){
            number_of_LDs_surveyed[i]=length(unique(DB2$ID.Lieu.dit[DB2$Date==dates_vector[i]])) #count how many unique LDs
          }
          number_of_LDs_surveyed # print vector
          DB2$number_of_LDs_surveyed = number_of_LDs_surveyed ## adding column on the number of LDs surveyed to dataframe
          protocoled_survey_dates = DB2$Date[DB2$number_of_LDs_surveyed>10] ## more than 10 LDs surveyed means the data is "protocoled". 
  	protocoled_survey_dates=unique(protocoled_survey_dates) # keep the dates for which 
          protocoled_survey_jour_annee = unique(DB2$Jour_de_l_année[DB2$number_of_LDs_surveyed>10])
          ## seems to be problems to find >10 LDs for 2007-2008
      
          ### Species observed that month and year
          species_observed = unique(DB2$Nom_latin)
          ## 34 Limosa limosa
          for (i in 1:length(vec_species_names))
          {
        
            # If the species is in the set for that month and year // I don't have included a condition 
            if (vec_species_names[i] %in% species_observed)
             {
                DB3=subset(DB2,DB2$Nom_latin==vec_species_names[i])
                # First find the protocoled dates for that species (if they can be found, otherwise use all of the available ones)
                # DB3$Date %in% protocoled_survey_dates
                # Wait: this will induce bias if some rare species are missed at the time of counting...
                # Better to do that above the species loop then
                if (length(protocoled_survey_dates)>0){ #if there's at least one protocoled survey
  	   sum_over_LDs=rep(0,length(protocoled_survey_dates))
               
                for (k in 1:length(protocoled_survey_dates))
                  {
                   sum_over_LDs[k] = sum(DB3$Nombre[DB3$Date==protocoled_survey_dates[k]])
                  }
                sum_over_LDs
                
                monthly_max = max(sum_over_LDs)
                index_max = which.max(sum_over_LDs)
            
                ### Now enter properly that info into the database
        
                new_row = DB3[1,] ## Initializing, no problem as values of interest will change
                new_row$Nombre = monthly_max
                new_row$Ref = 000000 # to show it has been modified, as it should have. 
                new_row$Lieu_dit = "USN00-Réserve ornithologique (générique)"
                new_row$Date = protocoled_survey_dates[index_max] # change the date to that when the max of sums was found
                new_row$Jour = as.POSIXlt(new_row$Date)$mday
                new_row$Jour_de_l_année = protocoled_survey_jour_annee[index_max]
                #new_row$Jour_de_l_année = unique(DB3$Jour_de_l_année[DB3$Date==protocoled_survey_dates[index_max]])
                #creates problems whenever there are positive observations out of the protocoled dates and zero at protocoled dates
                new_row$Mois = month
                new_row$Annee = year
                new_row$Commentaire = "Max of the month of summed abondance"
                new_row$Protocol=1
                new_row$number_of_LDs_surveyed=unique(DB2$number_of_LDs_surveyed[DB2$Date==new_row$Date]) #The fact that I'm looking at DB2 and not DB3 is because we can have no record for a species at a protocoled date, so sum will be 0, but the line unique(...) will return nothing #CP
  	    }  else  {  # Use all the dates in DB3
                protocoled_survey_dates_local=unique(DB3$Date) #CP : hereafter, protocoled_survey_dates was changed to protocoled_survey_dates_local and protocoled_survey_jour_annee -> protocoled_survey_jour_annee_local
                protocoled_survey_jour_annee_local=unique(DB3$Jour_de_l_année)
                
                sum_over_LDs=rep(0,length(protocoled_survey_dates_local))
                for (k in 1:length(protocoled_survey_dates_local))
                {
                  sum_over_LDs[k] = sum(DB3$Nombre[DB3$Date==protocoled_survey_dates_local[k]])
                }
                sum_over_LDs
                monthly_max = max(sum_over_LDs)
                index_max = which.max(sum_over_LDs)
            
                ### Now enter properly that info into the database
                new_row = DB3[1,] ## Initializing, no problem as values of interest will change
                new_row$Nombre = monthly_max
                new_row$Ref = 000001 # to show it has been modified, but we could not define a day for which all LDs where counted
                new_row$Lieu_dit = "USN00-Réserve ornithologique (générique)"
                new_row$Date = protocoled_survey_dates_local[index_max] # change the date to that when the max of sums was found
                new_row$Jour = as.POSIXlt(new_row$Date,format = "%d.%m.%Y")$mday
                new_row$Jour_de_l_année = protocoled_survey_jour_annee_local[index_max]
                #new_row$Jour_de_l_année = unique(DB3$Jour_de_l_année[DB3$Date==protocoled_survey_dates[index_max]])
                #creates problems whenever there are positive observations out of the protocoled dates and zero at protocoled dates
                new_row$Mois = month
                new_row$Annee = year
                new_row$Commentaire = "Max of the month of summed abondance"
                new_row$Protocol=2 # Not quite the perfect protocol
                new_row$number_of_LDs_surveyed=unique(DB2$number_of_LDs_surveyed[DB2$Date==new_row$Date]) #CP
              }
            # Again I am not sure this cannot induce some bias and it should be checked thoroughly later. 
          
            # Now we have to add the created row to the dataframe
            DB=rbind(DB,new_row)
            # Don't add stuff to DB3 or DB2 to avoid mistakes
            } # end of condition stating whether the species is observed that month or not
          } #end of loop on species
        } #end of condition on whether there is data this month. 
          Sys.sleep(0.1)
          print(month)
      }   # end of loop on months
    Sys.sleep(0.1)
    print(year)
    } #end of loop on years
  return (DB)
}
DB = CreateMonthly_photograph(DB)
write.csv(DB,file = "./Preliminary_analyses/IN/DBWithMonthlyPhotoTeich_protocol2corrected.csv")

### Plot the stuff now
vec_species_names=as.character(unique(DB$Nom_latin))
vec_species_french=as.character(unique(DB$Nom_espece))
### First cleaning all non-abundant species 
n_obs=rep(0,length(vec_species_names))
for (i in 1:length(vec_species_names)){
  n_obs[i]=sum(DB$Nom_latin==vec_species_names[i])
}
n_obs
n_obs[n_obs>75]
vec_species_names_frequent=vec_species_names[n_obs>75]
vec_species_french_frequent=vec_species_french[n_obs>75]

pdf("./Preliminary_analyses/OUT/photo/TrendsAbundance_Teich_monthly_photograph.pdf",width=8,height=12)
par(mfrow=c(5,1))#cex=1.5
for (i in 1:length(vec_species_names_frequent)){
  plot(DB$Date[(DB$Nom_latin==vec_species_names_frequent[i])&((DB$Protocol==1)|(DB$Protocol==2))],DB$Nombre[(DB$Nom_latin==vec_species_names_frequent[i])&((DB$Protocol==1)|(DB$Protocol==2))],ylab = c(vec_species_french_frequent[i],vec_species_names_frequent[i]),xlab="Date")
  #minor.tick(nx=10, tick.ratio=0.5)
  #axis(tck=-0.015,at=c(seq(from=min(as.Date(DB$Date),format = "%d.%m.%Y"),to=,by=100)))
  abline(v = as.Date("01.01.2007",format = "%d.%m.%Y"),col="red",lwd=2)
  lines(as.Date(DB$Date[(DB$Nom_latin==vec_species_names_frequent[i])&((DB$Protocol==1)|(DB$Protocol==2))],format = "%d.%m.%Y"),DB$Nombre[(DB$Nom_latin==vec_species_names_frequent[i])&((DB$Protocol==1)|(DB$Protocol==2))])
  points(DB$Date[(DB$Nom_latin==vec_species_names_frequent[i])&(DB$Protocol==2)],DB$Nombre[(DB$Nom_latin==vec_species_names_frequent[i])&(DB$Protocol==2)],col="blue")
}
dev.off()

### Just to check we get the same thing

pdf("./Preliminary_analyses/OUT/photo/TrendsAbundance_Teich_monthly_photograph.pdf",width=8,height=12)
par(mfrow=c(5,1))#cex=1.5
for (i in 1:length(vec_species_names_frequent)){
  plot(DB$Date[(DB$Nom_latin==vec_species_names_frequent[i])&((DB$Protocol==1)|(DB$Protocol==2))&(DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")],DB$Nombre[(DB$Nom_latin==vec_species_names_frequent[i])&((DB$Protocol==1)|(DB$Protocol==2))&(DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")],ylab = c(vec_species_french_frequent[i],vec_species_names_frequent[i]),xlab="Date")
  #minor.tick(nx=10, tick.ratio=0.5)
  #axis(tck=-0.015,at=c(seq(from=min(as.Date(DB$Date),format = "%d.%m.%Y"),to=,by=100)))
  abline(v = as.Date("01.01.2007",format = "%d.%m.%Y"),col="red",lwd=2)
  lines(as.Date(DB$Date[(DB$Nom_latin==vec_species_names_frequent[i])&((DB$Protocol==1)|(DB$Protocol==2))&(DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")],format = "%d.%m.%Y"),DB$Nombre[(DB$Nom_latin==vec_species_names_frequent[i])&((DB$Protocol==1)|(DB$Protocol==2))&(DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")])
  points(DB$Date[(DB$Nom_latin==vec_species_names_frequent[i])&(DB$Protocol==2)&(DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")],DB$Nombre[(DB$Nom_latin==vec_species_names_frequent[i])&(DB$Protocol==2)&(DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")],col="blue")
}
dev.off()

#### Start checking ##################
### Just try one species that increases and one that decreases
### Start with Limosa limosa (increases)

DB[(DB$Nom_latin=="Limosa limosa")&(DB$Annee==2008),]
DB[(DB$Nom_latin=="Limosa limosa")&(DB$Annee==2008)&(DB$Mois==1),]
# Be very careful of syntax when doing this - I managed to change the DB here by having = instead of == ... (damn R!)
#Why do we have Protocol=1? That does not sound right... But that's actually OK -> Limosa was counted once at that date but other species were counted in other LDs



  
############# Other suff to use ########################

### Check year 2007
## See how the data post-2007 is structured between dates and lieu-dit
table_sampling2007_lieux_dits_Teich=table(as.character(DB$Lieu_dit[DB$Annee==2007]),as.Date(DB$Date[DB$Annee==2007],format = "%d.%m.%Y"))
write.csv2(table_sampling2007_lieux_dits_Teich,file = "./Preliminary_analyses/OUT/photo/Table_sampling2007_lieux_dits_byDate.csv")

### Perhaps smarter ways to do this here
### http://stackoverflow.com/questions/3505701/r-grouping-functions-sapply-vs-lapply-vs-apply-vs-tapply-vs-by-vs-aggrega

#setwd(paste(DIRECTORY_ORIGIN,"IN/",sep=""))
################ New code using only the max of the summed abondance per month ############
# it's just for check the difference with the "montly photograph after 2007"
minyear=2007
maxyear=max(unique(DB$Annee))
for (year in minyear:maxyear) # Loop on years
{
  for (month in 1:12) # Loop on month
  {
    DB2=subset(DB,(DB$Annee==year)&(DB$Mois==month)) # Select a dataset for the month and year, for all species
    head(DB2)
    ### dates_vector=unique(DB2$Date) # This has drawbacks...
    dates_vector=DB2$Date
    length(dates_vector) #length of DB2
    if (nrow(DB2)>0) # If there is data collected this month
    {
      number_of_LDs_surveyed=rep(0,nrow(DB2)) # initialize vector
      for (i in 1:length(dates_vector)){
        number_of_LDs_surveyed[i]=length(unique(DB2$ID.Lieu.dit[DB2$Date==dates_vector[i]])) #count how many unique LDs
      }
      number_of_LDs_surveyed # print vector
      DB2$number_of_LDs_surveyed = number_of_LDs_surveyed ## adding column on the number of LDs surveyed to dataframe
      protocoled_survey_dates = DB2$Date[DB2$number_of_LDs_surveyed>10] ## more than 10 LDs surveyed means the data is "protocoled". 
      protocoled_survey_dates=unique(protocoled_survey_dates) # keep the dates for which 
      protocoled_survey_jour_annee = unique(DB2$Jour_de_l_année[DB2$number_of_LDs_surveyed>10])
      ## seems to be problems to find >10 LDs for 2007-2008
      
      ### Species observed that month and year
      species_observed = unique(DB2$Nom_latin)
      
      ## 34 Limosa limosa
      for (i in 1:length(vec_species_names))
      {
        # Old code - completely fucking wrong
        #DB2$Nombre[(DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")&(DB$Nom_latin==vec_species_names[i])]=
        #DB$Nombre[(DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")&(DB$Annee==year)&(DB$Mois==month)&(DB$Nom_latin==vec_species_names[i])]=max(DB2$Nombre[(DB2$Protocol==1)&(DB2$Nom_latin==)])### to implement - sum over all lieux dit
        
        # If the species is in the set for that month and year // I don't have included a condition 
        if (vec_species_names[i] %in% species_observed)
        {
          DB3=subset(DB2,DB2$Nom_latin==vec_species_names[i])
          # First find the protocoled dates for that species (if they can be found, otherwise use all of the available ones)
          # DB3$Date %in% protocoled_survey_dates
          # Wait: this will induce bias if some rare species are missed at the time of counting...
          # Better to do that above the species loop then
          
          if (0>1){ #to force ignorance
            sum_over_LDs=rep(0,length(protocoled_survey_dates))
            
            for (k in 1:length(protocoled_survey_dates))
            {
              sum_over_LDs[k] = sum(DB3$Nombre[DB3$Date==protocoled_survey_dates[k]])
            }
            sum_over_LDs
            
            monthly_max = max(sum_over_LDs)
            index_max = which.max(sum_over_LDs)
            
            ### Now enter properly that info into the database
            
            new_row = DB3[1,] ## Initializing, no problem as values of interest will change
            new_row$Nombre = monthly_max
            new_row$Ref = 000000 # to show it has been modified, as it should have. 
            new_row$Lieu_dit = "USN00-Réserve ornithologique (générique)"
            new_row$Date = protocoled_survey_dates[index_max] # change the date to that when the max of sums was found
            new_row$Jour = as.POSIXlt(new_row$Date,format="%d.%d.%Y")$mday
            new_row$Jour_de_l_année = protocoled_survey_jour_annee[index_max]
            #new_row$Jour_de_l_année = unique(DB3$Jour_de_l_année[DB3$Date==protocoled_survey_dates[index_max]])
            #creates problems whenever there are positive observations out of the protocoled dates and zero at protocoled dates
            new_row$Mois = month
            new_row$Annee = year
            new_row$Commentaire = "Max of the month of summed abondance"
            new_row$Protocol=1
          }  else  {  # Use all the dates in DB3
            protocoled_survey_dates=unique(DB3$Date)
            protocoled_survey_jour_annee=unique(DB3$Jour_de_l_année)
            
            sum_over_LDs=rep(0,length(protocoled_survey_dates))
            for (k in 1:length(protocoled_survey_dates))
            {
              sum_over_LDs[k] = sum(DB3$Nombre[DB3$Date==protocoled_survey_dates[k]])
            }
            sum_over_LDs
            monthly_max = max(sum_over_LDs)
            index_max = which.max(sum_over_LDs)
            
            ### Now enter properly that info into the database
            new_row = DB3[1,] ## Initializing, no problem as values of interest will change
            new_row$Nombre = monthly_max
            new_row$Ref = 000001 # to show it has been modified, but we could not define a day for which all LDs where counted
            new_row$Lieu_dit = "USN00-Réserve ornithologique (générique)"
            new_row$Date = protocoled_survey_dates[index_max] # change the date to that when the max of sums was found
            new_row$Jour = as.POSIXlt(new_row$Date)$mday
            new_row$Jour_de_l_année = protocoled_survey_jour_annee[index_max]
            #new_row$Jour_de_l_année = unique(DB3$Jour_de_l_année[DB3$Date==protocoled_survey_dates[index_max]])
            #creates problems whenever there are positive observations out of the protocoled dates and zero at protocoled dates
            new_row$Mois = month
            new_row$Annee = year
            new_row$Commentaire = "Max of the month of summed abondance"
            new_row$Protocol=2 # Not quite the perfect protocol
          }
          # Again I am not sure this cannot induce some bias and it should be checked thoroughly later. 
          
          # Now we have to add the created row to the dataframe
          DB=rbind(DB,new_row)
          # Don't add stuff to DB3 or DB2 to avoid mistakes
        } # end of condition stating whether the species is observed that month or not
      } #end of loop on species
    } #end of condition on whether there is data this month. 
    Sys.sleep(0.1)
    print(month)
  }   # end of loop on months
  Sys.sleep(0.1)
  print(year)
} #end of loop on years

write.csv(DB,file = "./Preliminary_analyses/IN/DBWithMonthlyPhotoTeich_maxsum.csv")
#setwd(paste(DIRECTORY_ORIGIN,"OUT/photo/",sep=""))

### Plot the stuff now


pdf("./Preliminary_analyses/OUT/photo/TrendsAbundance_Teich_max_sum_month_post2007.pdf",width=8,height=12)
par(mfrow=c(5,1))#cex=1.5
for (i in 1:length(vec_species_names_frequent)){
  x=DB$Date[(DB$Nom_latin==vec_species_names_frequent[i])&(DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")&((DB$Protocol==1)|(DB$Protocol==2))]
  y=DB$Nombre[(DB$Nom_latin==vec_species_names_frequent[i])&(DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")&((DB$Protocol==1)|(DB$Protocol==2))]
  plot(x[order(x)],y[order(x)],ylab = c(vec_species_french_frequent[i],vec_species_names_frequent[i]),xlab="Date")
  #minor.tick(nx=10, tick.ratio=0.5)
  #axis(tck=-0.015,at=c(seq(from=min(as.Date(DB$Date),format = "%d.%m.%Y"),to=,by=100)))
  abline(v = as.Date("01.01.2007",format = "%d.%m.%Y"),col="red",lwd=2)
  lines(x[order(x)],y[order(x)])
  points(DB$Date[(DB$Nom_latin==vec_species_names_frequent[i])&(DB$Protocol==2)&(DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")],DB$Nombre[(DB$Nom_latin==vec_species_names_frequent[i])&(DB$Protocol==2)&(DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")],col="blue")
}
dev.off()


# ############################# ADD NEW DATA IN "montly photograph after 2007" ################### ###################
## ADD 27/04/2018 by Christelle A.

# import file
#DB1<-read.csv(file="./IN/data_ROT20160324.csv",header=TRUE,sep="\t",dec=".")
DB1=DBori
DB1= DB1[,1:76]
DB1=DB1[grep("USN",DB1$Lieu_dit),] ## select all lieux-dits within the reserve
## new data to add
DBadd<-read.csv(file="./IN/Initial_files/export_30062016_095620.csv",header=TRUE,sep=",",dec=".")
colnames(DBadd)=colnames(DB1)#to have the same name of columns, it's the same information, I checked.
# solves the date problems
DB1$Date =as.Date(as.character(DB1$Date),format="%d.%m.%Y")
DBadd$Date = as.Date(paste(as.character(DBadd$Jour),"-",as.character(DBadd$Mois),"-",as.character(DBadd$Annee),sep=""),format="%d-%m-%Y")

DB1=DB1[grep("USN",DB1$Lieu_dit),] ## select all lieux-dits within the reserve
DB1_temp=DB1


DB=DB1
## small function, for find and delete duplicate
doublon =c()
for (i in 1:dim(DBadd)[1]){ # find indice of doublon, check if data in DBadd exists in DB
  temp = DBadd[i,]
  if(! dim(subset(DB,DB$Annee==temp$Annee& DB$Mois == temp$Mois
                  & DB$Jour == temp$Jour & DB$Nombre == temp$Nombre 
                  &  as.character(DB$Nom_latin) == as.character(temp$Nom_latin)))[1]==0){
    doublon=c(doublon,i)
  }
}
DBadd=DBadd[-doublon,] # remove all the duplicate

DB=rbind(DB1,DBadd) # merges the two datasets
DB=subset(DB, select=c("Ref", "Nom_espece","Nom_latin","Date","Jour","Mois","Annee","Jour_de_l_année","ID.Lieu.dit","Lieu_dit","Nombre","Commentaire","Remarque_privee")) # Select the columns of interest

label = "DBWithMonthlyPhotoTeich_completed" # new name of the database
DB = CreateMonthly_photograph(DB) #creation of the DB, with all the data
# write the DB
write.csv(DB,file = paste("./IN/",label,".csv",sep=""))



## For check the database, I aligned the old database and the new.
## No difference, before 2006 and after 2008.
## and between 2006 and 2008 they are for some species, a lot of new points (good !)
### pour faire le graphique.
#DB<-read.csv(file="/home/caluome/Documents/LabExCOTE/DATA/DATA/le_Teich/DBWithMonthlyPhotoTeich_completed.csv",header=TRUE,sep=",",dec=".")
DB1_2ssD = subset(DB,(DB$Annee>2005 & DB$Annee<2010) & (DB$Protocol==1 | DB$Protocol==2) & DB$Lieu_dit=="USN00-Réserve ornithologique (générique)")
DB1_2ssD$Date =as.Date(as.character(DB1_2ssD$Date),format="%Y-%m-%d")
#DB1<-read.csv(file="/home/caluome/Documents/LabExCOTE/DATA/Graphes/Axe1/TRANSFERT_LIMICOLES/IN/DBWithMonthlyPhotoTeich.csv",header=TRUE,sep=",",dec=".")
DB1<-read.csv(file="./Preliminary_analyses/IN/DBWithMonthlyPhotoTeich_protocol2corrected.csv",header=TRUE,sep=",",dec=".")


DB1 =  subset(DB1,(DB1$Annee>2005 & DB1$Annee<2010) & (DB1$Protocol==1 | DB1$Protocol==2) & DB1$Lieu_dit=="USN00-Réserve ornithologique (générique)")
DB1$Date =as.Date(as.character(DB1$Date),format="%Y-%m-%d")

pdf("./OUT/photo/comparison_MonthlyPhoto_with_data_addition.pdf",width=8,height=4)
AllSpeciesL = as.character(unique(DB1_2ssD$Nom_latin))

for(i in 1:length(AllSpeciesL)){
  print (AllSpeciesL[i])
  par(mar=c(4,4.5,3,2.5))
  par(oma = c(4.2, 0.5, 0.5, 0.5))
  data = subset(DB1_2ssD,as.character(DB1_2ssD$Nom_latin)==AllSpeciesL[i])
  plot(data$Date,data$Nombre,col="violet",type="o",cex=2,pch=16,
       xlim = c(as.Date("2006-01-01",format="%Y-%m-%d"),
       as.Date("2010-04-01",format="%Y-%m-%d")),ylim=c(0,max(data$Nombre,na.rm = TRUE)),
       xlab = "date",ylab = "Comptage", main = paste("Visualization of the years 2006-2008 for ",AllSpeciesL[i],sep=""))
  data = subset(DB1,as.character(DB1$Nom_latin)==AllSpeciesL[i])
  if(dim(data)[1]!=0){
    points(data$Date,data$Nombre,col="blue",cex=0.7,pch=16,type="o")
    for(d in 1:length(data$Date)){
      abline(v=data$Date[d],lwd=1,lty= 4, col="black")
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  
  legend("bottom",c("First MonthlyPhoto","Second MonthlyPhoto (adding data for 2006-2008), eliminating perfect duplicates."), col=c("blue","violet"),pch=c(16,16), lty=c(1,1), xpd = TRUE,
         horiz=FALSE,
         inset = c(0,0),bty = "n",cex=1,pt.cex=1.5,pt.lwd=1.5)
}
dev.off()
