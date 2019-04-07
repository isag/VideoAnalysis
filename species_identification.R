#-------------------------------------------------#
#--- EXAMPLE SCRIPT FOR SPECIES IDENTIFICATION ---# (example to identify Pau and Col in mixture)
#-------------------------------------------------#


#--- Clear variables
rm(list = ls())


#--- LOAD LIBRARIES
#------------------
#install.packages("devtools",dependencies = T)
#library(devtools)
# install_github("efronhofer/bemovi",  ref="experimental")
library(bemovi)

#install.packages("e1071",dependencies = T)
library(e1071)   # for svm model



#--- Directory and file names (this is an example with separate folder for monocultures used to train the model (machine learning process))
#----------------------------
video.description.folder = "0_video_description/"
video.description.file = "video_description.txt"
merged.data.folder = "5_merged_data/"
to.data.mono = "~/Documents/Recherche/PROJECTS/EW10_BLUE-GREEN/2_DATA/Monoculture_test_20160502/"  # Change this path to your own folder!!!
to.data.mixt = paste(getwd(),"/",sep="")    # Here it assumes that your data are in the current working directory



#--- Parameters (those used to obtain the traits)
#--------------
fps = 25 #... Video frame rate (in frames per second)
nsv = 5 #... Number of seconds per videos
measured_volume = 34.4 #... Volume in µL, for Leica M205 C with 1.6 fold magnification, sample height 0.5 mm and Hamamatsu Orca Flash 4
pixel_to_scale = 4.05 #... Size of a pixel in µm, for Leica M205 C with 1.6 fold magnification, sample height 0.5 mm and Hamamatsu Orca Flash 4
filter_min_net_disp = 20 #... minimum net displacement (µm)
filter_min_duration = 0.5 #... duration (s)
filter_detection_freq = 0.1 #... detection frequency
filter_median_step_length = 3 #... median step length [µm]




#------------------------------------------------#
#--- MODEL TO CLASSIFY INDIVIDUALS BY SPECIES ---#
#------------------------------------------------#

#--- Create the training dataset with monocultures of Paramecium and Colpidium
#-----------------------------------------------------------------------------

#--- Load the merged data of monocultures for Pau and Col
load(paste0(to.data.mono, merged.data.folder, "Master.RData"))

#--- Filter data: minimum net displacement, their duration, the detection frequency and the median step length
trajectory.data.filtered = filter_data(trajectory.data, filter_min_net_disp, filter_min_duration, filter_detection_freq, filter_median_step_length)

#--- summarize trajectory data to individual-based data 
morph_mvt0 = summarize_trajectories(trajectory.data.filtered, write = T, to.data, calculate.median=F, merged.data.folder)

#--- select the monoculture data
morph_mvt = subset(morph_mvt0, Community %in% c("Pau","Col"))  

      # or if you have a column Composition with modes "monculture" and "mixture" for instance: 
      morph_mvt = subset(morph_mvt0, Composition == "monoculture")  

      # you may want to filter out particles with trait data, for instance microflagelates in Pau monocultures
      morph_mvt = cbind(subset(morph_mvt0, Community == "Pau" & mean_minor > 45),subset(morph_mvt0, Community == "Col"))

#--- Create training dataset with only monocultures and no NAs in data
training_data = morph_mvt[complete.cases(morph_mvt), ]


#--- SVM model to compare all traits of individuals and discriminate species (you can also use random forest if you have balanced species densities)
#---------------------------------------------------------------------------
#(my column "Community" of the video description file give the name of the species in monoculture P or C)

svm1 = svm(factor(Community) ~   mean_grey + sd_grey + mean_area + sd_area + mean_perimeter +  mean_turning + sd_turning +
             sd_perimeter + mean_major + sd_major + mean_minor + sd_minor + mean_ar + sd_ar + duration +
             max_net  + net_disp + net_speed + gross_disp + max_step + min_step + sd_step +
             sd_gross_speed + max_gross_speed + min_gross_speed ,
           data=training_data, probability=T,na.action=na.pass)

# Alternatively, you may want to use random forest (if you have balanced densities among species) because it allows out of box test of the model
library(randomForest)

RFmodel = randomForest(factor(Community) ~   mean_grey + sd_grey + mean_area + sd_area + mean_perimeter +  mean_turning + sd_turning +
                       sd_perimeter + mean_major + sd_major + mean_minor + sd_minor + mean_ar + sd_ar + duration +
                       max_net  + net_disp + net_speed + gross_disp + max_step + min_step + sd_step +
                       sd_gross_speed + max_gross_speed + min_gross_speed ,
                     data=training_dat, proximity=T)



#--- Add confusion matrix and error to check the error rate (as to be lower)
#----------------------------------------------------------
confusion.matrix = table(svm1$fitted,training_data$Community)
confusion.matrix.nd = confusion.matrix
diag(confusion.matrix.nd) = 0
svm1$confusion = cbind(confusion.matrix,class.error=rowSums(confusion.matrix.nd)/rowSums(confusion.matrix))
svm1$confusion



#-------------------------------------#
#--- IDENTIFICATION IN THE SAMPLES ---#
#-------------------------------------#

species.names = c("Pau","Col")

#--- Load the merged data where species are in mixture
load(paste0(to.data.mixt, merged.data.folder, "Master.RData"))

#--- Filter data: minimum net displacement, their duration, the detection frequency and the median step length
trajectory.data.filtered = filter_data(trajectory.data, filter_min_net_disp, filter_min_duration, filter_detection_freq, filter_median_step_length)

#--- Summarize trajectory data to individual-based data
morph_mvt = summarize_trajectories(trajectory.data.filtered, write = T, to.data, calculate.median=F, merged.data.folder)[,which(colnames(morph_mvt)!="Col_manual")]

#--- Remove the individuals with NAs: 
#---NAs in "turning" traits are for less than 2 frame trajectories, we assume they are neglectible (otherwise it prevents svm from predicting)
data.to.predict = morph_mvt[complete.cases(morph_mvt),]


#--- Predict the species identity of individuals with svm1 model (it assign the species with the highest probability based on the traits to each individual particle)
p.id = predict(svm1, data.to.predict, type="response")     
data.to.predict$predicted_species = as.character(p.id)


#--- Population level data
pop.data = summarize_populations(trajectory.data.filtered, morph_mvt, write=T, to.data, merged.data.folder, video.description.folder, video.description.file,total_frame = fps * nsv)




#--- Function to add a column for the density of each species from a morph_mvt object having a predicted_species column
#----------------------------------------------------------------------------------------------------------------------
#--- sample_output is the population level data (line = a sample)
#--- indiv_predicted is the data.frame with predicted species of individuals, their sample, and the number of frames on which each appears (N_frames)
#--- species_names gives the names of the species (as specified in the predicted_species column)
#--- total_frames is the total number of frames on each videos (frames per second x number of seconds)

species.density = function(sample_output,indiv_predicted,species_names,total_frames,mv = measured_volume)
{
  #... get the names of the samples where there are individuals
  samples = unique(indiv_predicted$file)
  
  #... create the matrix of species densities
  sp.dens = matrix(0,nrow(sample_output),length(species_names))
  colnames(sp.dens) = species_names
  
  for(i in 1:length(samples))
  {
    #... select the data for each sample
    indiv = subset(indiv_predicted,file == samples[i])
    
    #... get the species present in there
    spec = unique(indiv$predicted_species)
    for(j in 1:length(spec))
    {
      #... select the data for one species in the sample
      all.indiv.sp = subset(indiv,predicted_species == spec[j])
      
      #... calculate its density from the total number of frames on which each individuals of the species is present
      dens = sum(all.indiv.sp$N_frames)/total_frames/mv
      sp.dens[which(sample_output$file == as.character(samples[i])) ,which(species_names == spec[j])] = dens
    }
  }
  return(cbind(sample_output,sp.dens)) 
}



#--- Add the population density for each species
output = species.density(pop.data,data.to.predict,species.names,total_frames = fps*nsv, mv = measured_volume)


### For the case you have a high error rate in species assignation, because of overlapping trait distributions among species, 
#-- some more complex procedure have been developped in the supplementary information of the following paper:
#-- Harvey* E., Gounand* I., Fronhofer E.A. & Altermatt F. (2018) Proceedings of the Royal Society B, 285: 20182441; doi: 10.1098/rspb.2018.2441
#-- Improvement was gained by adding information on trait distribution in the model (skewedness, kurtosis, quartiles), and by combining models discriminating better small versus large species

