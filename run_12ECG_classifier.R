library(randomForest)
library(signal)
library(e1071)

run_12ECG_classifier<-function(data,hea_data,classes,model){
  
  features<-c()
  features<-get_12ECG_features(data,hea_data)
  
  classes_model<-model[["classes"]]
  
  predicted_class<-predict(model,features,type="response")
  predicted_prob<-predict(model,features,type="prob")
  
  class_vector<-rep(0,length(classes))
  class_vector[grep(predicted_class,classes)]=1
  
  prob_vector<-rep(0,length(classes))
  
  for (i in 1:length(classes)) {
    pos<-grep(classes[i],classes_model)
    prob_vector[i]<-predicted_prob[pos]
  }
  
  return(list(class_vector,prob_vector))
  
}

def_peaks<-function(ecg,fs,gain){
  
        # The code uses an R code similar to Python Online and Offline ECG QRS Detector based 
        # on the Pan-Tomkins algorithm (https://github.com/c-labpl/qrs_detector). 
        # The code is a sample code for Physionet Challenge 2020 and not for any other experimental purposes. 
        # MIT License. Copyright (c) 2020. Andoni Elola (Universidad del Pais Vasco & Emory University).
  
  
        # Method responsible for extracting peaks from loaded ECG measurements data through measurements processing.
        # This implementation of a QRS Complex Detector is by no means a certified medical tool and should not be used in health monitoring. 
        # It was created and used for experimental purposes in psychophysiology and psychology.
        # You can find more information in module documentation:
        # https://github.com/c-labpl/qrs_detector
        # If you use these modules in a research project, please consider citing it:
        # https://zenodo.org/record/583770
        # If you use these modules in any other project, please refer to MIT open-source license.
        # If you have any question on the implementation, please refer to:
        # Michal Sznajder (Jagiellonian University) - technical contact (msznajder@gmail.com)
        # Marta lukowska (Jagiellonian University)
        # Janko Slavic peak detection algorithm and implementation.
        # https://github.com/c-labpl/qrs_detector
        # https://github.com/jankoslavic/py-tools/tree/master/findpeaks
        # 
        # MIT License
        # Copyright (c) 2017 Michal Sznajder, Marta Lukowska
        # 
        # Permission is hereby granted, free of charge, to any person obtaining a copy
        # of this software and associated documentation files (the "Software"), to deal
        # in the Software without restriction, including without limitation the rights
        # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        # copies of the Software, and to permit persons to whom the Software is
        # furnished to do so, subject to the following conditions:
        # The above copyright notice and this permission notice shall be included in all
        # copies or substantial portions of the Software.
        # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        # SOFTWARE.
  
  filter_lowcut<-0.001
  filter_highcut<-15
  filter_order<-1
  integration_window<-30
  findpeak_limit<-0.35
  findpeak_spacing<-100
  
  # Measurements filtering - 0-15 Hz band pass filter.
  filt<-butter(filter_order,c(2*filter_lowcut/fs,2*filter_highcut/fs),type = "pass")
  ecg_filtered<-filter(filt,ecg)

  ecg_filtered[1:5]<-ecg_filtered[5]
  
  # Derivative - provides QRS slope information.
  diff_ecg<-filter(c(1,-1),1,ecg_filtered)
  
  # Squaring - intensifies values received in derivative.
  squared_diff_ecg<-diff_ecg**2
  
  # Moving-window integration.
  integrated_ecg<-filter(rep(1,integration_window),1,squared_diff_ecg)
  
  # Fiducial mark - peak detection on integrated measurements.
  pks<-findpeaks(integrated_ecg, minpeakheight = findpeak_limit, minpeakdistance = findpeak_spacing)
  
  detected_peaks_indices<-sort(pks[,2])
  detected_peaks_values<-integrated_ecg[detected_peaks_indices]
  
  return(list(detected_peaks_indices,detected_peaks_values))
  
}


findpeaks <- function(x,nups = 1, ndowns = nups, zero = "0", peakpat = NULL,
                      minpeakheight = -Inf, minpeakdistance = 1,
                      threshold = 0, npeaks = 0, sortstr = FALSE)
{
  # Based on the code of the pracma package from CRAN
  
  # transform x into a "+-+...-+-" character string
  xc <- paste(as.character(sign(diff(x))), collapse="")
  xc <- gsub("1", "+", gsub("-1", "-", xc))
  # transform '0' to zero
  if (zero != '0') xc <- gsub("0", zero, xc)
  
  # generate the peak pattern with no of ups and downs
  if (is.null(peakpat)) {
    peakpat <- sprintf("[+]{%d,}[-]{%d,}", nups, ndowns)
  }
  
  # generate and apply the peak pattern
  rc <- gregexpr(peakpat, xc)[[1]]
  if (rc[1] < 0) return(NULL)
  
  # get indices from regular expression parser
  x1 <- rc
  x2 <- rc + attr(rc, "match.length")
  attributes(x1) <- NULL
  attributes(x2) <- NULL
  
  # find index positions and maximum values
  n <- length(x1)
  xv <- xp <- numeric(n)
  for (i in 1:n) {
    xp[i] <- which.max(x[x1[i]:x2[i]]) + x1[i] - 1
    xv[i] <- x[xp[i]]
  }
  
  # eliminate peaks that are too low
  inds <- which(xv >= minpeakheight & xv - pmax(x[x1], x[x2]) >= threshold)
  
  # combine into a matrix format
  X <- cbind(xv[inds], xp[inds], x1[inds], x2[inds])
  
  # eliminate peaks that are near by
  if (minpeakdistance < 1)
    warning("Handling 'minpeakdistance < 1' is logically not possible.")
  
  # sort according to peak height
  if (sortstr || minpeakdistance > 1) {
    sl <- sort.list(X[, 1], na.last = NA, decreasing = TRUE)
    X <- X[sl, , drop = FALSE]
  }
  
  # return NULL if no peaks
  if (length(X) == 0) return(c())
  
  # find peaks sufficiently distant
  if (minpeakdistance > 1) {
    no_peaks <- nrow(X)
    badpeaks <- rep(FALSE, no_peaks)
    
    # eliminate peaks that are close to bigger peaks
    for (i in 1:no_peaks) {
      ipos <- X[i, 2]
      if (!badpeaks[i]) {
        dpos <- abs(ipos - X[, 2])
        badpeaks <- badpeaks | (dpos > 0 & dpos < minpeakdistance)
      }
    }
    # select the good peaks
    X <- X[!badpeaks, , drop = FALSE]
  }
  
  return(X)
}

get_12ECG_features<-function(data,hea_data){

  tmp_hea<-unlist(strsplit(hea_data[1]," "))
  
  ptID<-tmp_hea[1]
  num_leads<-as.integer(tmp_hea[2])
  fs<-as.numeric(tmp_hea[3])
  
  gain_lead<-rep(0,num_leads)
  
  for (i in 1:num_leads) {
    tmp_hea<-unlist(strsplit(hea_data[i]," "))
    gain_lead[i]<-as.integer(unlist(strsplit(tmp_hea[3],"/"))[1])
    
  }
  
  # for testing, we included the mean age of 57 if the age is a NaN
  # This value will change as more data is being released
  age<-as.numeric(unlist(strsplit(hea_data[grepl("#Age", hea_data)],": "))[2])
  if (is.na(age)){age<-57}
  
  sex_string<-unlist(strsplit(hea_data[grepl("#Sex", hea_data)],": "))[2]
  if (sex_string=="Female") {
    
    sex<-1
    
  } else{
    
    sex<-0
    
  }
  
  #   We are only using data from lead1
  pks<-def_peaks(data[1,],fs,1)
  
  pks_indices<-pks[[1]]
  pks_values<-pks[[2]]
  
  rr<-pks_indices[2:length(pks_indices)]-pks_indices[1:length(pks_indices)-1]
  rr<-rr/fs
  
  feature_vector<-c()
  
  #Statistical moments
  feature_vector[1]<-mean(rr)
  feature_vector[2]<-mean(pks_values*gain_lead[1])
  
  feature_vector[3]<-sd(rr)
  feature_vector[4]<-sd(pks_values*gain_lead[1])
  
  feature_vector[5]<-skewness(rr)
  feature_vector[6]<-skewness(pks_values*gain_lead[1])
  
  feature_vector[5]<-kurtosis(rr)
  feature_vector[6]<-kurtosis(pks_values*gain_lead[1])
  
  #Age and sex
  feature_vector[7]<-age
  feature_vector[8]<-sex
  
  return(feature_vector)
  
}
