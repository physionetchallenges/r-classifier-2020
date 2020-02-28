library(R.matlab)

source(file.path(getwd(),"run_12ECG_classifier.R"))
source(file.path(getwd(),"load_12ECG_model.R"))

load_challenge_data<-function(filename){
  
  x<-readMat(filename)
  data<-x[["val"]]
  
  input_hea_file<-gsub(pattern = "\\.mat$", ".hea", filename)
  hea_data<-readLines(input_hea_file)
  
  return(list(data,hea_data))
  
}

get_classes<-function(input_directory,files){
  
  classes<-c()
  
  for (i in 1:length(files)) {
    
    hea_file<-gsub(pattern = "\\.mat$", ".hea", file.path(input_directory,files[i]))
    hea_act<-readLines(hea_file)
    
    hea_act<-hea_act[grepl("#Dx", hea_act)]
    tmp_hea<-unlist(strsplit(hea_act," "))[[2]]
    tmp_hea<-unlist(strsplit(tmp_hea,","))
    classes<-append(classes,trimws(tmp_hea))
    
  }
  
  return(sort(unique(classes)))
  
}

cinput<-commandArgs(trailingOnly = TRUE)

if (length(cinput)!=2){
  
	stop("Not enough input arguments")
  
} else{

  input_directory=cinput[1]
  output_directory=cinput[2]
  
}

files<-list.files(input_directory,pattern = "\\.mat$")

classes<-get_classes(input_directory,files)

#Load the model
print("Loading 12ECG model...")
model<-load_12ECG_model()

# Iterate over files.
print("Extracting 12ECG features...")
for (i in 1:length(files)) {
  #print(i)
  print(paste(i,length(files),sep = "/"))
  l<-load_challenge_data(file.path(input_directory,files[i]))
  
  data_act<-l[[1]]
  hea_act<-l[[2]]
  
  res_classifier<-run_12ECG_classifier(data_act,hea_act, classes, model)
  
  name_act<-gsub(pattern = "\\.mat$", ".csv", files[i])
  name_act<-file.path(output_directory,name_act)
  
  #Save results
  write(paste("#",gsub(pattern = "\\.mat$", " ", files[i]),sep = ""),file=name_act)
  write(classes,file=name_act,append = TRUE,ncolumns = length(classes),sep = ",")
  write(res_classifier[[1]],file=name_act,append = TRUE,ncolumns = length(classes),sep = ",")
  write(res_classifier[[2]],file=name_act,append = TRUE,ncolumns = length(classes),sep = ",")
  
  
}
