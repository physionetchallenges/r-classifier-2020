#Load model from disk
load_12ECG_model<-function(){
  
  load("modelo.RData")
  
  return(model)
  
}