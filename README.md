# Example prediction code for R for the PhysioNet/CinC Challenge 2019

## Contents

This classifier uses two scripts:

* `run_12ECG_classifier.R` makes classifications on 12-Leads ECG data.  Add your prediction code to the `run_12ECG_classifier` function. `load_12ECG_model.R` loads model weights, etc. for making classifications.  To reduce your code's run time, add any code to the `load_12ECG_model` function that you only need to run once, such as loading weights for your model.
* `driver.R` calls `load_12ECG_model` once and `run_12ECG_classifier` many times. It also performs all file input and output.  **Do not** edit this script -- or we will be unable to evaluate your submission.

Check the code in these files for the input and output formats for the `load_12ECG_model` and `run_12ECG_classifier` functions.

## Use

You can run this classifier by installing the packages in the `requirements.txt` file and running

    Rscript driver.R input_directory output_directory

where `input_directory` is a directory for input data files and `output_directory` is a directory for output classification files. The PhysioNet/CinC 2020 webpage provides a training database with data files and a description of the contents and structure of these files.

## Submission

The `driver.R`, `run_12ECG_classifier.R`, and `get_12ECG_features.R` scripts need to be in the base or root path of the Github repository. If they are inside a subfolder, then the submission will fail.

## Details
“The baseline classifiers are simple Random Forest. They use statistical moments of heart rate that we computed from the WFDB signal file (the `.mat` file) and demographic data taken directly from the WFDB header file (the `.hea` file) as predictors. 

The code uses an R code similar to Python Online and Offline ECG QRS Detector based on the Pan-Tomkins algorithm (https://github.com/c-labpl/qrs_detector). The code is a sample code for Physionet Challenge 2020 and not for any other experimental purposes. 
MIT License. Copyright (c) 2020. Andoni Elola (Universidad del Pais Vasco & Emory University).
