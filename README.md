# My Python Project
This is my script for the python course. The script takes a csv file and reads it. The relevant columns used are Well Label, Compound ID, and Network Burst Count. 
The script calculates the mean Network Burst Count for each Well Label.  After that, the second function will normalise the values. So each mean for each compound is divided by the mean of the baseline value in the corresponding well (example: value of well A1 of compound 1 conc. 1 is divided by value of well A1 of baseline and so on). The results will be saved in a new csv file.
From the normalized values, the function called get_SD_CV will calculate the mean, standard deviation and coefficient of variation for each concentration of the compound. Again, the results will be saved in a new csv file.
Then, I am using helper functions to fit a curve over the plots I am creating later on. The curve fitting (four_param_logistic) follows the hill-model (four parameter logistics) and is used to calculate the EC50 (fit_ec50). 
The last function will now create a scatter plot, where the concentration is on the y axis (logarithmically) and the the network burst count is on the x axis. The mean and SD of each concentration point are added to the plot as an errorbar plot. The curve_fit model of scipy is used to get the EC50 value. The curve shows up in the plot and the EC50 is added to the hue. The results of the curve fitting are also saved to a third csv file.
For each compound is a seperate plot created and all plots are saved in a folder called _plots_, where the plots are named like this: compound_normalized_NB_plot.png

## How to run the script
python get_NB_count_plot.py test_data_dose_response.csv

## Arguments
Pass the input file you want to use (test_data_dose_response.csv) to the command line.
To run the script please download the following modules:
* pandas
* matplotlib.pyplot
* matplotlib.ticker
* seaborn 
* sys, os, and re
* numpy
* curve_fit from scipy_optimize

## Test Data
I added a file called test_data_dose_response.csv to my project. Please make sure that you have the python script and the test_data_dose_response.csv in the same folder if you want to run the program.

## Output
The output is a folder called _plots_ that includes a normalized plot for each compound I tested. In addition to that, you will also get 3 new csv files: one for the mean Network Burst count + the normalization (called test_data_mean_NB_count.csv), one for the results of the mean, SD and %CV (called test_data_normalized_NB_stats.csv) and one for the curve fitting results (called test_data_NB_EC50_results.csv).