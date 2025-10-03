import pandas as pd # handles table of data
import matplotlib.pyplot as plt # makes plots and graphs
import matplotlib.ticker as mticker # used to fine-tune axis labels and tick marks
import seaborn as sns # higher-level plotting library built on top pf matplot --> nicer looking graphs 
import os # lets the script create folders
import re # regular expressions, used for pattern matching (here: extracting numbers and units from a string)
import sys # allows to interact with command-line arguments
import numpy as np
from scipy.optimize import curve_fit # used to fit a curve over a plot, in this case hill-fit model is used 

# get input file via the command line 
myfile = sys.argv[1]
df = pd.read_csv(myfile)

# generate a mean NB file that uses the same experiment name of input file but exchanges dose response to mean NB count 
mean_NB_file = myfile.replace("_dose_response.csv", "_mean_NB_count.csv")
# do the same for the EC50_file
EC50_file = myfile.replace("_dose_response.csv", "_NB_EC50_results.csv")

# define helper function to convert the concentration string into a float and transfers the units, so all concentration points are in the same unit (nM)
def to_nM(conc_str):
    if pd.isna(conc_str):
        return None
    conc_str = conc_str.replace(" ", "").lower()
    num = float(re.findall(r'[\d\.]+', conc_str)[0])
    if "µm" in conc_str or "um" in conc_str:
        return num * 1000
    elif "nm" in conc_str:
        return num
    
# calculate mean per well label 
def get_mean_per_well(data, output_file):
    # group data by well label and compund ID and calculate average number of network bursts for each well-compound pair
    mean_NB = data.groupby(["Well Label", "Compound ID"])["Network Burst Count"].mean().reset_index()
    # rename column
    mean_NB.rename(columns={"Network Burst Count": "Mean Network Burst Count"}, inplace=True)
   # save mean values to the output file
    mean_NB.to_csv(output_file, index=False)
    print(f"Mean values saved to {output_file}")
    return mean_NB

# take value of each well (for compounds) and divide by the value of the baseline_media well (e.g A1 compound 1 divided by A1 baseline)
def get_normalised_NB(mean_data, output_file):
    # extract baseline values per well: find rows where the compound ID is baseline_media and take the columns well label and mean network burst count
    baseline = mean_data[mean_data["Compound ID"] == "baseline_media"][["Well Label", "Mean Network Burst Count"]]
    # rename the mean network burst column to Baseline
    baseline.rename(columns={"Mean Network Burst Count": "Baseline"}, inplace=True)
    
  # Merge baseline with all the other rows for the same well label
    mean_data = mean_data.merge(baseline, on="Well Label", how="left")

    # avoid dividing by 0 by turning 0 values of baseline into NaN 
    mean_data.loc[mean_data["Baseline"] == 0, "Baseline"] = np.nan

    
    # Normalize: divide every measurement in the well by the baseline 
    mean_data["Normalized Network Burst Count"] = mean_data["Mean Network Burst Count"] / mean_data["Baseline"]
    
    # Drop Baseline column for cleaner CSV, otherwise it will be saved in the csv file as well
    mean_data = mean_data.drop(columns=["Baseline"])
    
    # save results in the same output file 
    mean_data.to_csv(output_file, index=False)
    print(f"Normalized values added and saved to '{output_file}'")
    
    return mean_data

# function that calculates mean, SD and %CV of each compound/concentration pair
def get_SD_CV(normalized_data, output_file):

    # create a copy to avoid modifying original data 
    temp = normalized_data.copy()
    # use a regular expression to split compound ID in compound name and concentration
    temp[['Compound_name', 'Concentration_str']] = temp['Compound ID'].str.extract(
        r'([a-zA-Z0-9]+)[\s_]+([\d\.]+\s*[µu]?M|[\d\.]+\s*nM)?', expand=True
    )
    # apply the to_nM helper function to the concentration string
    temp['Concentration_nM'] = temp['Concentration_str'].apply(to_nM)

    # Group by compound + concentration and calculate mean + SD
    summary = (
        temp.groupby(['Compound_name', 'Concentration_nM', 'Concentration_str'])
        ['Normalized Network Burst Count']
        .agg(['mean', 'std'])
        .reset_index()
    )

    # Calculate the %CV and convert infinitive %CV into NaN
    summary['CV_percent'] = (summary['std'] / summary['mean']) * 100
    summary['CV_percent'] = summary['CV_percent'].replace([np.inf, -np.inf], np.nan)

    # generate a summary stats file and cave data 
    summary_file = output_file.replace("_mean_NB_count.csv", "_normalized_NB_stats.csv")
    summary.to_csv(summary_file, index=False)
    print(f"Summary stats saved to '{summary_file}'")

    return summary

# create a EC50 fitting helper function 1: hill model for dose response curve
def four_param_logistic(x, bottom, top, logEC50, hill_slope):
    return bottom + (top - bottom) / (1 + 10**((logEC50 - np.log10(x)) * hill_slope))

# create EC50 fitting helper function 2: fit the EC50 curve to the data using curve_fit
def fit_EC50(xdata, ydata):
    # only continued if at least 4 values 
    if len(xdata) < 4:
        return None, None, None
    # initial guess for the curve fitting 
    p0 = [min(ydata), max(ydata), np.log10(np.median(xdata)), 1.0]
    # fit the curve and calculate EC50, return parameters or None if the fitting failed 
    try:
        popt, pcov = curve_fit(
            four_param_logistic,
            xdata,
            ydata,
            p0=p0,
            maxfev=5000
        )
        bottom, top, logEC50, hill_slope = popt
        EC50 = 10**logEC50
        return popt, pcov, EC50
    except RuntimeError:
        return None, None, None

# now plot normalised values for each compound --> create a folder where outputs can be saved in
def plot_each_compound(normalized_data, output_folder="plots", ec50_file=EC50_file):
    # Make the output folder if it doesn´t exist already
    os.makedirs(output_folder, exist_ok=True)
    
    # Exclude baseline and control rows
    plot_data = normalized_data[~normalized_data["Compound ID"].str.contains("baseline_media|CTRL_veh", case=False, na=False)].copy()

    # a regular expression is used to split the compound ID in the compound and concentration 
    plot_data[["Compound", "Concentration_str"]] = plot_data["Compound ID"].str.extract(
        r'([a-zA-Z0-9]+)[\s_]+([\d\.]+\s*[µu]?M|[\d\.]+\s*nM)', expand=True
    )
    # again: use to_nM to convert units 
    plot_data["Concentration_nM"] = plot_data["Concentration_str"].apply(to_nM)
    
    # Sort by concentration so lines connect properly
    plot_data = plot_data.sort_values(by="Concentration_nM")
    # extract the unique compounds and start a list for the EC50 results 
    compounds = plot_data["Compound"].unique()
    ec50_results = []
    
    # loop over the compounds to get a subset of compounds for the plotting
    for compound in compounds:
        subset = plot_data[plot_data["Compound"] == compound]

        # skip compounds where ALL mormalised values are 0 
        if (subset["Normalized Network Burst Count"] <= 0).all() or subset["Normalized Network Burst Count"].isna().all():
            print(f"Skipping compound {compound} (all values are 0 or NaN)")
            continue

        # only continue with valid wells --> where all values are >0
        valid_wells = subset.groupby("Well Label")["Normalized Network Burst Count"].transform(lambda x: (x > 0).any())
        subset = subset[valid_wells]
        
        # create a scatter plot for the valid wells with conc. (log) on x-axis, normalized burst count on y-axis
        plt.figure(figsize=(7,5))
        sns.scatterplot(
            data=subset, 
            x="Concentration_nM",
            y="Normalized Network Burst Count",
            hue="Well Label", #each well labale gets a different color
            marker="o", #object that marks the points in the plot
            linewidth=1,
            alpha=0.6
        )

        # Plot mean +/- SD as error bars 
        summary = (
            subset.groupby("Concentration_nM")["Normalized Network Burst Count"]
            .agg(["mean", "std"])
            .reset_index()
        )

        plt.errorbar(
            summary["Concentration_nM"],
            summary["mean"],
            yerr=summary["std"],
            fmt="o",
            capsize=3,
            label="Mean ± SD"
        )
        # use logarithmic scale for x axis 
        plt.xscale("log")
        
        # BUT: define custom ticks for better readability
        ticks = [1, 10, 100, 1000, 10000]  # nM values
        tick_labels = ["1 nM", "10 nM", "100 nM", "1 µM", "10 µM"] # rename ticks into the original concentrations
        plt.xticks(ticks, tick_labels)
        
        # name plot, axis etc.
        plt.title(f"{compound} - Normalized Network Burst Count")
        plt.xlabel("Concentration (nM)")
        plt.ylabel("Normalized Network Burst Count")
        
        # prepare data for EC50 fitting
        xdata = summary["Concentration_nM"].values
        ydata = summary["mean"].values
        mask = (xdata > 0) & (~np.isnan(ydata))
        xdata, ydata = xdata[mask], ydata[mask]
        popt, pcov, EC50 = fit_EC50(xdata, ydata)
        
        # only of EC50 is sucessfull (so not None): plot curve and store the parameters
        if EC50 is not None:
            xfit = np.logspace(np.log10(min(xdata)), np.log10(max(xdata)), 200)
            yfit = four_param_logistic(xfit, *popt)
            plt.plot(xfit, yfit, "r-", label=f"4PL fit (EC50={EC50:.1f} nM)")
            print(f"{compound} EC50 = {EC50:.2f} nM")
            bottom, top, logEC50, hill_slope = popt
            ec50_results.append({
                "Compound": compound,
                "Bottom": bottom,
                "Top": top,
                "HillSlope": hill_slope,
                "EC50_nM": EC50
            })
        # otherwise: store NaN and print message that it failed 
        else:
            print(f"EC50 fitting failed for {compound}")
            ec50_results.append({
                "Compound": compound,
                "Bottom": np.nan,
                "Top": np.nan,
                "HillSlope": np.nan,
                "EC50_nM": np.nan
            })
        # add a legend
        plt.legend(title="Well Label", bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.tight_layout()
        
        # Save each plot to the output folder and close it to save memory 
        plt.savefig(f"{output_folder}/{compound}_normalized_NB_plot.png")
        plt.close()
        print(f"Plot saved for compound: {compound}")

    # save EC50 results in a new file
    ec50_df = pd.DataFrame(ec50_results)
    ec50_df.to_csv(ec50_file, index=False)
    print(f"EC50 + fit parameters saved to {ec50_file}")

# execute workflow
mean_data = get_mean_per_well(df, mean_NB_file)
normalized_data = get_normalised_NB(mean_data, mean_NB_file)
summary = get_SD_CV(normalized_data, mean_NB_file)
plot_each_compound(normalized_data, output_folder="plots", ec50_file=EC50_file)