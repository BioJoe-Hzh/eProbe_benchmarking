import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import operator
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy import stats
import statsmodels.api as sm

# Load various dictionaries from files to hold different metrics (depth, GC content, Tm, etc.)
depth_dict = {}
with open("filtered.averaged_depths.norm.bed") as f1: # result from 4.6.1 in 4_statistic.sh
    for line in f1.readlines():
        pos = "-".join(line.split()[:3])
        depth_dict[pos] = float(line.split()[3])

GC_dict = {}
Tm_dict = {}
with open("all.probe.GC.Tm.txt") as f2: # result from probe generation with eProbe (SNP_filter.py); each row records the pos and feature for a probe
    for line in f2.readlines():
        pos = "-".join(line.split()[:3])
        GC_dict[pos] = float(line.split()[3])
        Tm_dict[pos] = float(line.split()[4])

Hp_dict = {}
with open("all.probe.hairpin.txt") as f3: # result from probe generation with eProbe (SNP_filter.py)
    for line in f3.readlines():
        pos = "-".join(line.split()[:3])
        Hp_dict[pos] = float(line.split()[3])

Dm_dict = {}
with open("all.probe.dimer.txt") as f4: # result from probe generation with eProbe (SNP_filter.py)
    for line in f4.readlines():
        pos = "-".join(line.split()[:3])
        Dm_dict[pos] = float(line.split()[3])

Dust_dict = {}
with open("all.probe.complexity.txt") as f5: # result from probe generation with eProbe (SNP_filter.py)
    for line in f5.readlines():
        pos = "-".join(line.split()[:3])
        Dust_dict[pos] = float(line.split()[3])

# Function to merge multiple dictionaries into a pandas DataFrame
def merge_dicts_to_df(dicts, col_names):
    """
    Merge dictionaries with shared keys into a DataFrame with specified column names.

    Parameters:
    - dicts: List of dictionaries to merge.
    - col_names: List of column names for the DataFrame including 'chr', 'start', 'end'.

    Returns:
    - A pandas DataFrame with combined data from dictionaries.
    """
    shared_keys = set(dicts[0].keys())
    for d in dicts[1:]:
        shared_keys.intersection_update(d.keys())

    # Split the shared keys into 'chr', 'start', 'end'
    chr_start_end = [key.split('-') for key in shared_keys]
    df = pd.DataFrame(chr_start_end, columns=['chr', 'start', 'end'])

    # Add the values from each dictionary as new columns
    for d, col_name in zip(dicts, col_names[3:]):
        df[col_name] = df.apply(lambda row: d[f"{row['chr']}-{row['start']}-{row['end']}"], axis=1)
    
    return df

# Function to remove rows that are outside of a specified range
def remove_rows_outside_threshold(df, column, min_threshold, max_threshold):
    """
    Remove rows from a DataFrame where the value in a specified column is outside a given range.
    """
    return df[(df[column] >= min_threshold) & (df[column] <= max_threshold)]

# Function to calculate Pearson correlation between two columns
def calculate_pearson_correlation(df, col_name1, col_name2):
    """
    Calculate the Pearson correlation coefficient between two columns in a DataFrame.
    """
    return df[[col_name1, col_name2]].corr(method='pearson').iloc[0, 1]

# Function to remove outliers using the IQR method
def remove_outliers(df, column):
    """
    Remove outliers from a DataFrame based on the IQR method.
    """
    Q1 = df[column].quantile(0.25)
    Q3 = df[column].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return df[(df[column] >= lower_bound) & (df[column] <= upper_bound)]

# Function to perform and plot polynomial regression
def polynomial_regression(df, x_col, y_col, degree=2):
    """
    Perform polynomial regression and plot the results.
    """
    x = df[x_col].values.reshape(-1, 1)
    y = df[y_col].values

    poly_features = PolynomialFeatures(degree=degree)
    x_poly = poly_features.fit_transform(x)

    # Fit polynomial regression
    poly_model = LinearRegression()
    poly_model.fit(x_poly, y)
    y_poly_pred = poly_model.predict(x_poly)

    r_squared = poly_model.score(x_poly, y)
    print(f'R-squared: {r_squared}')

    # Plot the original data and the polynomial regression line
    plt.scatter(x, y, alpha=0.5)
    sorted_zip = sorted(zip(x, y_poly_pred), key=operator.itemgetter(0))
    x, y_poly_pred = zip(*sorted_zip)
    plt.plot(x, y_poly_pred, color='m')
    plt.title('Polynomial Regression')
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.show()

    return poly_model

# Function to perform polynomial regression using Statsmodels and print detailed statistics
def polynomial_regression_statsmodels(df, x_col, y_col, degree=2):
    """
    Perform polynomial regression using statsmodels and print model summary.
    """
    x = df[x_col].values.reshape(-1, 1)
    y = df[y_col].values

    poly_features = PolynomialFeatures(degree=degree)
    x_poly = poly_features.fit_transform(x)

    # Add intercept for statsmodels
    x_poly_with_intercept = sm.add_constant(x_poly)

    # Fit OLS model
    model = sm.OLS(y, x_poly_with_intercept).fit()
    print(model.summary())

    return model

# Function to bin data, perform regression, and plot
def bin_and_regress(df, column1, column2, start, end, bin_size, min_samples, output, line_color, regression_type):
    """
    Perform binning, outlier removal, regression, and plot with options to customize scatter plot.
    """
    bins = np.arange(start, end, bin_size)
    bins = np.append(bins, end)

    df['bin'] = pd.cut(df[column1], bins=bins, include_lowest=True, right=False)

    bin_centers = []
    means = []

    for bin_left, group in df.groupby('bin'):
        if len(group) >= min_samples:
            q1 = group[column2].quantile(0.25)
            q3 = group[column2].quantile(0.75)
            iqr = q3 - q1
            filtered_group = group[(group[column2] >= q1 - 1.5 * iqr) & (group[column2] <= q3 + 1.5 * iqr)]
            if len(filtered_group) >= min_samples:
                bin_centers.append(bin_left.mid)
                means.append(filtered_group[column2].mean())

    if len(bin_centers) > 1:
        if regression_type == 'linear':
            slope, intercept, r_value, p_value, std_err = stats.linregress(bin_centers, means)
            r_squared = r_value ** 2

            x_pred = np.linspace(np.min(bin_centers), np.max(bin_centers), 50)
            y_pred = slope * x_pred + intercept

            # Scatter plot with linear regression line
            plt.figure(figsize=(10, 6))
            sns.scatterplot(x=bin_centers, y=means, s=150, color="grey")
            plt.plot(x_pred, y_pred, color=line_color, linewidth=3)
            plt.xlabel(column1)
            plt.ylabel(f'Mean of {column2}')
            plt.savefig(output, dpi=800)
            plt.show()

            print(f"R-squared: {r_squared:.4f}")
            return r_squared

        elif regression_type == 'polynomial':
            coeffs = np.polyfit(bin_centers, means, 2)
            polynomial = np.poly1d(coeffs)

            # Scatter plot with polynomial regression line
            plt.figure(figsize=(10, 6))
            plt.plot(bin_centers, polynomial(bin_centers), linewidth=3, color=line_color)
            sns.scatterplot(x=bin_centers, y=means, s=150, color="grey")
            plt.xlabel(column1)
            plt.ylabel(f'Mean of {column2}')
            r_squared_poly = r2_score(means, polynomial(bin_centers))
            print(f"R-squared (Polynomial): {r_squared_poly:.4f}")
            plt.savefig(output, dpi=800)
            plt.show()

            return None

    else:
        print("Not enough data points for regression.")
        return None

# Merge the dictionaries into a DataFrame
dict_list = [depth_dict, GC_dict, Tm_dict, Hp_dict, Dm_dict, Dust_dict]
col_name = ['chr', 'start', 'end', "Depth", "GC", "Tm", "Hp", "Dm", "Dust"]
merge_df = merge_dicts_to_df(dict_list, col_name)

# Filter data based on GC content between 40 and 60 (filter threshold)
filteredGC_merge_df = remove_rows_outside_threshold(merge_df, "GC", 40, 60)

# Perform regression and plot for various metrics
bin_and_regress(filteredGC_merge_df, "GC", "Depth", 40, 60, 2, 100, "GC_plot.pdf", '#CC5164', regression_type='linear')
bin_and_regress(filteredGC_merge_df, "Tm", "Depth", 71, 81, 1, 50, "Tm_plot.pdf", '#477CA1', regression_type='linear')
bin_and_regress(filteredGC_merge_df, "Hp", "Depth", 15, 55, 5, 100, "Hp_plot.pdf", '#43807C', regression_type='linear')
bin_and_regress(filteredGC_merge_df, "Dm", "Depth", 750, 2250, 200, 100, "Dm_plot.pdf", '#8C8CB4', regression_type='polynomial')
bin_and_regress(filteredGC_merge_df, "Dust", "Depth", 0.5, 1.7, 0.1, 100, "Dust_plot.pdf", '#DA9B60', regression_type='polynomial')

# Additional polynomial regression analysis
polynomial_regression(filteredGC_merge_df, "GC", "Depth", degree=1)
print(calculate_pearson_correlation(filteredGC_merge_df, "GC", "Depth"))
polynomial_regression_statsmodels(filteredGC_merge_df, "GC", "Depth", degree=1)

polynomial_regression(filteredGC_merge_df, "Tm", "Depth", degree=1)
print(calculate_pearson_correlation(filteredGC_merge_df, "Tm", "Depth"))
polynomial_regression_statsmodels(filteredGC_merge_df, "Tm", "Depth", degree=1)

polynomial_regression(filteredGC_merge_df, "Hp", "Depth", degree=1)
print(calculate_pearson_correlation(filteredGC_merge_df, "Hp", "Depth"))
polynomial_regression_statsmodels(filteredGC_merge_df, "Hp", "Depth", degree=1)

polynomial_regression(filteredGC_merge_df, "Dm", "Depth", degree=2)
print(calculate_pearson_correlation(filteredGC_merge_df, "Dm", "Depth"))
polynomial_regression_statsmodels(filteredGC_merge_df, "Dm", "Depth", degree=2)

polynomial_regression(filteredGC_merge_df, "Complexity", "Depth", degree=2)
print(calculate_pearson_correlation(filteredGC_merge_df, "Dust", "Depth"))
polynomial_regression_statsmodels(filteredGC_merge_df, "Dust", "Depth", degree=2)
