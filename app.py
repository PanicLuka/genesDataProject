import re
from flask import Flask, render_template, jsonify
import statistics
from collections import Counter
from decimal import Decimal, InvalidOperation
import pandas as pd
import numpy as np
import json
import math



app = Flask(__name__)


file_path = '9606_abund.txt' 
df1 = pd.read_csv(file_path, delim_whitespace=True)


file_path2 = '9606_gn_dom.txt' 
df2 = pd.read_csv(file_path2, delim_whitespace=True)

df2.rename(columns={'#Gn': 'Gn'}, inplace=True)

df1['Mean-copy-number'] = pd.to_numeric(df1['Mean-copy-number'], errors='coerce')




merged_df = pd.merge(df1, df2, on='Gn')

domain_avg_abundance = merged_df.groupby('Domain')['Mean-copy-number'].mean()

highest_avg_domain = domain_avg_abundance.idxmax()

highest_avg_abundance_value = domain_avg_abundance.max()


protein_stats = merged_df.groupby('Gn')['Mean-copy-number'].agg(['mean', 'std'])

merged_df['Percentile Rank'] = merged_df.groupby('Domain')['Mean-copy-number'].rank(pct=True)

output_file = 'combined_data_stats.csv'
merged_df.to_csv(output_file, index=False)





@app.route('/')
def index():
    return render_template('index.html')




def get_unique_values_with_one_occurrence(file_path, column_number):
    value_occurrences = Counter()
    with open(file_path, 'r') as file:
        next(file)
        for line in file:
            columns = re.split(r'\s+', line.strip())
            if len(columns) >= column_number:
                value = columns[column_number - 1]
                try:
                    value = Decimal(value)  
                    value_occurrences[value] += 1
                except InvalidOperation:
                    pass

    unique_values_with_one_occurrence = [value for value, count in value_occurrences.items() if count == 1]
    return unique_values_with_one_occurrence


def count_unique_values_with_one_occurrence(file_path, column_number):
    value_count = {}
    with open(file_path, 'r') as file:
        next(file)
        for line in file:
            columns = re.split(r'\s+', line.strip())
            if len(columns) >= column_number:
                value = columns[column_number - 1]
                value_count[value] = value_count.get(value, 0) + 1

    count_of_values_with_one_occurrence = sum(1 for count in value_count.values() if count == 1)
    return count_of_values_with_one_occurrence


def compute_mean_and_std_deviation(values):
    if len(values) > 0:
        mean = statistics.mean(values)
        std_dev = statistics.stdev(values)
        return mean, std_dev
    else:
        return None, None


@app.route('/numbers')
def get_count_unique_values_with_one_occurrence():
    file_path = '9606_abund.txt' 
    fourth_column_count = count_unique_values_with_one_occurrence(file_path, 4)
    return jsonify({"count": fourth_column_count})



@app.route('/get_stats')
def get_mean_and_std_deviation():
    file_path = '9606_abund.txt'  
    fourth_column_unique_values =  get_unique_values_with_one_occurrence(file_path, 4)
    mean, std_dev = compute_mean_and_std_deviation(fourth_column_unique_values)
    return jsonify({"mean": mean, "std_dev": std_dev})

@app.route('/get_highest_domain')
def get_highest_domain():
    return jsonify({"value": highest_avg_abundance_value})

@app.route('/get_mean_and_stdev_for_each_protein')
def get_mean_and_stdev_for_each_protein():
    protein_stats_dict = protein_stats.to_dict()

    if 'Protein Domain' in protein_stats.columns:
        for protein_domain, group in protein_stats.groupby('Protein Domain'):
            mean_value = group['Mean-copy-number'].mean()
            std_value = group['Mean-copy-number'].std()

            if isinstance(mean_value, float) and not math.isnan(mean_value):
                mean_value = round(mean_value, 2)
            else:
                mean_value = None

            if isinstance(std_value, float) and not math.isnan(std_value):
                std_value = round(std_value, 2)
            else:
                std_value = None

            protein_stats_dict[protein_domain] = {'mean': mean_value, 'std': std_value}

    json_response = json.dumps(protein_stats_dict, default=lambda x: None if x is None else float(x))

    response = app.response_class(
        response=json_response,
        status=200,
        mimetype='application/json'
    )

    return response


@app.route('/get_percentile_rank_for_each_domain')
def get_percentile_rank_for_each_domain():
    percentile_stats_dict = merged_df[['Gn', 'Domain', 'Mean-copy-number', 'Percentile Rank']].to_dict()
    return jsonify(percentile_stats_dict)



if __name__ == "__main__":
    app.run(debug=True)
