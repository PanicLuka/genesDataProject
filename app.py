import re
from flask import Flask, render_template, jsonify
import statistics
from collections import Counter
from decimal import Decimal, InvalidOperation
import pandas as pd
import numpy as np
import json
import math
from collections import Counter


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
    return jsonify({"value": highest_avg_domain})

@app.route('/get_mean_and_stdev_for_each_protein')
def get_mean_and_stdev_for_each_protein():
    gene_data = pd.read_csv('9606_abund.txt', delimiter='\t')

    domain_data = pd.read_csv('9606_gn_dom.txt', delimiter='\t')

    domain_data.rename(columns={'#Gn': 'Gn'}, inplace=True)

    gene_data['Mean-copy-number'] = pd.to_numeric(gene_data['Mean-copy-number'], errors='coerce')


    combined_data = pd.merge(gene_data, domain_data, on='Gn')

    grouped_data = combined_data.groupby(['Gn', 'Domain'])['Mean-copy-number'].mean().reset_index()
    grouped_data['Standard-Deviation'] = combined_data.groupby(['Gn', 'Domain'])['Mean-copy-number'].std().reset_index()['Mean-copy-number']


    grouped_data['Percentile-Rank'] = grouped_data.groupby('Domain')['Mean-copy-number'].rank(pct=True) * 100

    grouped_data_dict = grouped_data.to_dict()

    for index, row in grouped_data.iterrows():
        gene = row['Gn']
        domain = row['Domain']
        mean_abundance = row['Mean-copy-number']
        std_deviation = row['Standard-Deviation']
        percentile_rank = row['Percentile-Rank']
        
        print(f"Gene: {gene}\tDomain: {domain}\tMean Abundance: {mean_abundance:.2f}\tStd Deviation: {std_deviation:.2f}\tPercentile Rank: {percentile_rank:.2f}%")



@app.route('/get_percentile_rank_for_each_gene')
def get_percentile_rank_for_each_domain():
    data = pd.read_csv('9606_abund.txt', delimiter='\t')

    unique_data = data.drop_duplicates(subset=['Gn', 'Mean-copy-number'])

    unique_data['Rank'] = unique_data['Mean-copy-number'].rank(method='average', ascending=True)

    grouped = unique_data.groupby('Gn')['Rank'].max()
    total_genes = grouped.count()
    gene_percentiles = (grouped - 1) / (total_genes - 1) * 100

    gene_dict = gene_percentiles.to_dict()
    
    return jsonify(gene_dict)

   



if __name__ == "__main__":
    app.run(debug=True)
