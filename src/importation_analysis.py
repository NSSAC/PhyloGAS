# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.1
#   kernelspec:
#     display_name: escape
#     language: python
#     name: python3
# ---

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#import epiweeks as ew


# %%
#read in table
df = pd.read_csv('../data/importations/sarscov2_clusters_2024_11_12_filtered.tsv.gz', sep='\t')
df.head()
#region unique value counts
df['region'].value_counts()

# %%
df.head()

# %% [markdown]
#

# %%
df.columns

# %%
#create genbank_seed column by splitting samples column based on ',' then based on '|' and taking the second element
df['genbank_seed'] = df['samples'].str.split(',').str[0].str.split('|').str[1]
#create strain_seed column by splitting samples column based on ',' then based on '|' and taking the first element
df['strain_seed'] = df['samples'].str.split(',').str[0].str.split('|').str[0]
df.head()

# %%
from pango_aliasor.aliasor import Aliasor
from typing import List, Dict
#of the variants in base_variants make a dictionary of all the subvariants that point to them
#the recombinant option is based on a request and likely not to work well.
def make_variant_base_map(base_variants: List[str], recombinant=False) -> Dict[str, str]:
    namer = Aliasor()
    namer.enable_expansion()
    all_rules = namer.partition_focus(base_variants, recombinant=recombinant)
    lineage_base_map = {k: v for v, ks in all_rules.items() for k in ks}
    
    for b in base_variants:
        if b not in lineage_base_map:
            lineage_base_map[b] = b
            
    return lineage_base_map


# %%
#regularize the variant names
#target_list=['B.1.617.2','B.1.1.7']
target_list = ['B.1.1.7', 'B.1.617.2', 'BA.1', 'BA.2', 'BA.4', 'BA.5', 'XBB', 'XBB.1.5', 'XBB.1.16', 'XBB.1.9', 'BA.2.86', 'KP.3', 'LP.8', 'XEC']
replace_map = make_variant_base_map(target_list)
df['annotation_3']=df['annotation_2'].map(replace_map)

#make a week year column
df['earliest_date'] = pd.to_datetime(df['earliest_date'], errors='coerce')
df['week'] = df['earliest_date'].dt.strftime('%U')
df['year'] = df['earliest_date'].dt.strftime('%Y')
df['week_year'] = df['year'] + '-' + df['week']
df['week_year'] = pd.to_datetime(df['week_year'] + '-0', format='%Y-%U-%w').dt.strftime('%Y-%U')


# %%
group_by_df.head()

# %%
# Group by annotation_3 and create a boxplot for each group per region
for region in df['region'].unique():
    # Filter df for the region
    region_df = df[df['region'] == region]
    #create a bin size category based on sample_count
    region_df['bin_size'] = pd.cut(region_df['sample_count'], bins=[0, 10, 50, 100, 500, 1000, 5000, 10000, 50000], labels=['0-10', '10-50', '50-100', '100-500', '500-1000', '1000-5000', '5000-10000', '10000+'])
    #group by bin size and count the sum of sample_count per bin size.
    group_by_df = region_df.groupby(['bin_size','annotation_3'], as_index=False)['sample_count'].sum()    
    # Plot cluster size distribution as a boxplot grouped by annotation_3
    plt.figure(figsize=(12, 8))
    sns.boxplot(x='bin_size', y='sample_count', data=group_by_df, showfliers=False)
    plt.title(f'Cluster Size Distribution for {region}')
    plt.xlabel('Bin Size')
    plt.ylabel('Samples per variant wave')
    #plt.yscale('log')  # Log scale for better visualization
    plt.xticks(rotation=45, ha='right')
    
    # Show the plot
    plt.tight_layout()
    plt.show()


# %%
#sort by sample_count
df.sort_values('sample_count', ascending=False).head()
target_lineage="B.1.617.2"
#target_list = ['B.1.1.7', 'B.1.617.2', 'BA.1', 'BA.2', 'BA.4', 'BA.5', 'XBB', 'XBB.1.5', 'XBB.1.16', 'XBB.1.9', 'BA.2.86', 'KP.3', 'LP.8', 'XEC']
target_list = ['B.1.1.7', 'B.1.617.2', 'BA.1', 'BA.2', 'BA.4']
#for each region plot the cluster size distribution
for region in df['region'].unique():
    for target_lineage in target_list:
        #filter df for region
        region_df = df[(df['region'] == region) & (df['annotation_3'] == target_lineage)]
        #sort by sample_count
        region_df = region_df.sort_values('sample_count', ascending=False)
        #plot cluster size distribution
        plt.figure(figsize=(10, 6))
        sns.histplot(region_df['sample_count'], bins=30, kde=False, stat='count', color='blue')
        plt.title(f'Cluster Size Distribution for {region} - {target_lineage}')
        plt.xlabel('Cluster Size')
        plt.ylabel('Frequency')
        #log log scale
        plt.yscale('log')
       
        #plt.savefig(f'../figures/importations/cluster_size_distribution_{region}_{target_lineage}.png')
        #plt.close()
        plt.show()

    plt.show()

# %%
#sort by sample_count
df.sort_values('sample_count', ascending=False).head()
target_lineage="B.1.617.2"
target_list = ['B.1.1.7', 'B.1.617.2', 'BA.1', 'BA.2', 'BA.4']

# Create a 5x5 grid for the plots
fig, axes = plt.subplots(nrows=6, ncols=5, figsize=(20, 20))
axes = axes.flatten()  # Flatten the axes array for easier indexing

# Iterate over regions and target lineages
for idx, (region, target_lineage) in enumerate([(r, t) for r in df['region'].unique() for t in target_list]):
    if idx >= len(axes):  # Ensure we don't exceed the number of subplots
        break
    
    # Filter df for the region and target lineage
    region_df = df[(df['region'] == region) & (df['annotation_3'] == target_lineage)]
    
    # Create a bin size category based on sample_count
    region_df['bin_size'] = pd.cut(region_df['sample_count'], bins=[0, 10, 50, 100, 500, 1000, 5000, 10000, 50000], 
                                    labels=['0-10', '10-50', '50-100', '100-500', '500-1000', '1000-5000', '5000-10000', '10000+'])
    
    # Group by bin size and count the sum of sample_count per bin size
    group_by_df = region_df.groupby('bin_size', as_index=False)['sample_count'].sum()
    
    # Plot the bin size vs total sample count
    ax = axes[idx]
    sns.barplot(x='bin_size', y='sample_count', data=group_by_df, color='blue', ax=ax)
    ax.set_title(f'{region} - {target_lineage}')
    ax.set_xlabel('Cluster Size')
    ax.set_ylabel('Total Sample Count')
    ax.set_yscale('log')
    ax.tick_params(axis='x', rotation=45)

# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()


# %%
#for Region Washington and annotation_2 = B.1.617.2 plot the number of rows per week
def plot_importations(df, region, annotation):
    df2 = df[df['annotation_3'] == annotation]
    df2 = df2[df2['region'] == region]
    df3 = df2.groupby(['week_year']).size().reset_index(name='importations')
    df3['seq_count']=df2.groupby(['week_year'])['sample_count'].sum().reset_index(name='sample_count')['sample_count']
    importations = df3['importations'].sum()
    df3 = df3.sort_values(by='week_year')
    df3.plot(x='week_year', y='importations', kind='bar')
    plt.xlabel('Week')
    plt.ylabel('Number of Importations')
    #on a second y axis plot the seq_count using a line plot
    plt.twinx()
    plt.plot(df3['week_year'], df3['seq_count'], color='r')
    #make y axis log scale
    plt.yscale('log')
    plt.ylabel('Number of Sequences')
    #label the minor ticks on y
    plt.minorticks_on()
    #turn on number labels of minor ticks on y
    plt.yticks(minor=True)
    plt.title(f'{annotation} Importations in {region} per Week (total {importations})')
    plt.show()
    #print value counts for the number of regions in df3
    print(df2['region'].value_counts())
    #create a new column in df3 that has the sum of sample_count for each week
    #drop the rows with more than 1000 sample_count
    #wa_df3 = wa_df3[wa_df3['sample_count'] < 1000]
    #create scatter plot
    df3.plot(x='seq_count', y='importations', kind='scatter')
    plt.xlabel('Number of Sequences')
    plt.ylabel('Number of Importations')
    plt.title(f'{annotation} Importations in {region} per Week vs Seqeunces (total {importations})')
    plt.show()
    #write out importations


# %%
sample_dates.head()

# %%
df[(df['region'] == "Massachusetts") & (df['annotation_3'] == "B.1.1.7")].sort_values(by='earliest_date')

# %%
threshold_dates={'B.1.1.7':pd.to_datetime('2020-12-01'), 'B.1.617.2':pd.to_datetime('2021-03-01')}

def calc_earliest_date(sample_string,variant):
    #for each row split the sample string based on "," then "|" take the last value, to get an array of dates from each sample string
    sample_dates = pd.DataFrame({'sample_dates':[pd.to_datetime(i.split("|")[-1]) for i in sample_string.split(',')]})
    #get the sample_dates that are after the threshold date for B.1.1.7
    earliest_date=sample_dates[sample_dates['sample_dates'] > threshold_dates[variant]].min()
    return earliest_date

# %%


def calc_schedule(df, region, annotations):
    df2 = df[df['region'] == region]
    df2 = df2[df2['annotation_3'].isin(annotations)].sort_values(by='earliest_date')
    df2['earliest_date'] = df2.apply(lambda x: calc_earliest_date(x['samples'], x['annotation_3']) if x['earliest_date'] < threshold_dates[x['annotation_3']] else x['earliest_date'], axis=1)  #get min earliest date
    min_date = df2['earliest_date'].min()
    #drop rows with earliest_date that are na
    df2 = df2.dropna(subset=['earliest_date'])
    #creat tick column by converting earliest_date to datetime and subtracting min_date to create integer
    df2['tick'] = ((pd.to_datetime(df2['earliest_date']) - min_date).dt.days)
    df2['tick'] = df2['tick'].astype(int)
    df3 = df2.groupby(['tick','earliest_date','annotation_3']).size().reset_index(name='importations')
    #df3 rename earliest_date to date and annotation_3 to variant
    df3 = df3.rename(columns={'earliest_date':'date','annotation_3':'variant'})
    df3['seq_count']=df2.groupby(['tick'])['sample_count'].sum().reset_index(name='sample_count')['sample_count']
    #write schedule out to CSV
    df3.to_csv(f'../data/importations/schedules/{region}_schedule.csv', index=False)


# %%
annotations= ['B.1.617.2', 'B.1.1.7']
locations=["Washington","Virginia", "Minnesota", "Massachusetts","Georgia"]
for l in locations:
    calc_schedule(df, l, annotations)

# %% [raw] slideshow={"slide_type": "notes"} vscode={"languageId": "raw"}
# import math
# #for Region Washington and annotation_2 = B.1.617.2 plot the number of rows per week
# def calc_importations(df, regions, annotations):# Create subplots with 4 columns and n rows
#     num_regions = len(regions)
#     num_rows = 1
#     num_col = 2
#     if num_regions > num_col:
#         num_rows = math.ceil(num_regions / num_col)
#     
#
#     for annotation in annotations:
#         fig, axes = plt.subplots(nrows=num_rows, ncols=num_col, figsize=(20, 5*num_col))
#         # Iterate over each region
#         for i, region in enumerate(regions):
#             df2 = df[(df['region'] == region) & (df['annotation_3'] == annotation)]
#             df3 = df2.groupby(['week_year']).size().reset_index(name='importations')
#             df3['seq_count']=df2.groupby(['week_year'])['sample_count'].sum().reset_index(name='sample_count')['sample_count']
#             importations = df3['importations'].sum()
#             df3 = df3.sort_values(by='week_year')
#             cur_row=math.floor(i/num_col)
#             cur_col=i%num_col
#             # Plot "Importations per week" in the first column
#             axes[cur_row, cur_col].bar(df3['week_year'], df3['importations'])
#             axes[cur_row, cur_col].set_xlabel('Week')
#             axes[cur_row, cur_col].set_ylabel('Number of Importations')
#             axes[cur_row, cur_col].set_title(f'Importations per Week in {region}')
#             
#
#
#             # Adjust the spacing between subplots
#             plt.tight_layout()
#
#             # Show the plots
#             plt.show()

# %% [raw] vscode={"languageId": "raw"}
# annotations= ['B.1.617.2', 'B.1.1.7']
# #Massachusetts, Washington, Virginia, Minnesota, Georgia
# locations=["Washington","Virginia", "Massachusetts","Minnesota", "Georgia"]
# calc_importations(df, locations, annotations)

# %%
#plot the total number of sequences per week for all data
df3 = df.groupby(['week_year']).size().reset_index(name='count')
df3['sample_count']=df.groupby(['week_year'])['sample_count'].sum().reset_index(name='sample_count')['sample_count']
df3.plot(x='week_year', y='sample_count', kind='bar')
plt.xlabel('Week')
plt.ylabel('Number of Sequences')
plt.title('Sequences per Week')
#make plot wider
plt.gcf().set_size_inches(20, 8)
#make x-axis labels smaller font
plt.xticks(fontsize=8)
plt.show()


# %%
annotations= ['B.1.617.2', 'B.1.1.7']
locations=["Washington","Virginia", "Minnesota", "Massachusetts","Georgia"]
for a in annotations:
    for l in locations:
        plot_importations(df, l, a)
