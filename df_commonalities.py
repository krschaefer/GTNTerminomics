#########################################################################################
# Find common sequences in each excel sheet
# Retain desirable information
# Inputs:
#     Excel file 1: Cell surface N-termini datasets
#     Excel file 2: Domain, TM, Structure, etc. info
#			(Only retain info where if peptide sequences in both files are the same)
#           (Only retain info from Excel file 2)
# Outputs:
#     Excel file: Domain, TM, Structure, etc. info where sequence info matches with
#                 that in excel file 1
#########################################################################################

import pandas as pd
import numpy as np

curated = ''
features = ''

# If there are multiple excel sheets in the features excel file, collapse into one df
def collapse_features(feat_excel):
	# Get names of all the sheets in the excel file
	sheet_name_list = feat_excel.sheet_names
	print(sheet_name_list)
	
	df_list = []
	# Go through each sheet and add to df_list
	for i in range(len(sheet_name_list)):
		# Get values from current sheet
		sheet_df = feat_excel.parse(i)
		df_list.append(sheet_df)
	
	# Collapse into a single df
	all_features = pd.concat(df_list, axis = 0)
	return all_features

# Go through each sheet of the curated excel file
# Determine which sequences match to the features df
# Get corresponding data from features
def get_common(curat_excel, feature_df, output_name = 'curated.xlsx'):
	# New excel file writer
	commonalities_writer = pd.ExcelWriter(output_name)
	
	# Get names of all the sheets in the excel file
	sheet_name_list = curat_excel.sheet_names
	feature_cols = feature_df.columns
	
	# Go through each sheet and find commonalities
	for i in range(len(sheet_name_list)):
		sheet_df = curat_excel.parse(i)
		merged_df = pd.merge(sheet_df, feature_df,  how = 'inner', on = ['P4-P4\'', 'Accession'])
		merged_df[feature_cols].to_excel(commonalities_writer, index = False, sheet_name = sheet_name_list[i])
	
	commonalities_writer.save()
		
curated_excel = pd.ExcelFile(curated)
features_excel = pd.ExcelFile(features)
features_df = collapse_features(features_excel)

df_name = features[:-5] + '_curated.xlsx'
get_common(curated_excel, features_df, df_name)
