#########################################################################################
# Get how frequently each domain is cut
# Inputs:
#     Excel file: Information about which domain is cut in each protein
#
# Outputs:
#     Excel file: Frequency of in which each domain type is cut
#########################################################################################

import pandas as pd
import re
import numpy as np

file = ''

excel = pd.ExcelFile(file)
output_name = file[:-5] + '_freq.xlsx'
output_writer = pd.ExcelWriter(output_name)

# Get names of all the sheets in the excel file
sheet_name_list = excel.sheet_names

for i in range(len(sheet_name_list)):
	# Get current sheet
	curr_df = excel.parse(i)
	
	# Remove rows with no domain information
	curr_df.dropna(subset = ['Domain_range'], inplace = True)
	
	# Remove 'duplicates' (where cut is the same but the peptide differs)
	#      (or where cut is almost the same)
	curr_grouping = curr_df.groupby(['Accession', 'Cut_domain'])
	
	# Go through each group, and decide if it is a 'unique' cut
	cut_df_ls = []
	for gp in curr_grouping:
		curr_gp = gp[1].sort_values(by = ['Domain_dist', 'Intradomain_dist'], ascending = True)
		if curr_gp.shape[0] < 2:
			cut_df_ls.append(curr_gp.index[0])
		else:
			keep = []
			count = 0
			prev_inter_dist = -1
			prev_intra_dist = -1
			for ind, row in curr_gp.iterrows():
				if count == 0:
					keep.append(ind)				# Keep first row, it is unique
					prev_inter_dist = row['Domain_dist']
					prev_intra_dist = row['Intradomain_dist']
				else:
					curr_inter_dist = row['Domain_dist']
					# If interdomain distance is the same
					if curr_inter_dist == prev_inter_dist:
						curr_intra_dist = row['Intradomain_dist']
						# If intradomain_dist is more than 5 AA away, then it is unique
						if (curr_intra_dist > (prev_intra_dist + 5)) or \
						(curr_intra_dist < (prev_intra_dist - 5)):
							keep.append(ind)
							prev_intra_dist = curr_inter_dist
					else:
						# If interdomain_dist is more than 5 AA away, then it is unique
						if (curr_inter_dist > (prev_inter_dist + 5)) or \
						(curr_inter_dist < (prev_inter_dist - 5)):
							keep.append(ind)
							prev_inter_dist = curr_inter_dist
				count += 1
			cut_df_ls.extend(keep)
	
	unique_cut_df = curr_df.loc[cut_df_ls].sort_values(by = ['Cut_domain'], ascending = True)
	
	# Remove domain number from cut_domain
	unique_cut_df['Cut_domain'] = unique_cut_df['Cut_domain'].apply(lambda x: re.sub(' \d', '', x))
	
	# Group by cut_domain and count
	unique_domain_counts = unique_cut_df.groupby(['Cut_domain']).size().reset_index(name='counts')
	unique_domain_total = unique_domain_counts['counts'].sum()
	
	# Get frequency
	unique_domain_counts['freq'] = unique_domain_counts['counts']/unique_domain_total
	
	# Sort by frequency
	unique_domain_counts.sort_values(by = ['freq'], ascending = False, inplace = True)
	
	unique_domain_counts.to_excel(output_writer, index = False, sheet_name = sheet_name_list[i])
	
output_writer.save()	
