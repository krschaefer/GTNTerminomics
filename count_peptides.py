#########################################################################################
# Get how frequently each peptide (of a certain similarity) appears
# Inputs:
#     Excel file: Information about peptides and their start sites
#
# Outputs:
#     Excel file: Frequency of in which peptides appear
#########################################################################################

import pandas as pd
import re
import numpy as np

file = ''

excel = pd.ExcelFile(file)
output_name = file[:-5] + '_pep_freq.xlsx'
output_writer = pd.ExcelWriter(output_name)

# Get names of all the sheets in the excel file
sheet_name_list = excel.sheet_names

# Combine all sheets into a single df
df_ls = []
for i in range(len(sheet_name_list)):
	# Get current sheet
	curr_df = excel.parse(i)
	df_ls.append(curr_df)
	
df = pd.concat(df_ls).reset_index()
df.drop_duplicates(inplace = True)

# Keep only necessary columns
to_keep = ['Accession', 'Prev_AA', "P4-P4'", 
	'Peptide_start', 'Peptide_end', 'Protein_type']
df = df[to_keep]

# Clean up peptide starts in protein repeats
# Only keep first peptide start
clean_up = df[df['Peptide_start'].apply(lambda x: isinstance(x, str))].index
for i in clean_up:
	start_ls = df.iloc[i]['Peptide_start'].split(', ')
	df.iloc[i]['Peptide_start'] = int(start_ls[0])

# Group by protein
prot_gps = df.groupby('Accession')
sliding_dict = {}
SLIDING_WINDOW_SIZE = 3

# Get counts based on a sliding window
for prot_gp in prot_gps:
	curr = prot_gp[1]
	
	# If more than 1 peptide ID'd for this protein
	if len(curr) > 1:
		curr = curr.sort_values(by = 'Peptide_start')
		pep_to_compare = -1
		pep_start_to_compare = -1
		for ind, row in curr.iterrows():
			# If first peptide in protein, set this as first peptide to compare to
			if pep_start_to_compare < 0:
				pep_to_compare = ind
				pep_start_to_compare = row['Peptide_start']
				sliding_dict[pep_to_compare] = 1
			
			# Else, start comparing
			else:
				curr_pep = ind
				curr_pep_start = row['Peptide_start']
				
				# If current peptide start is within SLIDING_WINDOW_SIZE of last peptide
				# Add to the count of the last peptide
				# Consolidate and set count of current peptide to 0
				# Set current peptide start as new start to compare next peptide to
				if curr_pep_start - pep_start_to_compare <= SLIDING_WINDOW_SIZE:
					sliding_dict[pep_to_compare] += 1
					sliding_dict[curr_pep] = 0
					pep_start_to_compare = curr_pep_start
				# If current peptide start is much further down the protein
				# Set current peptide start as new start to compare next peptide to
				# Set current peptide as new peptide to compare to
				else:
					pep_start_to_compare = curr_pep_start
					pep_to_compare = curr_pep
					sliding_dict[curr_pep] = 1
	
	# If only 1 peptide ID'd for this protein
	else:
		ind = curr.index.values[0]
		sliding_dict[ind] = 1

fixed_dict = {}
FIXED_WINDOW_SIZE = 3
# Get counts based on a fixed window
for prot_gp in prot_gps:
	curr = prot_gp[1]
	
	# If more than 1 peptide ID'd for this protein
	if len(curr) > 1:
		curr = curr.sort_values(by = 'Peptide_start')
		curr_inds = curr.index.values
		curr_starts = curr['Peptide_start'].values
		
		# Group peptides by windows
		curr_window = [-1, -1]
		window_dict = {}
		for i in range(len(curr_starts)):
			s = curr_starts[i]
			# If current start is within the current window, keep adding to this group
			if s <= curr_window[1]:
				window_dict[curr_window[0]].append(i)
			# Else, make a new group
			else:
				curr_window = [s, s+(2*FIXED_WINDOW_SIZE)]
				window_dict[s] = [i]				
		# Go through the groups and pick the median as the peptide to add to
		# Set all other peptides to 0
		for k, v in window_dict.items():
			mid = k+FIXED_WINDOW_SIZE
			picked = -1
			diff = FIXED_WINDOW_SIZE*4
			for i in v:
				if picked < 0:
					picked = curr_inds[i]
					diff = abs(curr_starts[i] - mid)
				elif (curr_starts[i] - mid) < diff:
					diff = abs(curr_starts[i] - mid)
					fixed_dict[picked] = 0
					picked = curr_inds[i]
				else:
					fixed_dict[curr_inds[i]] = 0
			fixed_dict[picked] = len(v)		
				
	# If only 1 peptide ID'd for this protein
	else:
		ind = curr.index.values[0]
		fixed_dict[ind] = 1

df['Count_sliding'] = df.index.map(sliding_dict)
df['Count_fixed'] = df.index.map(fixed_dict)

df.sort_values(by = ['Accession', 'Peptide_start'], inplace = True)
#print(df[['Peptide_unmod', 'Peptide_start', 'Count_fixed']].head(20))

df.to_excel(output_writer, index = False)
	
output_writer.save()		
