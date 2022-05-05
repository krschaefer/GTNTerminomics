import pandas as pd
import numpy as np
import re
from scipy import stats

wollsch = 'wollscheid.csv'
wollsch_HEK = 'wollscheid_HEK.csv'
surfaceomics_file = ''
PEAKS_file = ''
Uniprot_file = 'topo_ls_3.csv'

uniprot_df = pd.read_csv(Uniprot_file)
uniprot_cols_to_keep = ['Entry', 'Mass', 'Sequence', 'Glycosylation', \
	'Domain [FT]', 'Lipidation', 'Topological domain', \
	'Signal peptide', 'Propeptide', 'Subcellular location [CC]', \
	'Beta strand', 'Helix', 'Turn', 'Gene ontology (molecular function)']
uniprot_df = uniprot_df[uniprot_cols_to_keep]

# Pick the columns from PEAKS_file to keep
PEAKS_df = pd.read_csv(PEAKS_file)

PEAKS_df_cols = PEAKS_df.columns
PEAKS_cols_to_keep = ['Accession', 'Peptide', 'PTM']
for c in PEAKS_df_cols:
	if 'Mass' in c:
		PEAKS_cols_to_keep.append(c)
	elif 'Ratio' in c:
		PEAKS_cols_to_keep.append(c)

PEAKS_df = PEAKS_df[PEAKS_cols_to_keep]

# Keep peptides that are abu tagged
PEAKS_df['PTM'] = PEAKS_df['PTM'].fillna('No PTM')

abu_tag_name = 'Abu'
if PEAKS_df['PTM'].str.contains('Amino-butyric Acid').any():
	abu_tag_name = 'Amino-butyric Acid'
elif PEAKS_df['PTM'].str.contains('Amino-butyric Acid').any():
	abu_tag_name = 'Amino-butyric Acid'
	
PEAKS_abu_df = PEAKS_df[PEAKS_df['PTM'].str.contains(abu_tag_name)].copy()

# Get unmodified peptide sequence
def get_unmod_pep(peptide):
	if type(peptide) == float:
		return
	unmod = ''.join(filter(str.isalpha, peptide))
	#if not peptide[1].isalpha():
	if peptide[1] == '.':
		#print(peptide)
		unmod = unmod[1:]
	return unmod
	
PEAKS_abu_df['Peptide_unmod'] = PEAKS_abu_df['Peptide'].apply(get_unmod_pep)

# Get Uniprot Entry
PEAKS_abu_df['Entry'] = PEAKS_abu_df['Accession'].apply(lambda x: x.split('|')[0])

# Rename Mass to peptide mass
if 'Mass' in PEAKS_abu_df.columns:
	PEAKS_abu_df.rename(columns = {'Mass': 'Peptide_mass'}, inplace = True)

# Get number of glycosites in the protein
def get_num_glycosites(entry):
	if type(entry) == str:
		glyco_regex = re.compile(r'CARBOHYD \d+')
		ls = glyco_regex.findall(entry)
		return len(ls)
	else:
		return 0

# Get AA start site of peptide
def get_peptide_start(protein, peptide):
	if type(protein) != str:
		return
	ls = protein.split(peptide)
	if len(ls) > 2:
		starts = []
		tmp_sum = 0
		for i in range(len(ls) - 1):
			tmp_sum += len(ls[i])
			starts.append(str(tmp_sum + 1 + (i * len(peptide))))
		return ', '.join(starts)
		
	else:
		return len(ls[0]) + 1

# Get AA end site of peptide
def get_peptide_end(protein, peptide):
	if type(protein) != str:
		return
	ls = protein.split(peptide)
	if len(ls) > 2:
		ends = []
		tmp_sum = 0
		for i in range(len(ls) - 1):
			tmp_sum += len(ls[i])
			ends.append(str(tmp_sum + ((i + 1) * len(peptide))))
		return ', '.join(ends)
	else:
		return len(ls[0]) + len(peptide)

# Get number of glycosites in each peptide
def get_pep_glyc(glyc, start, end):
	if type(glyc) != str:
		return 0
	glyco_regex = re.compile(r'CARBOHYD \d+')
	ls = glyco_regex.findall(glyc)
	
	if type(start) == str:
		g_count_ls = []
		starts = start.split(', ')
		ends = end.split(', ')
		for s in range(len(starts)):
			g_count = 0
			for i in ls:
				i_num = re.sub(r'CARBOHYD ', '', i)
				if int(i_num) >= int(starts[s]) and int(i_num) <= int(ends[s]):
					g_count += 1
			g_count_ls.append(g_count)
		return ', '.join(str(g_count_ls))
		
	g_count = 0
	for i in ls:
		i_num = re.sub(r'CARBOHYD ', '', i)
		if int(i_num) >= start and int(i_num) <= end:
			g_count += 1
	
	return g_count

# Get the 4 AA prior to peptide
def get_prev_AA(protein, peptide):
	if type(protein) != str:
		return
	ls = protein.split(peptide)
	if len(ls) > 2:
		AA_ls = []
		for i in ls:
			AA_prior = i[-4:]
			AA_ls.append(AA_prior)
		return ', '.join(AA_ls)
	else:
		return ls[0][-4:]

# Get which topological domain was cut
def get_topo_cut(topo, cut_site):
	if type(topo) != str:
		return
	if cut_site == None:
		return
	# Get all topo and corresponding AA ranges
	topo_regex = re.compile(r'(?:\d+\.+)?\d+;\s+/note=".*?"')
	topo_ls = topo_regex.findall(topo)
			
	if cut_site == str:
		cut_ls = cut_site.split(', ')
		cut_domain_ls = []
		for c in cut_ls:
			prev_dom = ''
			prev_end = -1
			found = False
			for t in topo_ls:
				num_regex = re.compile(r'\d+')
				num_ls = num_regex.findall(t)
				
				domain_regex = re.compile(r'note=".*?"')
				topo_domain = domain_regex.findall(t)
				
				# Determine if cut_site is in this AA_range
				if len(num_ls) > 1:
					if int(num_ls[0]) <= int(c) and int(num_ls[1]) >= int(c):
						cut_domain_ls.append(topo_domain[0][6:-1])
						found = True
						break
					else:
						if (prev_end != -1) and (topo_domain[0][6:-1] != prev_dom):
							if prev_end < cut_site and int(num_ls[0]) > cut_site:
								cut_domain_ls.append('Transmembrane')
								found = True
								break
						prev_end = int(num_ls[1])
						prev_dom = topo_domain[0][6:-1]
				else:
					if int(num_ls[0]) == int(c):
						cut_domain_ls.append(topo_domain[0][6:-1])
						found = True
						break
							
			if not found:
				cut_domain_ls.append('unknown')
		return ', '.join(cut_domain_ls)

	else:
		prev_dom = ''
		prev_end = -1
		for t in topo_ls:
			num_regex = re.compile(r'\d+')
			num_ls = num_regex.findall(t)
			
			domain_regex = re.compile(r'note=".*?"')
			topo_domain = domain_regex.findall(t)
			
			# Determine AA range of topo domain
			if len(num_ls) > 1:
				# Determine if cut_site is in this AA_range
				if int(num_ls[0]) <= cut_site and int(num_ls[1]) >= cut_site:
					return topo_domain[0][6:-1]
				# if not in range
				else:
					# if this is not the first topo domain listed and
					# if curr domain is different from prev domain, assume TM
					if (prev_end != -1) and (topo_domain[0][6:-1] != prev_dom):
						if prev_end < cut_site and int(num_ls[0]) > cut_site:
							return 'Transmembrane'
					prev_end = int(num_ls[1])
					prev_dom = topo_domain[0][6:-1]
					
			# If topo domain only consists of a single AA residue			
			else:
				if int(num_ls[0]) == cut_site:
					return topo_domain[0][6:-1]
	return

# Get type of membrane protein
def get_prot_type(entry):
	if type(entry) == str:
		e = entry.lower()
		type_regex = re.compile(';.*?pass .*?membrane protein')
		prot_types = type_regex.findall(e)
		
		if len(prot_types) > 0:
			for i in prot_types:
				i_split = i.split('; ')
				for s in i_split:
					if 'single' in s or 'multi' in s:
						return s
		
		if 'gpi' in e:
			return 'gpi-anchored'
		
		if 'secret' in e:
			return 'secreted'
	return
	
# Get subcellular localization (prettify)
def get_prot_loc(entry):
	if type(entry) != str:
		return 'unknown'
	e = entry.lower()
	e = re.sub(r'([\(\[\{]).*?([\)\]\}])', '', e)
	not_mem_set = set(['mitochondria', 'mitochondrion', \
		'golgi', 'endoplasmic', \
		'nucleolus', 'nuclear', 'nucleus', 'nucleoplasm', 'chromosome', 'dna', \
		'lysosome', 'endosome', 'peroxisome', 'melanosome', \
		'secreted', 'secrete', 'projection', \
		'junction', 'rafts', 'synapse', \
		'cytoplasm', 'cytosol', 'cytoskeleton'])
	mem_set = set(['peripheral', 'plasmamembrane', 'cellmembrane', \
	'mitochondria', 'mitochondrion', \
		'golgi', 'endoplasmic', 'nucleolus', 'nuclear', 'nucleus', \
		'vesicle', 'secreted', 'secrete', 'phagosome', \
		'lysosome', 'endosome', 'peroxisome', 'melanosome', \
		'junction', 'rafts', 'synapse', \
		'cytoplasm', 'cytosol'])
	PM_set = set(['peripheral', 'plasmamembrane', 'cellmembrane', 'junction', 'rafts', 'synapse'])
	mito_set = set(['mitochondria', 'mitochondrion'])
	sec_set = set(['secreted', 'secrete'])
	nuc_set = set(['nucleolus', 'nuclear', 'nucleus', 'nucleoplasm', 'chromosome', 'dna'])
	cyto_set = set(['cytoplasm', 'cytosol'])
	
	if 'membrane' in e:
		e = re.sub('plasma mem', 'plasmamem', e)
		e = re.sub('cell mem', 'cellmem', e)
		e_set = set(re.findall(r'[a-z]+', e))
		curr_set = mem_set.intersection(e_set)
		if len(curr_set) < 1:
			return 'membrane'
		else:
			in_PM = False
			in_mito = False
			in_sec = False
			in_nuc = False
			in_cyto = False
			if len(curr_set.intersection(PM_set)) > 0:
				in_PM = True
				curr_set = curr_set.difference(PM_set)
			if len(curr_set.intersection(mito_set)) > 0:
				in_mito = True
				curr_set = curr_set.difference(mito_set)
			if len(curr_set.intersection(sec_set)) > 0:
				in_sec = True
				curr_set = curr_set.difference(sec_set)
			if len(curr_set.intersection(nuc_set)) > 0:
				in_nuc = True
				curr_set = curr_set.difference(nuc_set)
			if len(curr_set.intersection(cyto_set)) > 0:
				in_cyto = True
				curr_set = curr_set.difference(cyto_set)
			
			all_loc = list(curr_set)
			if in_PM:
				all_loc.append('PM')
			if in_mito:
				all_loc.append('mitochondria')
			if in_sec:
				all_loc.append('secreted')
			if in_nuc:
				all_loc.append('nucleus')
			if in_cyto:
				all_loc.append('cytoplasm')
			locations = ', '.join(all_loc)
			locations = locations.replace('endoplasmic', 'ER')
			return locations
	else:
		e_set = set(re.findall(r'[a-z]+', e))
		curr_set = not_mem_set.intersection(e_set)
		if len(curr_set) < 1:
			return 'unknown'
		else:
			in_PM = False
			in_mito = False
			in_sec = False
			in_nuc = False
			in_cyto = False
			if len(curr_set.intersection(PM_set)) > 0:
				in_PM = True
				curr_set = curr_set.difference(PM_set)
			if len(curr_set.intersection(mito_set)) > 0:
				in_mito = True
				curr_set = curr_set.difference(mito_set)
			if len(curr_set.intersection(sec_set)) > 0:
				in_sec = True
				curr_set = curr_set.difference(sec_set)
			if len(curr_set.intersection(nuc_set)) > 0:
				in_nuc = True
				curr_set = curr_set.difference(nuc_set)
			if len(curr_set.intersection(cyto_set)) > 0:
				in_cyto = True
				curr_set = curr_set.difference(cyto_set)
			
			all_loc = list(curr_set)
			if in_PM:
				all_loc.append('PM')
			if in_mito:
				all_loc.append('mitochondria')
			if in_sec:
				all_loc.append('secreted')
			if in_nuc:
				all_loc.append('nucleus')
			if in_cyto:
				all_loc.append('cytoplasm')
			locations = ', '.join(all_loc)
			locations = locations.replace('endoplasmic', 'ER')
			return locations
	
# Get distance from signal sequence
# What are we calling "distance"? 
# Is the first AA after the signal sequence considered 1 away?
def dist_from_SS(SS, start):
	if type(SS) != str:
		return 'no SS'
	if start == None:
		return
	SS_regex = re.compile('SIGNAL \d+\.\.\d+')
	SS_range = SS_regex.findall(SS)
	if type(start) == str:
		starts = start.split(', ')
		dists = []
		for i in starts:
			SS_len = SS_range[0].split('..')[-1]
			dist = int(i) - int(SS_len) + 1
			dists.append(str(dist))
		
		return ', '.join(dists)
	
	if len(SS_range) > 0:
		SS_len = SS_range[0].split('..')[-1]
		dist = start - int(SS_len) + 1
		return dist
	else:
		return
	
# Determine if cut is in SS
def SS_in_topo(cut_topo, ss_dist):
	if type(ss_dist) != str and ss_dist != None:
		if ss_dist <= 0:
			if type(cut_topo) != str:
				return 'SS'
			else:
				return cut_topo + ', SS'

		elif ss_dist < 5:
			if type(cut_topo) != str:
				return 'Near SS'
			else:
				return cut_topo + ', Near SS'
	if type(cut_topo) == str:
		return cut_topo
		
# Get distance from propeptide
def dist_from_pro(propep, start):
	if type(propep) != str:
		return 'no Propeptide'
	if start == None:
		return
	propep_regex = re.compile('PROPEP \d+\.\.\d+')
	propep_range = propep_regex.findall(propep)
	
	if type(start) == str:
		starts = start.split(', ')
		dists = []
		for i in starts:
			if len(propep_range) > 0:
				propep_start, propep_end = propep_range[0].split('..')
				propep_start = propep_start.split(' ')[-1]
				if int(i) < int(propep_start):
					dist = int(i) - int(propep_start)
				elif int(i) > int(propep_end):
					dist = int(i) - int(propep_end)
				else:
					dist = 0
				dists.append(str(dist))
		return ', '.join(dists)
			
	if len(propep_range) > 0:
		propep_start, propep_end = propep_range[0].split('..')
		propep_start = propep_start.split(' ')[-1]
		if start < int(propep_start):
			return start - int(propep_start)
		elif start > int(propep_end):
			return start - int(propep_end)
		else:
			return 0
	else:
		return

# Determine if cut is in propeptide
############ Include type(start) == str???
def cut_in_pro(cut_topo, propep, start):
	if type(propep) != str or start == None:
		return cut_topo
	propep_regex = re.compile('PROPEP \d+\.\.\d+')
	propep_range = propep_regex.findall(propep)
	
	if type(start) == str:
		return
	if len(propep_range) > 0:
		propep_range = re.sub('PROPEP ', '', propep_range[0])
		propep_start, propep_end = propep_range.split('..')
		if start >= int(propep_start) and start <= int(propep_end):
			if type(cut_topo) == str:
				return cut_topo + ', Propeptide'
			else:
				return 'Propeptide'
		else:
			if type(cut_topo) == str:
				return cut_topo
	else:
		if type(cut_topo) == str:
			return cut_topo
	if type(cut_topo) == str:
		return cut_topo

# determine if protein is in wollscheid dataset
def in_wollsch(prot, wollsch):
	if prot in wollsch:
		return 'yes'
	else:
		return

# Get structural element that is cut
def cut_struct(beta, alpha, turn, cut_site):
	if (type(beta) != str) and (type(alpha) != str) and (type(turn) != str):
		return None
	if cut_site == None:
		return None
	if type(cut_site) == str:
		return None
	
	struct_dict = {}
	struct_regex = re.compile(r'(?:\d+\.+)?\d+;')
	# Get all beta strand regions
	if type(beta) == str:
		beta_ls = struct_regex.findall(beta)
		for i in beta_ls:
			struct_dict[i.strip(';')] = 'beta'
	# Get all alpha helix regions
	if type(alpha) == str:
		alpha_ls = struct_regex.findall(alpha)
		for i in alpha_ls:
			struct_dict[i.strip(';')] = 'alpha'
	# Get all turn regions
	if type(turn) == str:
		turn_ls = struct_regex.findall(turn)
		for i in turn_ls:
			struct_dict[i.strip(';')] = 'turn'
	
	# Deduce unstructured regions
	# And store all structural data in a new dict
	full_struct_dict = {}
	previous = -1
	for key in sorted(struct_dict):
		 if '..' in key:
		 	start, end = key.split('..')
		 else:
		 	start, end = key, key
		 if previous < 0:
		 	full_struct_dict[key] = struct_dict[key]
		 else:
		 	if int(start) > (previous + 1):
		 		new_range = '%d..%d'%(previous+1, int(start)-1)
		 		full_struct_dict[new_range] = 'unstructured'
		 		full_struct_dict[key] = struct_dict[key]
		 		
		 	else:
		 		full_struct_dict[key] = struct_dict[key]
		 previous = int(end)
	
	# Go through and find out what struct is cut
	for key in sorted(full_struct_dict):
		if '..' in key:
			start, end = key.split('..')
		else:
			start, end = key, key
		
		if cut_site >= int(start) and cut_site <= int(end):
			return full_struct_dict[key]
	
	return 'outside struct'
		 
# Get which domain was cut
def cut_domain(domains, cut_site):
	if type(domains) != str:
		return [None, None]
	if cut_site == None:
		return [None, None]
	if type(cut_site) == str:
		return [None, None]
	
	# Get all domains and corresponding AA ranges
	domain_regex = re.compile(r'(?:\d+\.+)?\d+;\s+/note=".*?"')
	domain_ls = domain_regex.findall(domains)
	
	prev = -1
	prev_range = ''
	prev_dom = ''
	for d in domain_ls:
		d_range, d_type = d.split(';  /note=')			# Get domain range and domain name
		if '..' in d_range:								# Get start and end of domain
			d_start, d_end = d_range.split('..')
		else:
			d_start = d_range
			d_end = d_range
		
		if prev == -1:									# If we are at the start
			if cut_site < int(d_start):					# If cut is before first domain
				d_range = 'before %s'%d_range
				return [d_range, 'before first domain']	# protein is whole
			elif cut_site >= int(d_start) and cut_site <= int(d_end):		# If cut is in first domain
				return [d_range, d_type]									# Return this
			
			# Keep track of this first domain
			prev = int(d_end)
			prev_range = d_range
			prev_dom = d_type
		
		else:
			if cut_site >= int(d_start) and cut_site <= int(d_end):		# If cut is in domain
				return [d_range, d_type]
			elif cut_site > prev and cut_site <= int(d_start):			# If cut is between 2 domains
				btwn_ranges = '%s and %s'%(prev_range, d_range)
				btwn_doms = '%s and %s'%(prev_dom, d_type)
				return [btwn_ranges, btwn_doms]
			else:
				prev = int(d_end)				# Keep track of this domain and
				prev_range = d_range			# Move on to analyze the next one
				prev_dom = d_type
	
	d_range = 'after %s'%prev_range
	return [d_range, 'after last domain']

# Determine which domains were removed by the cut
def removed_domains(domains, dom_range, prot_type):
	if type(domains) != str:
		return [None, None]
	if type(dom_range) != str:
		return [None, None]
	if ('before' in dom_range) or ('after' in dom_range):
		return [None, None] 
	if type(prot_type) != str:
		return [None, None]

	# Get all domains and corresponding AA ranges
	domain_regex = re.compile(r'(?:\d+\.+)?\d+;\s+/note=".*?"')
	domain_ls = domain_regex.findall(domains)
	
	#Determine if ECD is N-terminal or C-terminal
	ECD = ''
	if ('type i ' in prot_type) or ('type iii ' in prot_type) or ('gpi' in prot_type):
		ECD = 'N'
	if ('type ii ' in prot_type) or ('type iv ' in prot_type):
		ECD = 'C'
	
	domain_A = ''
	if ' and ' in dom_range:
		domain_A = dom_range.split(' and ')[0]
	else:
		domain_A = dom_range
		
	# Determine which domains are removed
	removed = []
	found = False
	for i in domain_ls:
		i_type = i.split('note=')[-1]
		if ECD == 'N':
			if domain_A in i:
				removed.append(i_type)
				found = True
				num_removed = len(removed)
				rem = ', '.join(removed)
				return [rem, num_removed]
			else:
				removed.append(i_type)
		elif ECD == 'C':
			if found == True:
				removed.append(i_type)
			elif domain_A in i:
				found = True
				removed.append(i_type)
			else:
				continue
				
	if ECD == 'C':
		num_removed = len(removed)
		rem = ', '.join(removed)
		return [rem, num_removed]
	return [None, None]

# Determine which domains remain after the cut	
def kept_domains(domains, dom_range, prot_type):
	if type(domains) != str:
		return [None, None]
	if type(dom_range) != str:
		return [None, None]
	if ('before' in dom_range) or ('after' in dom_range):
		return [None, None] 
	if type(prot_type) != str:
		return [None, None]

	# Get all domains and corresponding AA ranges
	domain_regex = re.compile(r'(?:\d+\.+)?\d+;\s+/note=".*?"')
	domain_ls = domain_regex.findall(domains)
	
	#Determine if ECD is N-terminal or C-terminal
	ECD = ''
	if ('type i ' in prot_type) or ('type iii ' in prot_type) or ('gpi' in prot_type):
		ECD = 'N'
	if ('type ii ' in prot_type) or ('type iv ' in prot_type):
		ECD = 'C'
	
	domain_A = ''
	if ' and ' in dom_range:
		domain_A = dom_range.split(' and ')[0]
	else:
		domain_A = dom_range
		
	# Determine which domains are removed
	kept = []
	found = False
	for i in domain_ls:
		i_type = i.split('note=')[-1]
		if ECD == 'C':
			if domain_A in i:
				found = True
				num_kept = len(kept)
				ke = ', '.join(kept)
				return [ke, num_kept]
			else:
				kept.append(i_type)
		elif ECD == 'N':
			if found == True:
				kept.append(i_type)
			elif domain_A in i:
				found = True
			else:
				continue
				
	if ECD == 'N':
		num_kept = len(kept)
		ke = ', '.join(kept)
		return [ke, num_kept]
	return [None, None]
	
# Determine the distance between domain and ID'd cut site
def nearest_domain(dom_range, cut_site):
	if type(dom_range) != str:
		return None
	if cut_site == None:
		return None
	if type(cut_site) == str:
		return None
	
	# If cut is before all domains
	if 'before' in dom_range:
		dom_start = dom_range[7:].split('..')[0]
		dist = cut_site - int(dom_start)
		return -dist
	
	# If cut is after all domains
	elif 'after' in dom_range:
		dom_end = dom_range.split('..')[-1]
		dist = cut_site - int(dom_end)
		return -dist
	
	# If cut is in a domain, distance is 0
	elif ' and ' not in dom_range:
		return 0
		
	dom_A, dom_B = dom_range.split(' and ')
	if '..' in dom_A:
		dom_A_end = dom_A.split('..')[-1]
	else:
		dom_A_end = dom_A
	if '..' in dom_B:
		dom_B_start = dom_B.split('..')[0]
	else:
		dom_B_start = dom_B
	
	dist_to_A = cut_site - int(dom_A_end)
	dist_to_B = cut_site - int(dom_B_start)
	
	if abs(dist_to_A) <= abs(dist_to_B):
		return -dist_to_A
	else:
		return -dist_to_B

# Determine distance within a domain of in which the ID'd cut site is
def domain_dist(dom_range, cut_site):
	if type(dom_range) != str:
		return None
	if cut_site == None:
		return None
	if type(cut_site) == str:
		return None
	
	# If cut is before all domains
	if 'before' in dom_range:
		return None
	
	# If cut is after all domains
	elif 'after' in dom_range:
		return None
	
	# If cut is in a domain
	elif ' and ' not in dom_range:
		if '..' in dom_range:
			start, end = dom_range.split('..')
		else:
			start = dom_range
		dist = cut_site - int(start)
		return -dist
	
	else:
		return None

# Get the nearest TM to the ID'd cut site
def nearest_TM(topo, cut_site):
	if type(topo) != str:
		return [None, None]
	if cut_site == None:
		return [None, None]
	
	# Get all topo and corresponding AA ranges
	topo_regex = re.compile(r'(?:\d+\.+)?\d+;\s+/note=".*?"')
	topo_ls = topo_regex.findall(topo)
	
	# Derive TMs
	prev = -1
	TM_ls = []
	for t in topo_ls:
		# Get start and end of each topological domain
		t_range = t.split(';')[0]
		if '..' in t_range:
			t_start, t_end = t_range.split('..')
		else:
			t_start = t_range
			t_end = t_range
		# Store the end of the first topological domain (this is start of first TM)
		if prev == -1:
			prev = int(t_end) + 1
		# Record start and end of TMs
		else:
			TM_ls.append('%d..%d'%(prev, int(t_start) - 1))
			prev = int(t_end) + 1
			
	if cut_site == str:
		cut_ls = cut_site.split(', ')
		cut_dist_ls = []
		cut_TM_ls = []
		for c in cut_ls:
			prev_dist = -1
			prev_end = -1
			closest_TM = ''
			found = False
			# Look through each TM
			for i in range(len(TM_ls)):
				t_start, t_end = TM_ls[i].split('..')
				dist = int(c) - int(t_start)			# Get distance from cut site to start of TM
				if i == 0:
					closest_TM = TM_ls[i]				# If this is the first TM
					if int(c) < int(t_start):			# If the cut site is before first TM
						cut_dist_ls.append(str(dist))	# It must be closest to first TM
						cut_TM_ls.append(closest_TM)
						found = True
						break
					else:
						prev_dist = dist				# Else, compare this dist to
						prev_end = int(t_end)			# Distance of next TM
				else:
					dist_to_prev = int(c) - prev_end
					if abs(prev_dist) < abs(dist_to_prev):	# If we are closer to start of the previous TM
						cut_dist_ls.append(str(-prev_dist))	# than to end of previous TM, return this
						cut_TM_ls.append(closest_TM)
						found = True
						break
					elif abs(dist_to_prev) < abs(dist):		# If we are closer to end of the previous TM
						cut_dist_ls.append(str(-dist_to_prev))	# than to start of current TM, return this
						cut_TM_ls.append(closest_TM)
						found = True
						break
					elif i == len(TM_ls) - 1:				# If exhausted all possibities
						cut_dist_ls.append(str(-dist))		# Must be closest to last TM
						cut_TM_ls.append(TM_ls[i])
						found = True
						break
					else:
						closest_TM = TM_ls[i]				# Else, save current TM, and continue
						prev_dist = dist
						prev_end = int(t_end)
			if not found:
				cut_dist_ls.append('unknown')
				cut_TM_ls.append('unknown')
		cut_dists = ', '.join(cut_dist_ls)
		cut_TMs = ', '.join(cut_TM_ls)
		return[cut_dists, cut_TMs]

	else:
		prev_dist = -1
		prev_end = -1
		closest_TM = ''
		# Look through each TM
		for i in range(len(TM_ls)):
			t_start, t_end = TM_ls[i].split('..')
			dist = cut_site - int(t_start)			# Get distance from cut site to start of TM
			if i == 0:
				closest_TM = TM_ls[i]				# If this is the first TM
				if cut_site < int(t_start):			# If the cut site is before first TM
					return [-dist, closest_TM]		# It must be closest to first TM
				else:
					prev_dist = dist				# Else, compare this dist to
					prev_end = int(t_end)			# Distance of next TM
			else:
				dist_to_prev = cut_site - prev_end
				#print(i, cut_site, prev_dist, dist_to_prev, dist)
				if abs(prev_dist) < abs(dist_to_prev):	# If we are closer to start of the previous TM
					#print(prev_dist)
					return [-prev_dist, closest_TM]		# than to end of previous TM, return this
				elif abs(dist_to_prev) < abs(dist):		# If we are closer to end of the previous TM
					#print(dist_to_prev)
					return [-dist_to_prev, closest_TM]	# than to start of current TM, return this
				elif i == len(TM_ls) - 1:				# If exhausted all possibities
					#print(dist)
					return [-dist, TM_ls[i]]				# Must be closest to last TM
				else:
					closest_TM = TM_ls[i]				# Else, save current TM, and continue
					prev_dist = dist
					prev_end = int(t_end)
			
	return [None, None]

# Get the nearest TM to the ID'd cut site
def nearest_lipid(lipid, cut_site):
	if type(lipid) != str:
		return [None, None, None]
	if cut_site == None:
		return [None, None, None]
	
	# Get all topo and corresponding AA ranges
	lipid_regex = re.compile(r'(?:\d+\.+)?\d+;\s+/note=".*?"')
	lipid_ls = lipid_regex.findall(lipid)
	
	if type(cut_site) == str:
		cut_ls = cut_site.split(', ')
		cut_dist_ls = []
		cut_lsite_ls = []
		cut_ltype_ls = []
		for c in cut_ls:
			dist = -1
			abs_dist = -1
			lipid_site = -1
			lipid_type = ''
			for i in lipid_ls:
				i_ls = i.split(';  /note=')
				l_site = int(i_ls[0])
				l_type = i_ls[1].strip('"')
			
				curr_dist = int(c) - l_site
			
				if abs(curr_dist) < abs_dist or abs_dist < 0:
					dist = curr_dist
					lipid_site = l_site
					lipid_type = l_type
					abs_dist = abs(curr_dist)
			cut_dist_ls.append(str(-dist))
			cut_dists = ', '.join(cut_dist_ls)
			cut_lsite_ls.append(str(lipid_site))
			cut_lsites = ', '.join(cut_lsite_ls)
			cut_ltype_ls.append(lipid_type)
			cut_ltypes = ', '.join(cut_ltype_ls)
		return [cut_dists, cut_lsites, cut_ltypes]
		
	else:
		dist = -1
		abs_dist = -1
		lipid_site = -1
		lipid_type = ''
		for i in lipid_ls:
			i_ls = i.split(';  /note=')
			l_site = int(i_ls[0])
			l_type = i_ls[1].strip('"')
			
			curr_dist = cut_site - l_site
			
			if abs(curr_dist) < abs_dist or abs_dist < 0:
				dist = curr_dist
				lipid_site = l_site
				lipid_type = l_type
				abs_dist = abs(curr_dist)
		
		return [-dist, lipid_site, lipid_type]
					
# Get the nearest glycosite to the ID'd cut site
def nearest_glycosite(glycosites, cut):
	if type(cut) == str:
		cuts = cut.split(', ')
		closest_ls = []
		glyco_ls = []
		if type(glycosites) == str:
			glyco_regex = re.compile(r'CARBOHYD \d+')
			ls = glyco_regex.findall(glycosites)
			
			for c in cuts:
				abs_dist = -1
				closest = -1
				closest_glycosite = -1
				for i in ls:
					glyco_num = int(i.split(' ')[1])
					curr_dist = int(c) - glyco_num			# Get dist of glycosite from cut site
					if curr_dist == 0:
						closest_ls.append(str(0))
						glyco_ls.append(c)
						break
			
					if abs_dist == -1:
						abs_dist = abs(curr_dist)
						closest = curr_dist
						closest_glycosite = glyco_num
				
					elif abs(curr_dist) < abs_dist:
						abs_dist = abs(curr_dist)
						closest = curr_dist
						closest_glycosite = glyco_num
				closest_ls.append(str(-closest))
				glyco_ls.append(str(closest_glycosite))
			closest_sites = ', '.join(closest_ls)
			closest_glycos = ', '.join(glyco_ls)
			return[closest_sites, closest_glycos]
		else:
			return [None, None]
			
	if type(glycosites) == str:
		glyco_regex = re.compile(r'CARBOHYD \d+')
		ls = glyco_regex.findall(glycosites)
		abs_dist = -1
		closest = -1
		closest_glycosite = -1
		for i in ls:
			glyco_num = int(i.split(' ')[1])
			curr_dist = cut - glyco_num				# Get dist of glycosite from cut site
			if curr_dist == 0:
				return [0, cut]
			
			if abs_dist == -1:
				abs_dist = abs(curr_dist)
				closest = curr_dist
				closest_glycosite = glyco_num
				
			elif abs(curr_dist) < abs_dist:
				abs_dist = abs(curr_dist)
				closest = curr_dist
				closest_glycosite = glyco_num
		return [-closest, closest_glycosite]
	else:
		return [None, None]

# Get molecular function (based on cytokine/receptor type)
def get_molec_func(entry):
	if type(entry) != str:
		return
	
	if ('cytokine' in entry) or ('receptor' in entry):
		receptor_func = []
		mol_functions = entry.split('; ')
		for i in mol_functions:
			if ('cytokine' in i) or ('receptor' in i):
				f = i.split(' [')[0]
				receptor_func.append(f)
		return '; '.join(receptor_func)
	else:
		return
	
# Get wilcoxon p-value (per peptide)
def get_pval(FC_ls):
	log2_FC_ls = []
	for i in FC_ls:
		if i != '-':
			log2_FC = np.log2(float(i))
			log2_FC_ls.append(log2_FC)	
	try:
		res = stats.wilcoxon(log2_FC_ls)[1]
		pval = -np.log10(res)
	except:
		pval = np.nan
	return pval

# Uncomment to add protein type annotation to wollscheid dataset
"""
wollsch_df = pd.read_csv(wollsch)
wollsch_df = wollsch_df.merge(uniprot_df, how = 'left', left_on = 'Protein', right_on = 'Entry')
wollsch_df['Protein_type'] = wollsch_df['Subcellular location [CC]'].apply(get_prot_type)
wollsch_df.drop(columns = uniprot_cols_to_keep, inplace = True)
wollsch_df.to_csv(wollsch, index = False)
"""

# Get set of proteins from wollscheid CSPA dataset
wollsch_df = pd.read_csv(wollsch)
wollsch_HEK_df = pd.read_csv(wollsch_HEK)
wollsch_set = set(wollsch_df['Protein'].tolist())
wollsch_HEK_set = set(wollsch_HEK_df['ID link'].tolist())

########################################################################################
# Combine PEAKS data and uniprot data
abu_df = PEAKS_abu_df.merge(uniprot_df, how = 'left', on = 'Entry')

########################################################################################
# Get log2(Experimental/Control)
abu_sample_cols = [x for x in abu_df.columns if 'Sample' in x and 'Ratio' in x]
abu_gp_cols = [x for x in abu_df.columns if 'Group' in x]
log2FC_ls = []
for i in abu_sample_cols:
	abu_df[i] = abu_df[i].replace(r'-', np.NaN, regex=True)
	abu_df[i] = abu_df[i].astype(float)
	gp_number = i.split('Sample ')[-1]
	abu_df['log2FC_%s'%gp_number] = np.log2(abu_df[i])
	log2FC_ls.append('log2FC_%s'%gp_number)

abu_df['log2FC_med'] = abu_df[log2FC_ls].median(axis = 1)
log2FC_ls.append('log2FC_med')
#abu_df['pval'] = abu_df.apply(lambda x: get_pval(x[abu_sample_cols]), axis = 1)

########################################################################################
# Grouping and concatenating peptides
abu_df['Peptide_start'] = abu_df.apply(
	lambda x: get_peptide_start(x['Sequence'], x['Peptide_unmod']), axis = 1)
abu_df['Peptide_end'] = abu_df.apply(
	lambda x: get_peptide_end(x['Sequence'], x['Peptide_unmod']), axis = 1)

pep_gps = abu_df.groupby(by = ['Accession', 'Peptide_start'])
df_ls = []
for i in pep_gps:
	if len(i[1]) > 1:
		curr = i[1].copy().sort_values(by = 'Peptide_end').reset_index()
		for FC in log2FC_ls:
			FC_med = curr[FC].median()
			curr.at[0, FC] = FC_med
		df_ls.append(curr.iloc[0].to_frame().transpose().drop(columns = 'index'))
	else:
		df_ls.append(i[1])

abu_df = pd.concat(df_ls).reset_index()

########################################################################################
# Get some info about abu dataset
abu_df['Molec_func'] = abu_df['Gene ontology (molecular function)'].apply(get_molec_func)
abu_df['Num_glycosites'] = abu_df['Glycosylation'].apply(get_num_glycosites)
abu_df['Pep_glycosites'] = abu_df.apply(lambda x: get_pep_glyc(x['Glycosylation'], x['Peptide_start'], x['Peptide_end']), axis = 1)
abu_df['Prev_AA'] = abu_df.apply(lambda x: get_prev_AA(x['Sequence'], x['Peptide_unmod']), axis = 1)
abu_df['P4-P4\''] = abu_df['Prev_AA'] + abu_df['Peptide_unmod'].str[:4]
abu_df['Cut_topo'] = abu_df.apply(lambda x: get_topo_cut(x['Topological domain'], x['Peptide_start']), axis = 1)
abu_df['Protein_type'] = abu_df['Subcellular location [CC]'].apply(get_prot_type)
abu_df['Location'] = abu_df['Subcellular location [CC]'].apply(get_prot_loc)
abu_df['SS_dist'] = abu_df.apply(lambda x: dist_from_SS(x['Signal peptide'], x['Peptide_start']), axis = 1)
abu_df['Cut_topo'] = abu_df.apply(lambda x: SS_in_topo(x['Cut_topo'], x['SS_dist']), axis = 1)
abu_df['Propep_dist'] = abu_df.apply(lambda x: dist_from_pro(x['Propeptide'], x['Peptide_start']), axis = 1)
abu_df['Cut_topo'] = abu_df.apply(lambda x: cut_in_pro(x['Cut_topo'], x['Propeptide'], x['Peptide_start']), axis = 1)
abu_df['Cut_struct'] = abu_df.apply(lambda x: cut_struct(x['Beta strand'], x['Helix'], x['Turn'], x['Peptide_start']), axis = 1)

tmp = abu_df.apply(lambda x: cut_domain(x['Domain [FT]'], x['Peptide_start']), axis = 1)
abu_df['Domain_range'] = tmp.apply(lambda x: x[0])
abu_df['Cut_domain'] = tmp.apply(lambda x: x[1])

tmp = abu_df.apply(lambda x: removed_domains(x['Domain [FT]'], x['Domain_range'], x['Protein_type']), axis = 1)
abu_df['Domains_removed'] = tmp.apply(lambda x: x[0])
abu_df['Num_domains_removed'] = tmp.apply(lambda x: x[1])

tmp = abu_df.apply(lambda x: kept_domains(x['Domain [FT]'], x['Domain_range'], x['Protein_type']), axis = 1)
abu_df['Domains_intact'] = tmp.apply(lambda x: x[0])
abu_df['Num_domains_intact'] = tmp.apply(lambda x: x[1])

abu_df['Domain_dist'] = abu_df.apply(lambda x: nearest_domain(x['Domain_range'], x['Peptide_start']), axis = 1)
abu_df['Intradomain_dist'] = abu_df.apply(lambda x: domain_dist(x['Domain_range'], x['Peptide_start']), axis = 1)

tmp = abu_df.apply(lambda x: nearest_TM(x['Topological domain'], x['Peptide_start']), axis = 1)
abu_df['TM_dist'] = tmp.apply(lambda x: x[0])
abu_df['nearest_TM'] = tmp.apply(lambda x: x[1])

tmp = abu_df.apply(lambda x: nearest_lipid(x['Lipidation'], x['Peptide_start']), axis = 1)
abu_df['lipid_dist'] = tmp.apply(lambda x: x[0])
abu_df['nearest_lipid'] = tmp.apply(lambda x: x[1])
abu_df['nearest_lipid_type'] = tmp.apply(lambda x: x[2])

tmp = abu_df.apply(lambda x: nearest_glycosite(x['Glycosylation'], x['Peptide_start']), axis = 1)
abu_df['glycosite_dist'] = tmp.apply(lambda x: x[0])
abu_df['nearest_glycosite'] = tmp.apply(lambda x: x[1])

abu_df['Wollscheid'] = abu_df.apply(lambda x: in_wollsch(x['Entry'], wollsch_set), axis = 1)
abu_df['Wollscheid_HEK'] = abu_df.apply(lambda x: in_wollsch(x['Entry'], wollsch_HEK_set), axis = 1)
abu_df['P1'] = abu_df['Prev_AA'].str[-1]

########################################################################################
# Comparing to Kevin's surfaceomic's dataset
KL_df = pd.read_csv(KL_file)
KL_cols = ['Background', 'log2 median enrichment', 'log10 pvalue']

#KL_df['Entry'] = KL_df['Accession'].apply(lambda x: x.split('|')[0])
#KL_keep_cols = ['Entry']

# Get log2(Experimental/Control)
#KL_gp_cols = [x for x in KL_df.columns if 'Ratio' in x]
#for i in KL_gp_cols:
#	if 'Group' in i:
#		gp_number = i.split('Group ')[-1]
#		KL_df['KL_log2FC_Gp_%s'%gp_number] = np.log2(KL_df[i])
#		KL_keep_cols.append('KL_log2FC_Gp_%s'%gp_number)
#	elif 'Sample' in i:
#		smpl_number = i.split('Sample ')[-1]
#		KL_df[i] = KL_df[i].replace(r'-', np.NaN, regex=True)
#		KL_df[i] = KL_df[i].astype(float)
#		KL_df['KL_log2FC_%s'%smpl_number] = np.log2(KL_df[i])
#		KL_keep_cols.append('KL_log2FC_%s'%smpl_number)

#KL_sample_cols = [x for x in KL_df.columns if ('Sample' in x) and ('Ratio' in x)]
#KL_df['KL_pval'] = KL_df.apply(lambda x: get_pval(x[KL_sample_cols]), axis = 1)
#KL_keep_cols.append('KL_pval')

#KL_df = KL_df[KL_keep_cols]
#final_df = abu_df.merge(KL_df, how = 'left', on = 'Entry')

final_df = abu_df.merge(KL_df, how = 'left', on = 'Entry')

########################################################################################
# Grouping by residue type and saving file
basic_aa = ['K', 'R']
aliphatic_aa = ['G', 'L', 'A', 'M', 'V', 'I', 'P', 'F', 'Y', 'W']
other_aa = ['S', 'T', 'Q', 'N', 'E', 'D', 'H', 'C']

# Remove unnecessary columns
cols_to_drop = ['Entry', 'Glycosylation', 'Topological domain', 'Sequence', \
	'Propeptide', 'Signal peptide', 'Subcellular location [CC]', 'Lipidation']
final_df = final_df.drop(columns = cols_to_drop)

# Separate df by P1 residue type
basic_df = final_df[final_df['P1'].isin(basic_aa)].copy()
aliphatic_df = final_df[final_df['P1'].isin(aliphatic_aa)].copy()
other_df = final_df[final_df['P1'].isin(other_aa)].copy()

# Get indices for extracellular proteins, gpi/secreted proteins, and intracellular proteins
all_ecd_ind = set(final_df.dropna(subset = ['Protein_type']).index)
sus_ecd_ind = set(final_df[final_df['Protein_type'].isin(['gpi-anchored', 'secreted'])].index)
real_ecd_ind = all_ecd_ind.difference(sus_ecd_ind)
cyto_ind = set(final_df.index).difference(all_ecd_ind)

# Get df with ECD proteins only
ecd_df = final_df.loc[list(real_ecd_ind)]
basic_ecd_df = ecd_df[ecd_df['P1'].isin(basic_aa)].copy()
aliphatic_ecd_df = ecd_df[ecd_df['P1'].isin(aliphatic_aa)].copy()
other_ecd_df = ecd_df[ecd_df['P1'].isin(other_aa)].copy()

# Get df with gpi/secreted proteins only
sus_df = final_df.loc[list(sus_ecd_ind)]
basic_sus_df = sus_df[sus_df['P1'].isin(basic_aa)].copy()
aliphatic_sus_df = sus_df[sus_df['P1'].isin(aliphatic_aa)].copy()
other_sus_df = sus_df[sus_df['P1'].isin(other_aa)].copy()

# Get df with only cytoplasmic proteins
cyto_df = final_df.loc[list(cyto_ind)]
basic_cyto_df = cyto_df[cyto_df['P1'].isin(basic_aa)].copy()
aliphatic_cyto_df = cyto_df[cyto_df['P1'].isin(aliphatic_aa)].copy()
other_cyto_df = cyto_df[cyto_df['P1'].isin(other_aa)].copy()

base_name = PEAKS_file[:-4]+'_v2'
base_cols = ['Peptide_unmod', 'Accession', 'log2FC_med', 'SS_dist', \
	'Prev_AA', 'P4-P4\'', \
	'Peptide_start', 'Peptide_end',	'Cut_topo', 'Num_glycosites', 'Protein_type', 'Location', \
	'Propep_dist', 'Wollscheid', 'Wollscheid_HEK', 'Cut_struct', 'Molec_func']
	
base_cols.extend(log2FC_ls)
base_cols.extend(KL_cols)

all_name = base_name + '_All.xlsx'
ecd_name = base_name + '_ECD.xlsx'
sus_name = base_name + '_other_ECD.xlsx'
icd_name = base_name + '_ICD.xlsx'

TM_cols = ['Peptide', 'Accession', 'log2FC_med', 'Cut_topo', 'SS_dist', \
	'Peptide_start', 'Prev_AA', 'P4-P4\'', 'P1', 'TM_dist', 'nearest_TM', \
	'lipid_dist', 'nearest_lipid', 'nearest_lipid_type', \
	'glycosite_dist', 'nearest_glycosite', 'Protein_type', 'Location']
TM_name = base_name + '_TM_dist.xlsx'

domain_cols = ['Peptide','P4-P4\'', 'Accession', 'Domain_range', 'Cut_domain', \
	'Domains_removed', 'Num_domains_removed', \
	'Domains_intact', 'Num_domains_intact', \
	'Domain_dist', 'Intradomain_dist']
dom_name = base_name + '_dom_dist.xlsx'


ecd_writer = pd.ExcelWriter(ecd_name)
basic_ecd_df[base_cols].to_excel(ecd_writer, index = False, sheet_name = 'basic')
aliphatic_ecd_df[base_cols].to_excel(ecd_writer, index = False, sheet_name = 'aliphatic')
other_ecd_df[base_cols].to_excel(ecd_writer, index = False, sheet_name = 'other')
ecd_writer.save()

sus_writer = pd.ExcelWriter(sus_name)
basic_sus_df[base_cols].to_excel(sus_writer, index = False, sheet_name = 'basic')
aliphatic_sus_df[base_cols].to_excel(sus_writer, index = False, sheet_name = 'aliphatic')
other_sus_df[base_cols].to_excel(sus_writer, index = False, sheet_name = 'other')
sus_writer.save()


TM_writer = pd.ExcelWriter(TM_name)
basic_df[TM_cols].to_excel(TM_writer, index = False, sheet_name = 'basic')
aliphatic_df[TM_cols].to_excel(TM_writer, index = False, sheet_name = 'aliphatic')
other_df[TM_cols].to_excel(TM_writer, index = False, sheet_name = 'other')
TM_writer.save()

domain_writer = pd.ExcelWriter(dom_name)
basic_df[domain_cols].to_excel(domain_writer, index = False, sheet_name = 'basic')
aliphatic_df[domain_cols].to_excel(domain_writer, index = False, sheet_name = 'aliphatic')
other_df[domain_cols].to_excel(domain_writer, index = False, sheet_name = 'other')
domain_writer.save()
