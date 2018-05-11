#######################################################################
# Postprocessing for Spel et al. (2018)                               #
# Elmer Stickel, Netherlands Cancer Institute, 2018.                  #
# Contact: elmer.stickel [at] posteo.net                              #
#######################################################################

# Import libraries
import pandas as pd
import numpy as np
import random

# Load results from R into Pandas Dataframe
raw = pd.read_csv("intermediate_results.csv", sep=',', names=['Gene_sgRNA', 'baseMean', 'log2FoldChange', 'stat', 'pvalue', 'padj'], skiprows=1, quotechar='"')

# Split dataframe in targeting and non-targeting guides
results = raw[~raw.Gene_sgRNA.str.match('NonTargetingControlGuideForHuman')]
nontargeting = raw[raw.Gene_sgRNA.str.match('NonTargetingControlGuideForHuman')]

# Randomly assign the NonTargetingGuides over 167 hypothetical NonTargetingControlGenes (=6 Guides per Gene)
r = [i for i in range(1,168)]*6
random.shuffle(r)
nontargeting.reset_index(inplace=True, drop=True)
nontargeting = nontargeting.join(pd.DataFrame(r[:1000]))
nontargeting = nontargeting.assign(Gene=('NonTargetingControlGuideForHuman-' + nontargeting[0].astype(str)))
nontargeting = nontargeting.assign(sgRNA=(nontargeting.Gene_sgRNA.str.split("_", expand=False).str.get(2)))
nontargeting = nontargeting.drop([0, 'Gene_sgRNA'], 1)

# Split Gene_sgRNA string into separate fields
results = results.assign(Gene=(results.Gene_sgRNA.str.split("_", expand=False).str.get(0)))
results = results.assign(sgRNA=(results.Gene_sgRNA.str.split("_", expand=False).str.get(1)))
results = results.drop('Gene_sgRNA', 1)

# Merge
results = pd.concat([results, nontargeting])

hits = []
counts = results[results.padj<=0.05]['Gene'].value_counts() # Number of sgRNA per gene matching the p-value cutoff
for i in list(counts[counts>2].keys()): # Loop genes having more than 2 (ie. at least 3) significant guides and check of the Log2FoldChange of their guides have the same sign
    log2FoldChanges = results[(results.padj<=0.05) & (results.Gene == i)].log2FoldChange.values
    hits.append(i) if all(item >= 0 for item in log2FoldChanges) or all(item < 0 for item in log2FoldChanges) else None


results['hit'] = np.where(np.in1d(results.Gene, hits), True, False)
results.to_csv('results.csv', sep='\t', header=True, index=False, quoting=2)