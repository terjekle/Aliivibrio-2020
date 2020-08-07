#!/usr/local/bin/python3
from Bio import SeqIO
from io import BytesIO
from difflib import SequenceMatcher
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import time

start_time = time.time()
fileInput = "input.fasta"
ymaxvalue = 25 # KDE plot Y max limit value

### Function for removal of pairwise gaps and calculation of sequence variance.
def rawSeqInp(sequence_a, sequence_b):
# A while loop continously search for an existing gap in sequence_a and if found
# remove its location from both sequences.
	gaps = sequence_a.count('-')
	count = 0
	while count <= int(gaps):
		count += 1
		result = sequence_a.find('-')
		if result == -1:
			continue
		else:
			sequence_a = sequence_a[:result] + sequence_a[result+1:]
			sequence_b = sequence_b[:result] + sequence_b[result+1:]
	gaps = sequence_b.count('-')
	count = 0

# A while loop continously search for an existing gap in sequence_b and if found
# remove its location from both sequences.
	while count <= int(gaps):
		count += 1
		result = sequence_b.find('-')
		if result == -1:
			continue
		else:
			sequence_a = sequence_a[:result] + sequence_a[result+1:]
			sequence_b = sequence_b[:result] + sequence_b[result+1:]

# SequenceMatcher calculate the ratio of difference between sequence a and b.
	differ = SequenceMatcher(None, sequence_a.upper(), sequence_b.upper(),
                             autojunk=False).ratio()
	sequence_set = {'seq_1': sequence_a, 'seq_2': sequence_b,
                    'difference': differ}
	return sequence_set

### End of function \\\



### Import & management
print("\n\n\tInfo:\n\tIUPAC characters other than ATGC will not be " + \
      "accounted for.")

alignment = {} # Headers and sequences stored as a dictionary.
headers   = [] # List of headers
sequences = [] # List of sequences
pairwiseD = [] # List of paired distances for printing (Header1 dist Header2)
pooled_identities = [] # List of all identity values between sequence pairs.

for seq_record in SeqIO.parse(fileInput, "fasta"):
	headers.append(seq_record.description)
	sequences.append(seq_record.seq)
	alignment[seq_record.description] = seq_record.seq

print('\n\tNumber of sequences to process:\t' + str(len(headers)) + \
      '\n\t! Processing many long sequences can take a while.\n')

# Creating an empty data frame to append calculated identity values.

d = {}
df = pd.DataFrame(d)

# Adding headers in first column of data frame.

df['strains']=headers

# A while loop runs a diminishing process on all sequences, popping out headers.
# The header corresponds to a sequence stored in the dictionary, while remaining
# headers correspond to sequences withdrawn for comparison.
while len(headers) > 0:
	header = headers.pop(0)
	print('\tProcessing:\t' + header)
	elapsed_time = time.time() - start_time
	print("\tTime passed:\t" + str(time.strftime("%H:%M:%S", time.gmtime(elapsed_time))))
	identities = []
	sequence_a = alignment[header]

	counter = 0
	for ramaining_header in headers:
		sequence_b = alignment[ramaining_header]
		sequence_set = rawSeqInp(sequence_a, sequence_b)
		identities.append(sequence_set['difference'])
		pooled_identities.append(sequence_set['difference'])
		pairwiseD.append(header + '\t' + \
                         str(format(sequence_set['difference'] * 100, '7.2f')) + \
                         '\t' + ramaining_header)

	while len(identities) < len(sequences):
		identities.insert(0, 'NaN')
	df[header]=identities

### End of import & management \\\



### Output writer
print("\t" + str(round((min(pooled_identities)),5) * 100))

# Writing all paired distances to file, 'Header1 distance Header2'.
f = open("out_pooled_identities.txt","w+")

# Summing up avg identity, min and max.
f.write("\nAverage identity\tLowest identity\tHighest identity\n")
f.write(str(round((sum(pooled_identities)/len(pooled_identities)),5) * 100) + "\t" + \
        str(round((min(pooled_identities)),5) * 100) + "\t" + str(round(max(pooled_identities),5) * 100) + \
		"\n\n")
for pairs in pairwiseD:
	f.write(pairs + '\n')

print("\n\tPooled identities exported to 'out_pooled_identities.txt'")

# Writinng a matrix of all pairwise distances to file.
export_csv = df.to_csv (r'out_identities.csv', index = None, header=True)

print("\tData frame exported to 'out_identities.csv'")

# Generating kernel density estimate (KDE) distribution plot of stored values.
plt.style.use('seaborn-darkgrid')
plt.xlabel('Identity')
plt.ylabel('Density')
plt.title('Kernel density estimate')
plt.figure(1)

sns.kdeplot(pooled_identities, bw=.2, shade=True, cut=0)

axes = plt.gca()
axes.vlines(x=0.97, ymin=0, ymax=50, colors='red', linestyles='dashed')
axes.set_xlim([0,1.00])
axes.set_ylim([0,ymaxvalue])

# Exporting plot as 300 DPI figure
f = BytesIO()
plt.savefig('out_identity_distribution_KDA.svg', format="svg", dpi=300)

print("\t\nKernel density estimate (KDE) plot exported to" + \
      "'out_identity_distribution_KDA.svg'\n\n")

### End of output writer \\\
