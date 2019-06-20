import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# make a dotplot x-slen, y-qcovs
f = 'data-PflB.xlsx'
df_bbh = pd.read_excel(f, sheet_name='BBH')
df_bbh = df_bbh.loc[:, ['slen', 'qcovs']]
x_bbh = df_bbh['slen'].apply(lambda x: x/820)
y_bbh = df_bbh['qcovs']

df_psi = pd.read_excel(f, sheet_name='psi')
df_psi = df_psi.loc[:, ['slen', 'qcovs']]
x_psi = df_psi['slen'].apply(lambda x: x/820)
y_psi = df_psi['qcovs']

#~~~~~~~~~~~~~~~~~~~~~~~
# MAKE THE SCATTER PLOT:
#~~~~~~~~~~~~~~~~~~~~~~~

#make the plot
fig, ax = plt.subplots()

ax.scatter(x_bbh, y_bbh, color='#005bd3', marker='o', edgecolor='k', linewidth=0.7, label='bidirectional blastp hits')
ax.scatter(x_psi, y_psi, color='C2', marker='^', edgecolor='k', linewidth=0.7, label='psiblast hits')


# axis labels
plt.xlabel('hit length/ query length ratio')
plt.ylabel('% query coverage')

#set xticks
plt.xticks(np.arange(0.5, 1.1, 0.1))
ax.set_xlim(0.47, 1.13)
plt.yticks(range(60, 110, 10))
ax.set_ylim(58)

ax.xaxis.set_minor_locator(AutoMinorLocator(n=2))
ax.yaxis.set_minor_locator(AutoMinorLocator(n=2)) # turn on minor ticks to be for every 1 bin

#turn on y-axis gridlines
ax.yaxis.grid(b=True, which='major', linestyle='--', color='0.9')
ax.yaxis.grid(b=True, which='minor', linestyle='--', color='0.9')
ax.set_axisbelow(True)


#figure size
fig.set_size_inches(6, 4)

#position legend, show plot
#plt.margins(1)
plt.legend()
plt.tight_layout()
plt.savefig('PflB_slen_qcovs_scatter', dpi=300)
plt.show()
