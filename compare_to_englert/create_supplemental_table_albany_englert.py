import pandas as pd
import sys

gene_sets_up_albany_f = sys.argv[1]
gene_sets_down_albany_f = sys.argv[2]

df_up = pd.read_csv(gene_sets_up_albany_f, sep='\t')
df_down = pd.read_csv(gene_sets_down_albany_f, sep='\t')

df_up.insert(len(df_up.columns), 'de_genes', ['up_in_covid_ards' for i in range(len(df_up))])
df_down.insert(len(df_down.columns), 'de_genes', ['down_in_covid_ards' for i in range(len(df_down))])

df_up = df_up.set_index(['gene_set', 'de_genes'])
df_down = df_down.set_index(['gene_set', 'de_genes'])

df_final = pd.concat([df_up, df_down])

df_final.to_csv('Supplemental_Table_5_albany_vs_englert.tsv', sep='\t')


