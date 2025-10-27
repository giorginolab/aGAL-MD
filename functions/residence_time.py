import pandas as pd
import numpy as np

structures = ['DGJ', 'DGJ_N215S', 'DGJ_R301Q']
replicas = ['1', '2', '3']

data_tot = []

for s in structures:
    for r in replicas:
        for resi in [1, 2]:
            df = pd.read_csv(f'../results/tables/{s}_{r}_lig_{resi}_rmsd.csv')

            # rmsd >= 5
            t_exit = df.loc[df['rmsd'] >= 5, 'time']
            r_time = t_exit.iloc[0] if not t_exit.empty else np.nan

            data = {
                'structure': s,
                'replica': r,
                'DGJ': resi,
                'time (ns)': r_time
            }
            data_tot.append(data)

df_results = pd.DataFrame(data_tot)
df_results.to_csv('../results/tables/residence_time.csv', index=False)