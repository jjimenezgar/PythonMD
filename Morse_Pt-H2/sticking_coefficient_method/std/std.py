import pandas as pd
import numpy as np


df1 = pd.read_csv("/home/johnconnor/Documentos/PtH_MD/pigeon_morse/std/filter/sco_final.txt", sep= ",")

df1.columns= ["Ek","De1", "De2", "Alfa1","Alfa2" ,"Re2","SC"]

# Agrupar por los valores de las primeras cinco columnas y calcular la media y la desviación estándar de "SC"
df_grouped = df1.groupby(["Ek", "De1", "De2", "Alfa1", "Alfa2","Re2"])["SC"].agg(['mean', 'std']).reset_index()

# Renombrar las columnas de la media y la desviación estándar
df_grouped.rename(columns={'mean': 'SC_mean', 'std': 'SC_std'}, inplace=True)


df_grouped.loc[df_grouped["Ek"] == 33, "error_exp"] = abs (df_grouped.loc[df_grouped["Ek"] == 33, "SC_mean"] - 0.05)
df_grouped.loc[df_grouped["Ek"] == 44, "error_exp"] = abs (df_grouped.loc[df_grouped["Ek"] == 44, "SC_mean"] - 0.30)
df_grouped.loc[df_grouped["Ek"] == 54, "error_exp"] = abs (df_grouped.loc[df_grouped["Ek"] == 54, "SC_mean"] - 0.50)
#df_grouped.loc[df_grouped["Ek"] == 63, "error_exp"] = abs (df_grouped.loc[df_grouped["Ek"] == 63, "SC_mean"] - 0.60)

grupos = df_grouped.groupby(['De1','De2', 'Alfa1', 'Alfa2', 'Re2'])


for i, (grupo_name, grupo) in enumerate(grupos):
    
    mean_exp = (np.mean(grupo["error_exp"]))
    mean_std = (np.mean(grupo["SC_std"]))

    if mean_std < 0.06:
        grupo.to_csv("/home/johnconnor/Documentos/PtH_MD/pigeon_morse/std/filter/{}.csv".format(i), index=False,header= None, sep=" ")



#print(df_grouped)