#!/usr/bin/env python3

import pandas as pd

excel = "data/supplementary_materials.xlsx"
df_dict = pd.read_excel(excel, sheet_name=None)

# saving each sheet as a separate TSV file
for sheet in df_dict:
    df_dict[sheet].to_csv(f"data/tables/{sheet}.tsv", index=False)