import pandas as pd

# Pfad zur Tabelle
table_path = './data/RBP2GO/Table_HS_RBP.txt'

# Tabelle mit pandas parsen
df = pd.read_csv(table_path, delimiter='\t')

# Ausgabe der Tabelle
print(df.head(1))