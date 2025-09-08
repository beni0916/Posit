import pandas as pd
from sklearn.model_selection import train_test_split

url = "https://raw.githubusercontent.com/selva86/datasets/master/AmesHousing.csv"
df = pd.read_csv(url)

X = df.drop(columns=["SalePrice"])
y = df["SalePrice"]

df_all = pd.concat([X, y.rename("SalePrice")], axis=1)

df_train, df_test = train_test_split(df_all, test_size=0.2, random_state=42)


df_train.to_csv("AHD_train.txt", sep='\t', index=False)
df_test.to_csv("AHD_test.txt", sep='\t', index=False)

print(" AHD_train.txt  AHD_test.txt")
