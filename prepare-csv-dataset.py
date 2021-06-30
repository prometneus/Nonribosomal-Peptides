import pandas as pd
import numpy as np
from Bio import SeqIO


DOWNLOAD_TAR = False
TAR_URL = "https://dl.secondarymetabolites.org/mibig/mibig_gbk_2.0.tar.gz"
TAR_FILENAME = TAR_URL.split('/')[-1]
TAR_FILEPATH = f"./data/{TAR_FILENAME}"

if DOWNLOAD_TAR:
    import requests
    import tarfile
    from tqdm import tqdm

    r = requests.get(TAR_URL, stream=True)

    with open(TAR_FILEPATH, "wb") as f:
        for data in tqdm(r.iter_content()):
            f.write(data)

    tar = tarfile.open(TAR_FILEPATH)
    tar.extractall("./data/gbk/")
    tar.close()


# Pandas display options
pd.set_option('display.max_colwidth', 30)
pd.set_option('display.width', 2000)
pd.set_option('display.max_columns', 10)

# Reading combined CSV with several sources
source_df = pd.read_csv("./data/A_domains.csv")
source_df = source_df.replace("$\s*-\s*^", np.nan, regex=True)

# Getting data from Prieto Sequence CSV (with sequences)
df_final_prieto = pd.read_csv("./data/Prieto_Adomain_Substrate.csv")


# Getting sequence data from MIBIG tar.gz GenBank files
mibig_df = source_df[source_df["Source"] == "MIBiG 2.0"]


def get_dna_from_row(row):
    gb_filename = row["Name"] + ".gbk"
    print(gb_filename)
    try:
        result_seq = ""
        for seq_record in SeqIO.parse(f"./data/gbk/{TAR_FILENAME.replace('.tar.gz', '')}/{gb_filename}", "genbank"):
            result_seq += seq_record.seq
            break
        start, end = row["Domain location"].split("-")
        dna = "".join(result_seq[int(start):int(end)])
        return dna
    except FileNotFoundError:
        return None


mibig_df['Sequence'] = mibig_df.apply(lambda row: get_dna_from_row(row), axis=1)
df_final_mibig = mibig_df[["Name", "Specificity", "Sequence"]]
df_final_mibig.columns = ["ID", "Substrate", "Sequence"]

print(df_final_prieto.columns)
print(df_final_mibig.columns)

final_df = pd.concat([df_final_prieto, df_final_mibig], axis=0, ignore_index=True)
final_df = final_df.dropna()

print(final_df.head())
final_df.to_csv("data/final.csv")
