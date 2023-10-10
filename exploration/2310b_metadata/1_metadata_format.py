# %%
import pandas as pd
import matplotlib.pyplot as plt

# %%

mtd_1 = "../../config/datasets/ST131_refseq/metadata.csv"
mtd_2 = "../../config/datasets/ST131_full/metadata.csv"
# Load CSV file into pandas dataframe
df1 = pd.read_csv(mtd_1, index_col=0)
df2 = pd.read_csv(mtd_2, index_col=0)
df = pd.concat([df1, df2], axis=0)

substitution_dict = {
    "year": {"missing": None},
    "host": {
        "missing": None,
        "pantropical spotted dolphin": "dolphin",
        "Homo sapiens; 55 year old": "human",
        "Pig": "pig",
        "Homo sapiens": "human",
    },
    "geo_loc_name": {
        "unknown": None,
        "missing": None,
        "not applicable": None,
    },
    "isolation_source": {
        "missing": None,
        "unknown": None,
        "clinical": "clinical",
        "blood": "blood",
        "urine": "urine",
        "human": "clinical",
        "wastewater": "wastewater",
        "environmental: river near WWTP effluent": "environment",
        "feces": "faeces",
        "environmental: river/lake": "environment",
        "environmental-WWTP": "wastewater",
        "human, urinary tract infection": "UTI",
        "water": "environment",
        "Urinary tract": "UTI",
        "wastewater treatment plant efflux": "wastewater",
        "sputum": "sputum",
        "cecal sample": "faeces",
        "wastewater treatment plant effluent": "wastewater",
        "Hospital Universitario Central de Asturias (Oviedo)": "clinical",
        "Hospital Universitario Gregorio Maranon (Madrid)": "clinical",
        "stool": "faeces",
        "rectal swab": "clinical",
        "vaginal swab": "UTI",
        "wastewater efflux": "wastewater",
        "UTI": "UTI",
        "dog, faeces": "faeces",
        "Hera General Hopsital, patient": "clinical",
        "abscess": "clinical",
        "Endometrium": "clinical",
        "Rectal swab": "clinical",
        "ascites fluid": "clinical",
        "drain drainage": "clinical",
        "human clinical specimen": "clinical",
        "urine (asymptomatic UTI)": "UTI",
        "button": None,
        "ICU patient rectal swab": "clinical",
        "ICU staff rectal swab": "clinical",
        "bile duct discharge": "clinical",
        "Chronic wound": "clinical",
        "hospital": "clinical",
        "urine (recurrent UTI)": "UTI",
        "Swab from Wound": "clinical",
        "Rectal": "clinical",
        "not collected": None,
        "body fluids": "clinical",
    },
}

continents = {
    "USA": "North America",
    "Switzerland": "Europe",
    "Sweden": "Europe",
    "China": "Asia",
    "Japan": "Asia",
    "United Kingdom": "Europe",
    "Czech Republic": "Europe",
    "Australia": "Australia",
    "Spain": "Europe",
    "Argentina": "South America",
    "Armenia": "Asia",
    "Canada": "North America",
    "Thailand": "Asia",
    "Germany": "Europe",
    "Ghana": "Africa",
    "Italy": "Europe",
    "New Zealand": "Australia",
    "Mexico": "North America",
    "Kenya": "Africa",
    "Cambodia": "Asia",
    "South Korea": "Asia",
    "Kazakhstan": "Asia",
    "Hong Kong": "Asia",
    "United Arab Emirates": "Asia",
    "Saudi Arabia": "Asia",
    "Lao": "Asia",
    "Brazil": "South America",
    "South Africa": "Africa",
    "Colombia": "South America",
    "Viet Nam": "Asia",
    "Singapore": "Asia",
    "Morocco": "Africa",
    "Netherlands": "Europe",
}

# only consider country for geo-location
df["geo_loc_name"] = df["geo_loc_name"].str.split(":", expand=True)[0]

for k, d in substitution_dict.items():
    df[k] = df[k].replace(d)

df["continent"] = df["geo_loc_name"].replace(continents)
# %%
