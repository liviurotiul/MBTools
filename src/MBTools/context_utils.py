import pandas as pd
from tqdm import tqdm

def merge_overlapping_genes(df):
    # Sort DataFrame based on "Q_START" column

    def select_tag(tags):

        if "MGE" in tags:
            if "BR" in tags:
                return "BR"
            if "AR" in tags:
                return "AR"
            return "MGE"

        if "VF" in tags and "AR" in tags:
            return "VF"

        if "PAI" in tags and "VF" in tags:
            return "PAI"

        if "BR" in tags and "AR" in tags:
            return "BR"

        return tags[0]

    df = df.sort_values(by="Q_START")

    merged_predictions = []

    current_list = [df.iloc[0]['TAG']]
    current_end = df.iloc[0]["Q_END"]

    for index in range(1, len(df)):
        start_next = df.iloc[index]["Q_START"]
        end_next = df.iloc[index]["Q_END"]
        tag_next = df.iloc[index]['TAG']

        if current_end >= start_next:
            # Overlapping genes
            current_list.append(tag_next)
            current_end = max(current_end, end_next)
        else:
            # Non-overlapping genes
            merged_predictions.append(select_tag(current_list))
            current_list = [tag_next]
            current_end = end_next

    # Add the last merged prediction
    merged_predictions.append(select_tag(current_list))
    return merged_predictions


def concatenate_overlapping_genes(df):
    # Sort DataFrame based on "Q_START" column

    df = df.sort_values(by="Q_START")

    merged_predictions = []

    current_list = [df.iloc[0]['TAG']]
    current_end = df.iloc[0]["Q_END"]

    for index in range(1, len(df)):
        start_next = df.iloc[index]["Q_START"]
        end_next = df.iloc[index]["Q_END"]
        tag_next = df.iloc[index]['TAG']

        if current_end >= start_next:
            # Overlapping genes
            current_list.append(tag_next)
            current_end = max(current_end, end_next)
        else:
            # Non-overlapping genes
            merged_predictions.append('+'.join(current_list))
            current_list = [tag_next]
            current_end = end_next

    # Add the last merged prediction
    merged_predictions.append('+'.join(current_list))
    return merged_predictions


def concatenate_overlapping_gene_names(df):
    # Sort DataFrame based on "Q_START" column

    df = df.sort_values(by="Q_START")

    merged_predictions = []

    current_list = [df.iloc[0]['SHORT_ANNOTATION']]
    current_end = df.iloc[0]["Q_END"]

    for index in range(1, len(df)):
        start_next = df.iloc[index]["Q_START"]
        end_next = df.iloc[index]["Q_END"]
        tag_next = df.iloc[index]['SHORT_ANNOTATION']

        if current_end >= start_next:
            # Overlapping genes
            current_list.append(tag_next)
            current_end = max(current_end, end_next)
        else:
            # Non-overlapping genes
            merged_predictions.append('+'.join(current_list))
            current_list = [tag_next]
            current_end = end_next

    # Add the last merged prediction
    merged_predictions.append('+'.join(current_list))
    return merged_predictions


def create_zone_column(df, window):

    # Remember index
    old_index = df.index

    # The dataframe must be sorted by sample id, contig and start postion
    df = df.sort_values(by=["SAMPLE", "CONTIG", "Q_START"])

    zones = {}
    zone = []
    zone_id = 0
    zone_ids = []

    # Use tqdm to show progress
    # for index, annotation in enumerate(df.iterrows()):
    for index in tqdm(range(len(df))):
        zone_ids.append(zone_id)

        if index == len(df)-1:
            break
 
        if df.iloc[index]["CONTIG"] != df.iloc[index+1]["CONTIG"] or df.iloc[index]["SAMPLE"] != df.iloc[index+1]["SAMPLE"] or df.iloc[index]["Q_END"] + window < df.iloc[index+1]["Q_START"]:
            zone_id += 1
    
    df["ZONE_ID"] = zone_ids
    df = df.sort_index()

    return df["ZONE_ID"]


def drop_duplicate_annotations(df):

    # Make checks for column names
    # Make settings for id and coverage; criterions for droping

    df = df.drop_duplicates(subset=["Sample", "Contigs", "Q_START", "Q_END"])

    return df


def generate_tag(predictor):

    return predictor_dict.get(predictor.upper(), "UKN")


def create_tag_column(df):
    
    return df["PREDICTOR"].apply(generate_tag)


dbcolor_mapping = {
    "CGEVF": "purple",
    "VFDB": "purple",
    "BacAntRes": "red",
    "ResFinder": "red",
    "CGERF": "red",
    "CARD": "red",
    "BacMet": "red",
    "PAIDB": "black",
    "CGEMGE": "green",
    "MGEGoth": "green",
    "ISFinder": "green",
    "BacAntIn": "green",
    "BacAntTn": "green",
    "ICEberg2": "orange",
    "Phaster": "blue",
    "BacAntRep": "grey",
    "CGEPF": "grey",
    "CGESF": "grey"
}


predictor_dict = {
    "CARD": "AR",
    "VFDB": "VF",
    "BACMET": "BR",
    "PHASTER": "PHG",
    "ISFINDER": "IS",
    "ICEBERG2": "ICE",
    "PAIDB": "PAI",
    "CGEVF": "VF",
    "CGESF": "SF",
    "BACANTREP": "PLS",
    "CGESEROT": "SER",
    "CGMGE": "MGE",
    "MGEGOTH": "MGE",
    "CGEPF": "PLS",
    "BACANTTN": "TN",
    "BACANTRES": "AR",
    "CGERF": "AR",
    "CGEMGE": "MGE",
    "BACANTIN": "IS",
    "TADB_AT": "AT",
    "TADB_T": "T",
    "BACANTREP": "PLS"
}


def translator(diamond_annotation, predictor):

    # BacAntIn In511_AY551331_1
    # Ab@BacAntRep rep10_2_CDS1(pIM13)_M13761_1
    # Ab@BacAntRep IncL_1__JN626286_1
    # BacAntRes ncbi~~~abaF~~~CP000521.1~~~FOSFOMYCIN_1
    # BacMet BAC0016|adeH|tr|Q2FD80|Q2FD80_ACIBA
    # CARD HM|gb|AGV28567.1|ARO:3000559|adeN
    # CGEMGE Tn6172|1|KU744946|StrA
    # CGERF catA1_1_V00622_1
    # ICEberg2 ICEberg|1023
    # ISFinder html.2020//ISAba1
    # MGEGoth 2277_qacEdelta_CT025799.2_1
    # PAIDB gi_482862065_ref_AGK36638.1__OprE3_[Acinetobacter_baumannii]
    # Phaster PHAGE_Acinet_YMC11/11/R3177_NC_041866-gi_100018_ref_YP_009593347.1__hypothetical_protein_[Acinetobacter_phage_YMC11/11/R3177]
    # VFDB VFG037506(gb_WP_002004304)_(BJAB0715_RS05245)_transferrin-binding_protein-like_solute_binding_protein_[HemO_cluster_(VF0636)_-_Nutritional/Metabolic_factor_(VFC0272)]_[Acinetobacter_baumannii_BJAB0715]
    # CGEVF	acm_1_CP003351.1_1
    # CGESF O-Type_wzt_1_AB010150_O8_1
    # CGESeroT O-Type_wzt_1_AB010150_O8_1
    # CGEPF IncHI1B(pNDM-CIT)_1__JX182975_1
    # BacAntTn Tn602_AH000951_1
    # TADB_AT AT10068_WP_001288210.1_NZ_CP018664:c926440-926183_[Acinetobacter_baumannii]
    # TADB_T T3593_YP_001085049.1_NC_009085:c2346189-2346031_[Acinetobacter_baumannii_ATCC_17978]

    if predictor == "BacAntIn":
        return diamond_annotation.split("_")[0]

    if predictor == "BacAntRep":
        return diamond_annotation.split("_")[0]
    
    if predictor == "BacAntRes":
        return diamond_annotation.split("~~~")[1]
    
    if predictor == "BacMet":
        return diamond_annotation.split("|")[1]
    
    if predictor == "CARD":
        return diamond_annotation.split("|")[3]
    
    if predictor == "CGEMGE":
        return diamond_annotation.split("|")[0]
    
    if predictor == "CGERF":
        return diamond_annotation.split("_")[0]
    
    if predictor == "ICEberg2":
        return diamond_annotation

    if predictor == "ISFinder":
        return diamond_annotation.split("//")[1]
    
    if predictor == "MGEGoth":
        return diamond_annotation.split("_")[1]

    if predictor == "PAIDB":
        return diamond_annotation.split("__")[1].split("_")[0]

    if predictor == "Phaster":
        return "PHAGE" + diamond_annotation.split("_")[2]
    
    if predictor == "VFDB":
        return diamond_annotation.split("_")[3].replace("(", "").replace(")", "")
    
    if predictor in ["CGEVF", "CGESF", "CGESeroT", "CGEPF"]:
        return '_'.join(diamond_annotation.split("_")[0:1])
    
    if predictor == "BacAntTn":
        return diamond_annotation.split("_")[0]
    
    if predictor in ["TADB_AT", "TADB_T"]:
        return diamond_annotation.split("_")[0]