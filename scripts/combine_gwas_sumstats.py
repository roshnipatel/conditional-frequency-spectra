import pandas as pd

gwas_dir = "/oak/stanford/groups/pritch/users/roshni/stabilizing_selection/data/gwas/"

# load quantitative trait sum stats from Yuval (downloaded via Slack)
quant_gwas = pd.read_csv(gwas_dir + "raw/hits_for_roshni.csv")[["SNP", "refA", "b", "trait"]]
quant_gwas["trait"] = quant_gwas.trait - 1
quant_trait_table = pd.read_csv(gwas_dir + "raw/trait_names_for_roshni.csv", names=["full_trait"])
quant_gwas = pd.merge(quant_gwas, quant_trait_table, left_on="trait", right_index=True)

# load binary trait sum stats from Julie (from Oak directory)
binary_traits = {"adhd": "Attention deficit hyperactivity disorder", 
                 "ad": "Alzheimer's disease", 
                 "bip": "Bipolar disorder",
                 "brca": "Breast cancer",
                 "cad": "Coronary artery disease",
                 "ibd": "Inflammatory bowel disease",
                 "ms": "Multiple sclerosis",
                 "pd": "Parkinson's disease",
                 "prca": "Prostate cancer",
                 "ra": "Rheumatoid arthritis",
                 "scz": "Schizophrenia",
                 "t2d": "Type 2 diabetes"}
binary_gwas_list = []
for trait_abbv in binary_traits:
    curr_gwas = pd.read_csv(gwas_dir + "raw/" + trait_abbv + ".txt", sep=' ')[["SNP", "effect_allele", "b_cojo"]]
    curr_gwas["full_trait"] = binary_traits[trait_abbv]
    binary_gwas_list.append(curr_gwas)
binary_gwas = pd.concat(binary_gwas_list)

# need to match effect sizes to alternate allele
quant_gwas["alt"] = quant_gwas.apply(lambda x: x.SNP.split(':')[-1], axis=1)
quant_gwas["ref"] = quant_gwas.apply(lambda x: x.SNP.split(':')[-2], axis=1)
quant_gwas_effect_alt = quant_gwas[quant_gwas.refA.str.upper() == quant_gwas.alt]
quant_gwas_effect_alt["effect"] = quant_gwas_effect_alt.b
quant_gwas_effect_ref = quant_gwas[quant_gwas.refA.str.upper() == quant_gwas.ref]
quant_gwas_effect_ref["effect"] = -quant_gwas_effect_ref.b
quant_gwas = pd.concat([quant_gwas_effect_alt[["SNP", "effect", "full_trait"]],
                        quant_gwas_effect_ref[["SNP", "effect", "full_trait"]]])
quant_gwas["binary"] = 0

binary_gwas["alt"] = binary_gwas.apply(lambda x: x.SNP.split(':')[-1], axis=1)
binary_gwas["ref"] = binary_gwas.apply(lambda x: x.SNP.split(':')[-2], axis=1)
binary_gwas_effect_alt = binary_gwas[binary_gwas.effect_allele.str.upper() == binary_gwas.alt]
binary_gwas_effect_alt["effect"] = binary_gwas_effect_alt.b_cojo
binary_gwas_effect_ref = binary_gwas[binary_gwas.effect_allele.str.upper() == binary_gwas.ref]
binary_gwas_effect_ref["effect"] = -binary_gwas_effect_ref.b_cojo
binary_gwas = pd.concat([binary_gwas_effect_alt[["SNP", "effect", "full_trait"]],
                         binary_gwas_effect_ref[["SNP", "effect", "full_trait"]]])
binary_gwas["binary"] = 1

all_gwas = pd.concat([quant_gwas, binary_gwas])
all_gwas["trait_idx"] = all_gwas.groupby("full_trait").ngroup()

all_gwas.to_csv(gwas_dir + "combined_gwas.txt", sep='\t', index=False)
all_gwas[["trait_idx", "full_trait", "binary"]].drop_duplicates().to_csv(gwas_dir + "combined_trait_table.txt", sep='\t', index=False)