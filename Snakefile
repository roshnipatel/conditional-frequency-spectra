import pandas as pd

shell.executable("/usr/bin/bash")
shell.prefix("source /home/users/rpatel7/.bashrc; ")

CHROMS=range(1, 23)
POPS=["CHB", "YRI"]

# Input files needed:
### "data/CADD_bestfit" - B-statistics downloaded from Sella lab GitHub
### "data/combined_gwas.txt" - GWAS summary statistics including SNP, effect, and trait ID
### "data/snp_ID_conversion_table.txt" - file linking SNP and rsID
### "data/1KG_sample_info.txt" - file linking 1K Genomes sample IDs to populations
### "data/1KGenomes" - 1K Genomes hg19/GRCh37 VCFs
### "data/variation_feature.txt.gz" - Ensembl variation_feature table specifying ancestral allele states
### "data/freq_WB" - SNP frequencies in UK Biobank White British 

rule all:
    input:
        expand("data/distributions/empirical/{pref}_gwas_pvalues.txt", 
               pref=["all", "binary0_joint", "noCpG", "binary1_negative", "binary1_positive"] + 
                    ["trait" + str(i) + "_negative" for i in [90, 96, 63, 77, 7, 16, 51, 35, 27, 86, 93, 30, 23, 42]] +
                    ["trait" + str(i) + "_positive" for i in [90, 96, 63, 77, 7, 16, 51, 35, 27, 86, 93, 30, 23, 42]]),
        expand("data/simulations/generation2k_ancestral1e4/{sim}/{s}/ancestral_freq{freq}.txt",
               sim=["modern1e4_growth0",
                    "modern1e4_growth0.001",
                    "modern1e3_growth0",
                    "modern1e3_growth0.001"],
               s=["h0.5_s0.0",
                  "h0.5_s-1.0e-3",
                  "h0.5_s+1.0e-3",
                  "h5e6_s-1.0e-10"],
               freq=["%.2f" % (i/100) for i in range(1, 100)] + ["0.0001", "0.0005", "0.001", "0.002", "0.005", "0.015", "0.025"]),
        expand("data/simulations/jouganous_wo_migration/{s}/ancestral_freq{freq}.txt",
               s=["h0.5_s0.0",
                  "h0.5_s+1.0e-3", "h0.5_s+1.0e-4", "h0.5_s+5.0e-4",
                  "h0.5_s-1.0e-3", "h0.5_s-1.0e-4", "h0.5_s-5.0e-4",
                  "h5e6_s-5.0e-11", "h5e6_s-1.0e-11", "h5e6_s-1.0e-10"],
               freq=["%.2f" % (i/100) for i in range(1, 100)] + ["0.00003", "0.0001", "0.0005", "0.001", "0.002", "0.005", "0.015", "0.025"]),
        expand("data/distributions/{sim}/{mode}_ancestor_conditional_p2.txt",
               sim=["generation2k_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0",
                    "generation2k_ancestral1e4.p1_ne1e4_growth0.p2_ne1e3_growth0.001",
                    "generation2k_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0.001",
                    "generation2k_ancestral1e4.p1_ne1e4_growth0.p2_ne1e4_growth0"],
               mode=["h0.5_s0.0",
                     "h0.5_s-1.0e-3",
                     "h0.5_s+1.0e-3",
                     "h5e6_s-1.0e-10"]),
        expand("data/distributions/{sim}/{mode}_ancestor_conditional_ceu.txt",
               sim=["dtwf-jouganous_wo_migration"],
               mode=["h0.5_s0.0",
                     "h0.5_s+1.0e-3", "h0.5_s+5.0e-4", "h0.5_s+1.0e-4",
                     "h0.5_s-1.0e-3", "h0.5_s-5.0e-4", "h0.5_s-1.0e-4",
                     "h5e6_s-1.0e-10", "h5e6_s-5.0e-11", "h5e6_s-1.0e-11"]),
        expand("data/distributions/{sim}/{mode}_ancestor_conditional_ceu.txt",
               sim=["jouganous_wo_migration"],
               mode=["h0.5_s0.0",
                     "h0.5_s+1.0e-3", "h0.5_s+5.0e-4", "h0.5_s+1.0e-4",
                     "h0.5_s-1.0e-3", "h0.5_s-5.0e-4", "h0.5_s-1.0e-4",
                     "h5e6_s-1.0e-10", "h5e6_s-5.0e-11", "h5e6_s-1.0e-11"]),
        expand("data/distributions/dtwf-generation{gen}_ancestral{ancestral}.p1_ne{p1}_growth0.p2_ne{p2}/{mode}_ancestor_conditional_p2.txt",
               mode=["h0.5_s0.0",
                     "h0.5_s+1.0e-3",
                     "h0.5_s-1.0e-3",
                     "h5e6_s-1.0e-10"],
               gen=["2e3"],
               ancestral=["1e4"],
               p1=["1e4"],
               p2=["1e3_growth0", "1e4_growth0", "1e4_growth0.0005", "1e4_growth0.001", "1e3_growth0.001", "3e3_growth0"])

rule process_ukb_data: 
    input:
        "data/freq_WB/chr{chr}.afreq"
    output:
        "data/chr{chr}/all_ukb_snps.txt"
    shell:
        """
        cut -f2,5 {input} > {output}
        """

rule sort_rsID_table:
    input:
        "data/snp_ID_conversion_table.txt"
    output:
        "data/snp_ID_conversion_table_sorted.txt"
    shell:
        """
        sort -k1,1 {input} > {output}
        """

rule merge_snps_rsID:
    input:
        snps="data/chr{chr}/all_ukb_snps.txt",
        table=rules.sort_rsID_table.output
    output:
        "data/chr{chr}/all_ukb_snps_rsID.txt"
    shell:
        """
        sort -k1,1 {input.snps} | join -j 1 -t'\t' {input.table} - > {output}
        """

rule resort_snps:
    input:
        rules.merge_snps_rsID.output
    output:
        "data/chr{chr}/all_ukb_snps_rsID_sorted.txt"
    run:
        snps = pd.read_csv(input[0], sep='\t', names=["SNP", "rsID", "alt_freq"])
        snps = snps.drop_duplicates()
        snps["position"] = snps.apply(lambda row: row.SNP.split(':')[1], axis=1)
        snps = snps.sort_values(by="position")
        snps.to_csv(output[0], sep='\t', index=False, columns=["position", "SNP", "rsID", "alt_freq"])   

rule map_B_vals:
    input:
        B="data/CADD_bestfit/chr{chr}.bmap.txt",
        snps=rules.resort_snps.output
    output:
        "data/chr{chr}/all_ukb_snps_rsID_sorted_B.txt"
    shell:
        """
        conda activate py3
        python scripts/map_b_vals.py --snps {input.snps} \
            --bvals {input.B} \
            --out {output}
        conda deactivate
        """

rule process_variation_file:
    input:
        "data/variation_feature.txt.gz"
    output:
        "data/variation_feature.txt"
    shell:
        """
        zcat {input} | cut -f8,9 | grep -v 'N' > {output}
        """

rule aggregate_UKB:
    input:
        expand("data/chr{chr}/all_ukb_snps_rsID_sorted_B.txt", chr=CHROMS)
    output:
        "data/all_ukb_snps_rsID_sorted_B.txt"
    shell:
        """
        head -n 1 {input[0]} > {output}
        tail -n +2 -q {input} >> {output}
        """

rule subset_variation_file:
    input:
        snps=rules.aggregate_UKB.output,
        table="data/variation_feature.txt"
    output:
        "data/variation_feature_subsetted.txt"
    shell:
        """
        conda activate pandas-env
        python scripts/subset_big_table.py --snps {input.snps} \
            --table {input.table} --out {output}
        conda deactivate
        """

rule sort_variation_file:
    input:
        "data/variation_feature_subsetted.txt"
    output:
        "data/variation_feature_sorted.txt"
    shell:
        """
        sort -k2,2 {input} > {output}
        """

rule fetch_ancestral_states:
    input:
        snps=rules.aggregate_UKB.output,
        table=rules.sort_variation_file.output
    output:
        "data/all_ukb_snps_ancestral_allele.txt"
    shell: 
        """
        cut -f2-5 {input.snps} | sort -k2,2 | join -j 2 -t'\t' {input.table} - > {output}
        """ 

rule generate_control_snps:
    input:
        gwas="data/combined_gwas.txt",
        snps=rules.fetch_ancestral_states.output
    output:
        gwas="data/all_gwas.txt",
        control="data/all_control.txt"
    shell:
        """
        conda activate pandas-env
        python scripts/generate_controls.py --snps {input.snps} \
            --gwas {input.gwas} \
            --out_gwas {output.gwas} \
            --out_control {output.control}
        conda deactivate
        """

rule create_ref_filter:
    input:
        # Using the sample info from hg38 download to map individuals to pops; actual VCFs used are hg19/GRCh37
        "data/1KG_sample_info.txt"
    output:
        "data/filter_{grp}.txt"
    shell:
        """
        grep {wildcards.grp} {input} | cut -f1 > {output}
        """

rule split_chrom:
    input:
        "data/{snp_type}.txt"
    output:
        "data/chr{chr}/chr{chr}_{snp_type}.txt"
    shell:
        """
        grep ^{wildcards.chr}: {input} > {output}
        """

rule create_pos_file:
    input:
        "data/chr{chr}/chr{chr}_{snp_type}.txt"
    output:
        "data/chr{chr}/{snp_type}_positions.txt"
    run:
        df = pd.read_csv(input[0], sep='\t', names=["SNP", "rsID", "alt_freq", "ancestral", "B"])
        df[["chrom", "pos", "ref", "alt"]] = df.apply(lambda row: row.SNP.split(':'), axis=1, result_type="expand")
        df = df.sort_values(by="pos")
        df.to_csv(output[0], columns=["chrom", "pos"], header=False, index=False, sep='\t')

rule extract_positions:
    input:
        pos="data/chr{chr}/{snp_type}_positions.txt",
        vcf="data/1KGenomes/chr{chr}/chr{chr}.filt.1kg.phase3.v5a.biSNPs.vcf.gz"
    output:
        "data/chr{chr}/chr{chr}.extracted_{snp_type}.bcf.gz"
    shell:
        """
        conda activate bcftools-env
        bcftools view -R {input.pos} -T {input.pos} -Ob -o {output} {input.vcf}
        conda deactivate
        """

rule filter_ref:
    input:
        ref=rules.extract_positions.output,
        filter="data/filter_{grp}.txt"
    output:
        bcf="data/chr{chr}/chr{chr}.pop{grp}.extracted_{snp_type}.bcf.gz"
        # idx="data/chr{chr}/chr{chr}.pop{grp}.extracted_{snp_type}.bcf.gz.csi"
    shell:
        """
        conda activate bcftools-env
        bcftools view -S {input.filter} --force-samples -Ou {input.ref} | \
            bcftools view --genotype ^miss --phased -Ob -o {output.bcf}
        conda deactivate
        """

rule calculate_frequency:
    input:
        rules.filter_ref.output.bcf
    output:
        "data/chr{chr}/chr{chr}.pop{grp}.extracted_{snp_type}.frq"
    params:
        prefix="data/chr{chr}/chr{chr}.pop{grp}.extracted_{snp_type}"
    shell:
        """
        conda activate vcftools-env
        vcftools --bcf {input} --freq --out {params.prefix}
        conda deactivate
        """

rule compute_DAF:
    input:
        freq=expand("data/chr{{chr}}/chr{{chr}}.pop{grp}.extracted_{{snp_type}}.frq", grp=POPS),
        anc=rules.split_chrom.output
    output:
        "data/chr{chr}/{snp_type}_DAF.txt"
    shell:
        """
        conda activate pandas-env
        python scripts/merge_frequency_tables.py --ancestral_file {input.anc} \
            --frequency_files {input.freq} \
            --out {output}
        conda deactivate
        """

rule aggregate_snps:
    input:
        expand("data/chr{chr}/{{snp_type}}_DAF.txt", chr=CHROMS)
    output:
        "data/pooled_{snp_type}_DAF.txt"
    shell:
        """
        head -n 1 {input[0]} > {output}
        tail -n +2 -q {input} >> {output}
        """

rule drop_CpG:
    input:
        "data/pooled_all_{snp_type}_DAF.txt"
    output:
        "data/pooled_noCpG_{snp_type}_DAF.txt"
    shell:
        """
        head -n 1 {input} > {output}
        awk -F'\t' '(($8 == "C" || $8 == "T") && ($9 == "G" || $9 == "A")) || \
        (($8 == "G" || $8 == "A") && ($9 == "T" || $9 == "C"))' {input} >> {output}
        """

rule filter_gwas_by_effect_sign:
    input:
        gwas="data/combined_gwas.txt",
        snps="data/pooled_all_gwas_DAF.txt"
    output:
        "data/pooled_{sign_pref,[a-z]+}tive_gwas_DAF.txt"
    run:
        gwas = pd.read_csv(input.gwas, sep='\t')[["SNP", "effect"]] # pull out SNP and effect of alternate allele
        snps = pd.read_csv(input.snps, sep='\t')
        merged = pd.merge(gwas, snps, on="SNP")

        if wildcards.sign_pref == "posi":
            merged = merged[((merged.effect > 0) & (merged.alt == merged.derived)) | ((merged.effect < 0) & (merged.alt == merged.ancestral))] 
        elif wildcards.sign_pref == "nega":
            merged = merged[((merged.effect < 0) & (merged.alt == merged.derived)) | ((merged.effect > 0) & (merged.alt == merged.ancestral))] 
    
        merged = merged.drop(columns=["effect"]).drop_duplicates()
        merged.to_csv(output[0], sep='\t', index=False)

rule filter_gwas_by_trait:
    input:
        gwas="data/combined_gwas.txt",
        snps="data/pooled_all_gwas_DAF.txt"
    output:
        negative="data/pooled_trait{trait}_negative_gwas_DAF.txt",
        positive="data/pooled_trait{trait}_positive_gwas_DAF.txt"
    run:
        gwas = pd.read_csv(input.gwas, sep='\t')[["SNP", "trait_idx", "effect"]] # pull out SNP and effect of alternate allele
        gwas = gwas[gwas.trait_idx == int(wildcards.trait)]
        snps = pd.read_csv(input.snps, sep='\t')
        merged = pd.merge(gwas, snps, on="SNP")

        negative = merged[((merged.effect < 0) & (merged.alt == merged.derived)) | ((merged.effect > 0) & (merged.alt == merged.ancestral))] 
        positive = merged[((merged.effect > 0) & (merged.alt == merged.derived)) | ((merged.effect < 0) & (merged.alt == merged.ancestral))] 

        negative = negative.drop(columns=["trait_idx", "effect"])
        negative.to_csv(output.negative, sep='\t', index=False)
        positive = positive.drop(columns=["trait_idx", "effect"])
        positive.to_csv(output.positive, sep='\t', index=False)

rule filter_gwas_by_trait_type:
    input:
        gwas="data/combined_gwas.txt",
        snps="data/pooled_all_gwas_DAF.txt"
    output:
        negative="data/pooled_binary{idx}_negative_gwas_DAF.txt",
        positive="data/pooled_binary{idx}_positive_gwas_DAF.txt",
        all="data/pooled_binary{idx}_joint_gwas_DAF.txt"
    run:
        gwas = pd.read_csv(input.gwas, sep='\t')[["SNP", "binary", "effect"]] # pull out SNP and effect of alternate allele
        gwas = gwas[gwas.binary == int(wildcards.idx)]
        snps = pd.read_csv(input.snps, sep='\t')
        merged = pd.merge(gwas, snps, on="SNP")

        all_merged = merged.drop(columns=["effect", "binary"]).drop_duplicates()
        all_merged.to_csv(output.all, sep='\t', index=False)

        negative = merged[((merged.effect < 0) & (merged.alt == merged.derived)) | ((merged.effect > 0) & (merged.alt == merged.ancestral))] 
        positive = merged[((merged.effect > 0) & (merged.alt == merged.derived)) | ((merged.effect < 0) & (merged.alt == merged.ancestral))] 

        negative = negative.drop(columns=["binary", "effect"]).drop_duplicates()
        negative.to_csv(output.negative, sep='\t', index=False)
        positive = positive.drop(columns=["binary", "effect"]).drop_duplicates()
        positive.to_csv(output.positive, sep='\t', index=False)

def pick_control(wildcards):
    if ("tive" in wildcards.pref) or ("trait" in wildcards.pref) or ("bin" in wildcards.pref):
        return(expand("data/pooled_all_control_DAF.txt"))
    else:
        return(expand("data/pooled_{pref}_control_DAF.txt", pref=wildcards.pref))

rule empirical_probabilities:
    input:
        gwas="data/pooled_{pref}_gwas_DAF.txt",
        control=pick_control
    output:
        "data/distributions/empirical/{pref}_gwas_pvalues.txt",
        "data/distributions/empirical/{pref}_gwas_summary.txt",
        "data/distributions/empirical/{pref}_matched_summary.txt",
        expand("data/distributions/empirical/{{pref}}_gwas_{grp}_cfs.txt",
               grp=["CHB", "YRI"]),
        expand("data/distributions/empirical/{{pref}}_matched_{grp}_cfs.txt",
               grp=["CHB", "YRI"])
    params:
        prefix="data/distributions/empirical/{pref}"
    shell:
        """
        mkdir -p data/distributions/empirical
        conda activate pandas-env
        python scripts/empirical_cfs.py \
            --gwas {input.gwas} \
            --control {input.control} \
            --out {params.prefix}
        conda deactivate
        """

rule ancestral_sfs:
    output:
        "data/distributions/ancestral/{sim}_h{h_coeff}_s{s_coeff}_ancestral_sfs_count_probs.npy"
    params:
        dir="data/distributions/ancestral"
    shell:
        """
        mkdir -p {params.dir}
        conda activate sm-debug
        python scripts/ancestral_sfs.py \
            --s_coeff={wildcards.s_coeff} \
            --h_coeff={wildcards.h_coeff} \
            --sim {wildcards.sim} \
            --out {output}
        conda deactivate
        """

def find_sims(wildcards):
    if wildcards.model == "jouganous_wo_migration":
        n_sims = 2000
    else:
        n_sims = 1000
    return(expand("data/simulations/{model}/{s}/ancestor_{freq}/output_{idx}.txt", 
                  model=wildcards.model,
                  s=wildcards.s,
                  freq=wildcards.freq,
                  idx=range(n_sims)))

rule concatenate_simulations:
    input:
        find_sims 
    output:
        "data/simulations/{model}/{s}/ancestral_freq{freq}.txt"
    shell:
        """
        head -n 1 -q {input[0]} > {output}
        tail -n +2 -q {input} >> {output}
        """

rule ooa_probabilities:
    input:
        sims=expand("data/simulations/{{sim}}/h{{h}}_s{{s}}/ancestral_freq{freq}.txt",
                    freq=["%.2f" % (i/100) for i in range(1, 100)] + ["0.00003", "0.0001", "0.0005", "0.001", "0.002", "0.005", "0.015", "0.025"]),
        ancestral=ancient("data/distributions/ancestral/{sim}_h{h}_s{s}_ancestral_sfs_count_probs.npy")
    output:
        expand("data/distributions/{{sim}}/h{{h}}_s{{s}}_{dist}.txt",
               dist=["yri_conditional_ancestor",
                     "ceu_conditional_ancestor",
                     "ancestor_conditional_ceu",
                     "ancestor_conditional_yri",
                     "yri_conditional_ceu"])
    params:
        prefix="data/distributions/{sim}/h{h}_s{s}_"
    shell:
        """
        conda activate py3
        python scripts/ooA_simulation_cfs.py \
            --simulated_data {input.sims} \
            --ancestral_sfs {input.ancestral} \
            --out {params.prefix}
        conda deactivate
        """

rule two_pop_probabilities:
    input:
        p1_sims=expand("data/simulations/generation{{t}}_ancestral{{ancestral}}/modern{{p1_ne}}_growth{{p1_g}}/h{{h}}_s{{s}}/ancestral_freq{freq}.txt",
                    freq=["%.2f" % (i/100) for i in range(1, 100)] + ["0.0001", "0.0005", "0.001", "0.002", "0.005", "0.015", "0.025"]),
        p2_sims=expand("data/simulations/generation{{t}}_ancestral{{ancestral}}/modern{{p2_ne}}_growth{{p2_g}}/h{{h}}_s{{s}}/ancestral_freq{freq}.txt",
                    freq=["%.2f" % (i/100) for i in range(1, 100)] + ["0.0001", "0.0005", "0.001", "0.002", "0.005", "0.015", "0.025"]),
        ancestral=ancient("data/distributions/ancestral/ancestral{ancestral}_h{h}_s{s}_ancestral_sfs_count_probs.npy")
    output:
        expand("data/distributions/generation{{t}}_ancestral{{ancestral}}.p1_ne{{p1_ne}}_growth{{p1_g}}.p2_ne{{p2_ne}}_growth{{p2_g}}/h{{h}}_s{{s}}_{dist}.txt",
               dist=["p1_conditional_ancestor",
                     "p2_conditional_ancestor",
                     "ancestor_conditional_p2",
                     "ancestor_conditional_p1",
                     "p1_conditional_p2",
                     "p2_conditional_p1"])
    params:
        prefix="data/distributions/generation{t}_ancestral{ancestral}.p1_ne{p1_ne}_growth{p1_g}.p2_ne{p2_ne}_growth{p2_g}/h{h}_s{s}_",
        dir="data/distributions/generation{t}_ancestral{ancestral}.p1_ne{p1_ne}_growth{p1_g}.p2_ne{p2_ne}_growth{p2_g}",
        p1_final_ne=lambda wildcards: float(wildcards.p1_ne) * (1 + float(wildcards.p1_g)) ** (int(wildcards.t[:-1]) * 1000),
        p2_final_ne=lambda wildcards: float(wildcards.p2_ne) * (1 + float(wildcards.p2_g)) ** (int(wildcards.t[:-1]) * 1000)
    shell:
        """
        mkdir -p {params.dir}
        conda activate py3
        python scripts/two_pop_simulation_cfs.py \
            --p1_sims {input.p1_sims} \
            --p2_sims {input.p2_sims} \
            --n_p1 {params.p1_final_ne} \
            --n_p2 {params.p2_final_ne} \
            --ancestral_sfs {input.ancestral} \
            --out_prefix {params.prefix}
        conda deactivate
        """

rule generate_dtwf_transition_no_growth:
    input:
        ancient("data/distributions/ancestral/ancestral{ancestral}_h{h}_s{s}_ancestral_sfs_count_probs.npy")
    output:
        "data/dtwf_chunks/chunk{chunk}_ancestral{ancestral}_modern{modern}_growth0_gen{gen}_h{h}_s{s}.npy"
    params:
        dir="data/dtwf_chunks"
    shell:
        """
        mkdir -p {params.dir}
        conda activate sm-debug
        python scripts/dtwf_transition.py \
            --ancestral {input} \
            --n_ancestral {wildcards.ancestral} \
            --n_modern {wildcards.modern} \
            --n_gen {wildcards.gen} \
            --chunk {wildcards.chunk} \
            --h_coeff={wildcards.h} \
            --s_coeff={wildcards.s} \
            --out {output}
        conda deactivate
        """

rule generate_dtwf_transition:
    input:
        ancient("data/distributions/ancestral/ancestral{ancestral}_h{h}_s{s}_ancestral_sfs_count_probs.npy")
    output:
        "data/dtwf_chunks/chunk{chunk}_ancestral{ancestral}_modern{modern}_growth0.00{growth}_gen{gen}_h{h}_s{s}.npy"
    params:
        dir="data/dtwf_chunks"
    shell:
        """
        mkdir -p {params.dir}
        conda activate sm-debug
        python scripts/dtwf_transition_w_growth.py \
            --ancestral {input} \
            --n_ancestral {wildcards.ancestral} \
            --n_modern {wildcards.modern} \
            --growth 0.00{wildcards.growth} \
            --n_gen {wildcards.gen} \
            --chunk {wildcards.chunk} \
            --h_coeff={wildcards.h} \
            --s_coeff={wildcards.s} \
            --out {output}
        conda deactivate
        """

rule merge_dtwf_transitions:
    input:
        expand("data/dtwf_chunks/chunk{chunk}_ancestral{{ancestral}}_modern{{modern}}_growth{{growth}}_gen{{gen}}_h{{h}}_s{{s}}.npy", 
               chunk=range(21))
    output:
        "data/distributions/dtwf-raw/ancestral{ancestral}_modern{modern}_growth{growth}_gen{gen}_h{h}_s{s}.npy"
    params:
        input_prefix="data/dtwf_chunks/chunk",
        input_suffix="_ancestral{ancestral}_modern{modern}_growth{growth}_gen{gen}_h{h}_s{s}.npy"
    shell:
        """
        mkdir -p data/distributions/dtwf-raw
        conda activate sm-debug
        python scripts/munge_dtwf_transitions.py \
            --prefix {params.input_prefix} \
            --suffix {params.input_suffix} \
            --n_ancestral {wildcards.ancestral} \
            --h_coeff={wildcards.h} \
            --s_coeff={wildcards.s} \
            --out {output}
        conda deactivate
        """

rule dtwf_two_pop_probabilities:
    input:
        anc=ancient("data/distributions/ancestral/ancestral{ancestral}_h{h}_s{s}_ancestral_sfs_count_probs.npy"),
        p1=ancient("data/distributions/dtwf-raw/ancestral{ancestral}_modern{p1}_growth{p1_growth}_gen{gen}_h{h}_s{s}.npy"),
        p2=ancient("data/distributions/dtwf-raw/ancestral{ancestral}_modern{p2}_growth{p2_growth}_gen{gen}_h{h}_s{s}.npy")
    output:
        expand("data/distributions/dtwf-generation{{gen}}_ancestral{{ancestral}}.p1_ne{{p1}}_growth{{p1_growth}}.p2_ne{{p2}}_growth{{p2_growth}}/h{{h}}_s{{s}}_{dist}.txt",
               dist=["p1_conditional_ancestor",
                     "p2_conditional_ancestor",
                     "ancestor_conditional_p2",
                     "ancestor_conditional_p1",
                     "p1_conditional_p2",
                     "p2_conditional_p1"])
    params:
        dir="data/distributions/dtwf-generation{gen}_ancestral{ancestral}.p1_ne{p1}_growth{p1_growth}.p2_ne{p2}_growth{p2_growth}",
        prefix="data/distributions/dtwf-generation{gen}_ancestral{ancestral}.p1_ne{p1}_growth{p1_growth}.p2_ne{p2}_growth{p2_growth}/h{h}_s{s}_"
    shell:
        """
        mkdir -p {params.dir}
        conda activate py3
        python scripts/two_pop_dtwf_cfs.py \
            --ancestral {input.anc} \
            --p1_forward {input.p1} \
            --p2_forward {input.p2} \
            --out_prefix {params.prefix}
        conda deactivate
        """

rule dtwf_jouganous_wo_migration_probabilities:
    input:
        anc=ancient("data/distributions/ancestral/jouganous_wo_migration_h{h}_s{s}_ancestral_sfs_count_probs.npy"),
        anc_yri="data/distributions/dtwf-raw/ancestral23721_modern23721_growth0_gen4103_h{h}_s{s}.npy",
        anc_ooa="data/distributions/dtwf-raw/ancestral23721_modern2831_growth0_gen2517_h{h}_s{s}.npy",
        ooa_ceu="data/distributions/dtwf-raw/ancestral2831_modern2512_growth0.0016_gen1586_h{h}_s{s}.npy",
        ooa_chb="data/distributions/dtwf-raw/ancestral2831_modern1019_growth0.0026_gen1586_h{h}_s{s}.npy",
    output:
        expand("data/distributions/dtwf-jouganous_wo_migration/h{{h}}_s{{s}}_{dist}.txt",
               dist=["yri_conditional_ancestor",
                     "ceu_conditional_ancestor",
                     "ancestor_conditional_ceu",
                     "ancestor_conditional_yri",
                     "yri_conditional_ceu"])
    params:
        dir="data/distributions/dtwf-jouganous_wo_migration",
        prefix="data/distributions/dtwf-jouganous_wo_migration/h{h}_s{s}_"
    shell:
        """
        mkdir -p {params.dir}
        conda activate py3
        python scripts/ooA_dtwf_cfs.py \
            --ancestral {input.anc} \
            --ancestor_yri {input.anc_yri} \
            --ancestor_ooa {input.anc_ooa} \
            --ooa_ceu {input.ooa_ceu} \
            --ooa_chb {input.ooa_chb} \
            --out_prefix {params.prefix}
        conda deactivate
        """
