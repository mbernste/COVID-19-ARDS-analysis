##

configfile: 'config.json'

rule all:
    input:
        '{}/COVID_ICU_VS_NONICU_DE_GSEA.json'.format(config['output_location']),
        '{}/NONCOVID_ICU_VS_NONICU_DE_GSEA.json'.format(config['output_location']),
        '{}/ICU_COVID_VS_NONCOVID_DE_GSEA.json'.format(config['output_location']),
        '{}/NONICU_COVID_VS_NONCOVID_DE_GSEA.json'.format(config['output_location'])                
               
 
###############################################################
#   Data preprocessing
###############################################################
rule remove_duplicate_patient:
    input:
        ec='{}/genes.ec.no_hg.tab'.format(
            config['raw_data_location']
        ),
        tpm='{}/genes.tpm.no_hg.tab'.format(
            config['raw_data_location']
        )
    output:
        ec='{}/genes.ec.no_hg.no_C054.tab'.format(
            config['intermediate_data_location']
        ),
        tpm='{}/genes.tpm.no_hg.no_C054.tab'.format(
            config['intermediate_data_location']
        )
    run:
        cmds=[
            'python remove_duplicate.py {input.ec} {output.ec}',
            'python remove_duplicate.py {input.tpm} {output.tpm}'
        ]
        for cmd in cmds:
            shell('echo "{}"'.format(cmd))
            shell(cmd)

###############################################################
#   Run GSVA
###############################################################
rule run_gsva_GO_biological_process:
    input:
        '{}/genes.ec.no_hg.no_C054.tab'.format(
            config['intermediate_data_location']
        )
    output:
        '{}/gsva.no_hg.no_C054.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python {gsva}/run_gsva.py {{input}} {gsva}/gene_sets/c5.bp.v7.1.symbols.gmt -o {{output}}'.format(
            gsva=config['gsva_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule run_gsva_GO_molecular_function:
    input:
        '{}/genes.ec.no_hg.no_C054.tab'.format(
            config['intermediate_data_location']
        )
    output:
        '{}/gsva_go_mf.no_hg.no_C054.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python {gsva}/run_gsva.py {{input}} {gsva}/gene_sets/c5.mf.v7.1.symbols.gmt -o {{output}}'.format(
            gsva=config['gsva_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)


rule run_gsva_KEGG:
    input:
        '{}/genes.ec.no_hg.no_C054.tab'.format(
            config['intermediate_data_location']
        )
    output:
        '{}/gsva_kegg.no_hg.no_C054.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python {gsva}/run_gsva.py {{input}} {gsva}/gene_sets/c2.cp.kegg.v7.1.symbols.gmt -o {{output}}'.format(
            gsva=config['gsva_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule run_gsva_canonical:
    input:
        '{}/genes.ec.no_hg.no_C054.tab'.format(
            config['intermediate_data_location']
        )
    output:
        '{}/gsva_canonical_pathways.no_hg.no_C054.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python {gsva}/run_gsva.py {{input}} {gsva}/gene_sets/c2.cp.v7.1.symbols.gmt -o {{output}}'.format(
            gsva=config['gsva_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule run_gsva_hallmark:
    input:
        '{}/genes.ec.no_hg.no_C054.tab'.format(
            config['intermediate_data_location']
        )
    output:
        '{}/gsva_hallmark.no_hg.no_C054.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python {gsva}/run_gsva.py {{input}} {gsva}/gene_sets/h.all.v7.1.symbols.gmt -o {{output}}'.format(
            gsva=config['gsva_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule run_gsva_immunological:
    input:
        '{}/genes.ec.no_hg.no_C054.tab'.format(
            config['intermediate_data_location']
        )
    output:
        '{}/gsva_immunological.no_hg.no_C054.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python {gsva}/run_gsva.py {{input}} {gsva}/gene_sets/c7.all.v7.1.symbols.gmt -o {{output}}'.format(
            gsva=config['gsva_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)

###############################################################
#   Dimensionality reduction plots
###############################################################
# Run dimension reduction on only COVID samples
rule run_dimension_reduction_expression_COVID:
    input:
        tpm='{}/genes.tpm.no_hg.no_C054.tab'.format(
            config['intermediate_data_location']
        ),
        meta='{}/DCD_v2.tsv'.format(
            config['raw_data_location']
        )
    output:
        '{}/no_hg.no_C054.only_COVID.PCA_hospital_free.pdf'.format(
            config['output_location']
        )
    run:
        cmd='python dimensionality_reduction.py -l -s COVID {{input.tpm}} {{input.meta}} -o {}/no_hg.no_C054.only_COVID'.format(config['output_location'])
        shell('echo "{}"'.format(cmd))
        shell(cmd)


# Run dimension reduction on only COVID samples
rule run_dimension_reduction_gsva_COVID:
    input:
        gsva='{}/gsva.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        meta='{}/DCD.tab'.format(
            config['raw_data_location']
        )
    output:
        '{}/no_hg.no_C054.gsva_GO_biological_process.only_COVID.PCA_hospital_free.pdf'.format(
            config['output_location']
        )
    run:
        cmd='python dimensionality_reduction.py -l -s COVID {{input.gsva}} {{input.meta}} -o {}/no_hg.no_C054.gsva_GO_biological_process.only_COVID'.format(config['output_location'])
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule run_dimension_reduction_gsva_GO_MF_COVID:
    input:
        gsva='{}/gsva_go_mf.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        meta='{}/DCD.tab'.format(
            config['raw_data_location']
        )
    output:
        '{}/no_hg.no_C054.gsva_GO_molecular_function.only_COVID.PCA_hospital_free.pdf'.format(
            config['output_location']
        )
    run:
        cmd='python dimensionality_reduction.py -l -s COVID {{input.gsva}} {{input.meta}} -o {}/no_hg.no_C054.gsva_GO_molecular_function.only_COVID'.format(config['output_location'])
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule run_dimension_reduction_gsva_kegg_COVID:
    input:
        gsva='{}/gsva_kegg.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        meta='{}/DCD.tab'.format(
            config['raw_data_location']
        )
    output:
        '{}/no_hg.no_C054.gsva_KEGG.only_COVID.PCA_hospital_free.pdf'.format(
            config['output_location']
        )
    run:
        cmd='python dimensionality_reduction.py -l -s COVID {{input.gsva}} {{input.meta}} -o {}/no_hg.no_C054.gsva_KEGG.only_COVID'.format(config['output_location'])
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule run_ElasticNet_regression_GSVA_COVID:
    input:
        gsva='{}/gsva.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        meta='{}/DCD_v2.tsv'.format(
            config['raw_data_location']
        )
    output:
        '{}/EN_GSVA_GO_biological_processes.only_COVID.kept_features.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python elasticnet_hospital_free.py {{input.gsva}} {{input.meta}} -s COVID -o {out}/EN_GSVA_GO_biological_processes.only_COVID'.format(
            out=config['output_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule run_ElasticNet_regression_expression_COVID:
    input:
        tpm='{}/genes.tpm.no_hg.tab'.format(
            config['raw_data_location']
        ),
        meta='{}/DCD_v2.tsv'.format(
            config['raw_data_location']
        )
    output:
        '{}/EN_genes.only_COVID.kept_features.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python elasticnet_hospital_free.py {{input.tpm}} {{input.meta}} -l -s COVID -o {out}/EN_genes.only_COVID'.format(
            out=config['output_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule gsea_elastic_net_expression_covid:
    input:
        pos='{}/EN_genes.only_COVID.positive_features.tsv'.format(
            config['output_location']
        ),
        neg='{}/EN_genes.only_COVID.negative_features.tsv'.format(
            config['output_location']
        )
    output:
        pos='{}/EN_genes.only_COVID.GSEA_positive.tsv'.format(
            config['output_location']
        ),
        neg='{}/EN_genes.only_COVID.GSEA_negative.tsv'.format(
            config['output_location']
        )
    run:
        cmds=[
            'python gsea_regression_output.py {input.neg} -o {output.neg}',
            'python gsea_regression_output.py {input.pos} -o {output.pos}'
        ]
        for c in cmds:
            shell('echo "{}"'.format(c))
            shell(c)

NEG_HALLMARK_SETS = [
    'HALLMARK_COAGULATION',
    'HALLMARK_INTERFERON_GAMMA_RESPONSE',
    'HALLMARK_INTERFERON_ALPHA_RESPONSE',
    'HALLMARK_INFLAMMATORY_RESPONSE'
]
rule elastic_net_negative_hallmark_heatmaps:
    output:
        [
            '{out}/COVID_elastic_net_heatmaps/{g_set}.heatmap.pdf'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in NEG_HALLMARK_SETS
        ]
    run:
        cmds = [
            'python heatmap_from_regression.py ../gsvapy/gene_sets/h.all.v7.1.symbols.gmt {g_set} {out}/EN_genes.only_COVID.positive_features.tsv {out}/EN_genes.only_COVID.negative_features.tsv AHNMJYDMXX_rsem/DCD_v2.tsv AHNMJYDMXX_rsem/genes.tpm.no_hg.no_C054.tab -o {out}/COVID_elastic_net_heatmaps -d'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in NEG_HALLMARK_SETS
        ]
        for c in cmds:
            shell('echo "{}"'.format(c))
            shell(c)

NEG_GO_SETS = [
    'GO_HUMORAL_IMMUNE_RESPONSE',
    'GO_INNATE_IMMUNE_RESPONSE',
    'GO_INFLAMMATORY_RESPONSE',
    'GO_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND',
    'GO_WOUND_HEALING',
    'GO_REGULATION_OF_RESPONSE_TO_CYTOKINE_STIMULUS',
    'GO_PLATELET_DEGRANULATION',
    'GO_MODULATION_BY_HOST_OF_VIRAL_PROCESS'
]
rule elastic_net_negative_go_heatmaps:
    output:
        [
            '{out}/COVID_elastic_net_heatmaps/{g_set}.heatmap.pdf'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in NEG_GO_SETS
        ]
    run:
        cmds = [
            'python heatmap_from_regression.py ../gsvapy/gene_sets/c5.bp.v7.1.symbols.gmt {g_set} {out}/EN_genes.only_COVID.positive_features.tsv {out}/EN_genes.only_COVID.negative_features.tsv AHNMJYDMXX_rsem/DCD_v2.tsv AHNMJYDMXX_rsem/genes.tpm.no_hg.no_C054.tab -o {out}/COVID_elastic_net_heatmaps -d'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in NEG_GO_SETS
        ]
        for c in cmds:
            shell('echo "{}"'.format(c))
            shell(c)

# Positive genes
POS_REACTOME_SETS = [
    'REACTOME_AMYLOID_FIBER_FORMATION',
    'REACTOME_TELOMERE_MAINTENANCE'
]
rule elastic_net_positive_cannonical_heatmaps:
    output:
        [
            '{out}/COVID_elastic_net_heatmaps/{g_set}.heatmap.pdf'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in POS_REACTOME_SETS
        ]
    run:
        cmds = [
            'python heatmap_from_regression.py ../gsvapy/gene_sets/c2.cp.v7.1.symbols.gmt {g_set} {out}/EN_genes.only_COVID.positive_features.tsv {out}/EN_genes.only_COVID.negative_features.tsv AHNMJYDMXX_rsem/DCD_v2.tsv AHNMJYDMXX_rsem/genes.tpm.no_hg.no_C054.tab -o {out}/COVID_elastic_net_heatmaps -u'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in POS_REACTOME_SETS
        ]
        for c in cmds:
            shell('echo "{}"'.format(c))
            shell(c)



#######################################################
#   Create barplots for all enrichment scores
#######################################################
rule gsva_plots_EN_genes:
    input:
        gsva_go_bp='{}/gsva.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        gsva_go_mf='{}/gsva_go_mf.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        gsva_kegg='{}/gsva_kegg.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        gsva_canon='{}/gsva_canonical_pathways.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        gsva_hall='{}/gsva_hallmark.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),    
        gsva_immun='{}/gsva_immunological.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        meta='{}/DCD_v2.tsv'.format(
            config['raw_data_location']
        ),
        gsea='{}/EN_genes.only_COVID.GSEA.json'.format(
            config['output_location']
        )
    run:
        cmds=[
            'mkdir -p {}/GSVA_figures'.format(config['output_location']),
            'python plot_gsva_scores.py {{input.gsva_go_bp}} {{input.meta}} {{input.gsea}} -o {out}/GSVA_figures'.format(
                out=config['output_location']
            ),
            'python plot_gsva_scores.py {{input.gsva_go_mf}} {{input.meta}} {{input.gsea}} -o {out}/GSVA_figures'.format(
                out=config['output_location']
            ),
            'python plot_gsva_scores.py {{input.gsva_kegg}} {{input.meta}} {{input.gsea}} -o {out}/GSVA_figures'.format(
                out=config['output_location']
            ),
            'python plot_gsva_scores.py {{input.gsva_canon}} {{input.meta}} {{input.gsea}} -o {out}/GSVA_figures'.format(
                out=config['output_location']
            ),
            'python plot_gsva_scores.py {{input.gsva_hall}} {{input.meta}} {{input.gsea}} -o {out}/GSVA_figures'.format(
                out=config['output_location']
            ),
            'python plot_gsva_scores.py {{input.gsva_immun}} {{input.meta}} {{input.gsea}} -o {out}/GSVA_figures'.format(
                out=config['output_location']
            )
        ]
        for cmd in cmds:
            shell('echo "{}"'.format(cmd))
            shell(cmd)

rule gsva_plots_EN_GSVA:
    input:
        gsva_go_bp='{}/gsva.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        meta='{}/DCD_v2.tsv'.format(
            config['raw_data_location']
        ),
        gsea='{}/EN_GSVA_GO_biological_processes.only_COVID.kept_features.tsv'.format(
            config['output_location']
        )
    run:
        cmds=[
            'mkdir -p {}/GSVA_figures'.format(config['output_location']),
            'python plot_gsva_scores.py {{input.gsva_go_bp}} {{input.meta}} {{input.gsea}} -o {out}/GSVA_figures'.format(
                out=config['output_location']
            )
        ]
        for cmd in cmds:
            shell('echo "{}"'.format(cmd))
            shell(cmd)
            

#######################################################
#   FET on DE genes
#######################################################

# COVID ICU vs. non-ICU
rule covid_ICU_DE:
    output:
        up='{}/COVID_ICU_VS_NONICU.DE_up_GSEA.tsv'.format(config['output_location']),
        down='{}/COVID_ICU_VS_NONICU.DE_down_GSEA.tsv'.format(config['output_location'])
    run:
        cmds=[
            'python gsea.py Albany/ebseq_COVID.ICU.v.NO_ICU/Up.Genes.pp99.txt Albany/ebseq_COVID.ICU.v.NO_ICU/Up.Genes.pp99.Normed.tsv AHNMJYDMXX_rsem/genes.tpm.no_hg.no_C054.tab -o {output.up}',
            'python gsea.py Albany/ebseq_COVID.ICU.v.NO_ICU/Down.Genes.pp99.txt Albany/ebseq_COVID.ICU.v.NO_ICU/Down.Genes.pp99.Normed.tsv AHNMJYDMXX_rsem/genes.tpm.no_hg.no_C054.tab -o {output.down}'
        ]
        for cmd in cmds:
            shell(cmd)

#HERE!

# non-COVID ICU vs. non-ICU
rule noncovid_ICU_DE:
    output:
        '{}/NONCOVID_ICU_VS_NONICU_DE_GSEA.tsv'.format(config['output_location'])
    run:
        cmd='python gsea.py NONCOVID_ICU_ebseq/Up.Genes.pp95.txt,NONCOVID_ICU_ebseq/Down.Genes.pp95.txt -o {output}'
        shell('echo "{}"'.format(cmd))
        shell(cmd)

# ICU COVID vs. non-COVID
rule icu_covid_vs_noncovid_DE:
    output:
        '{}/ICU_COVID_VS_NONCOVID_DE_GSEA.tsv'.format(config['output_location'])
    run:
        cmd='python gsea.py ebseq_ICU.COVID.v.NON_COVID/Up.Genes.pp95.txt,ebseq_ICU.COVID.v.NON_COVID/Down.Genes.pp95.txt -o {output}'
        shell('echo "{}"'.format(cmd))
        shell(cmd)

# non-ICU COVID vs. non-COVID
rule nonicu_covid_vs_noncovid_DE:
    output:
        '{}/NONICU_COVID_VS_NONCOVID_DE_GSEA.tsv'.format(config['output_location'])
    run:
        cmd='python gsea.py ebseq_NO_ICU.COVID.v.NON_COVID/Up.Genes.pp95.txt,ebseq_NO_ICU.COVID.v.NON_COVID/Down.Genes.pp95.txt -o {output}'
        shell('echo "{}"'.format(cmd))
        shell(cmd)

# ICU COVID vs. Englart ARDS
rule icu_covid_vs_englart_ards_down:
    output:
        up='{}/ICU_COVID_VS_Englart_ARDS_up_DE_GSEA.tsv'.format(config['output_location']),
        down='{}/ICU_COVID_VS_Englart_ARDS_down_DE_GSEA.tsv'.format(config['output_location'])
    run:
        cmds=[
            'python gsea.py ebseq_ARDS.COVID.v.NO_HSCT/Up.Genes.pp99.txt ebseq_ARDS.COVID.v.NO_HSCT/Up.Genes.pp99.Normed.tsv -o {output.up}',
            'python gsea.py ebseq_ARDS.COVID.v.NO_HSCT/Down.Genes.pp99.txt ebseq_ARDS.COVID.v.NO_HSCT/Down.Genes.pp99.Normed.tsv -o {output.down}'
        ]
        for cmd in cmds:
            shell(cmd)

#################################################################
#   Heatmaps
#################################################################

# Up in Englart
UP_ENGLART_REACTOME_SETS = [
    'REACTOME_INTERFERON_ALPHA_BETA_SIGNALING',
    'REACTOME_SIRT1_NEGATIVELY_REGULATES_RRNA_EXPRESSION',
    'KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS',
    'REACTOME_PRE_NOTCH_EXPRESSION_AND_PROCESSING',
    'REACTOME_AMYLOID_FIBER_FORMATION',
    'REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE',
    'REACTOME_ESTROGEN_DEPENDENT_GENE_EXPRESSION',
    'REACTOME_METALLOPROTEASE_DUBS',
    'REACTOME_RHO_GTPASES_ACTIVATE_PKNS',
    'REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP',
    'REACTOME_TELOMERE_MAINTENANCE'
]
rule englert_up_Reactome_heatmaps:
    output:
        [
            '{out}/ICU_COVID_VS_Englart_ARDS_DE_heatmaps/Up_Englert/{g_set}.heatmap.pdf'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in UP_ENGLART_REACTOME_SETS
        ]
    run:
        cmds = [
            'python heatmap_from_de.py ../gsvapy/gene_sets/c2.cp.v7.1.symbols.gmt {g_set} ebseq_ARDS.COVID.v.NO_HSCT/Up.Genes.pp99.Normed.tsv ebseq_ARDS.COVID.v.NO_HSCT/Down.Genes.pp99.Normed.tsv AHNMJYDMXX_rsem/DCD_v2.tsv englert_meta.tsv AHNMJYDMXX_rsem/genes.tpm.no_hg.no_C054.tab ARDS/g.tpm.tab -o {out}/ICU_COVID_VS_Englart_ARDS_DE_heatmaps/Up_Englert -d'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in UP_ENGLART_REACTOME_SETS
        ]
        for c in cmds:
            shell('echo "{}"'.format(c))
            shell(c)

UP_ENGLART_HALLMARK_SETS = [
    'HALLMARK_INTERFERON_ALPHA_RESPONSE',
    'HALLMARK_INTERFERON_GAMMA_RESPONSE'
]
rule englert_up_hallmark_heatmaps:
    output:
        [
            '{out}/ICU_COVID_VS_Englart_ARDS_DE_heatmaps/Up_Englert/{g_set}.heatmap.pdf'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in UP_ENGLART_HALLMARK_SETS
        ]
    run:
        cmds = [
            'python heatmap_from_de.py ../gsvapy/gene_sets/h.all.v7.1.symbols.gmt {g_set} ebseq_ARDS.COVID.v.NO_HSCT/Up.Genes.pp99.Normed.tsv ebseq_ARDS.COVID.v.NO_HSCT/Down.Genes.pp99.Normed.tsv AHNMJYDMXX_rsem/DCD_v2.tsv englert_meta.tsv AHNMJYDMXX_rsem/genes.tpm.no_hg.no_C054.tab ARDS/g.tpm.tab -o {out}/ICU_COVID_VS_Englart_ARDS_DE_heatmaps/Up_Englert -d'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in UP_ENGLART_HALLMARK_SETS
        ]
        for c in cmds:
            shell('echo "{}"'.format(c))
            shell(c)

# Up in Albany dataset
UP_GO_SETS = [
    'GO_POSITIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY',
    'GO_RESPONSE_TO_LEUKEMIA_INHIBITORY_FACTOR',
    'GO_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY',
    'GO_VIRAL_GENE_EXPRESSION'
]
rule albany_up_GO_heatmaps:
    output:
        [
            '{out}/ICU_COVID_VS_Englart_ARDS_DE_heatmaps/Up_Albany{g_set}.heatmap.pdf'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in UP_GO_SETS
        ]
    run:
        cmds = [
            'python heatmap_from_de.py ../gsvapy/gene_sets/c5.bp.v7.1.symbols.gmt {g_set} ebseq_ARDS.COVID.v.NO_HSCT/Up.Genes.pp99.Normed.tsv ebseq_ARDS.COVID.v.NO_HSCT/Up.Genes.pp99.Normed.tsv AHNMJYDMXX_rsem/DCD_v2.tsv englert_meta.tsv AHNMJYDMXX_rsem/genes.tpm.no_hg.no_C054.tab ARDS/g.tpm.tab -o {out}/ICU_COVID_VS_Englart_ARDS_DE_heatmaps/Up_Albany -u'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in UP_GO_SETS
        ]
        for c in cmds:
            shell('echo "{}"'.format(c))
            shell(c)

UP_CANON_SETS = [
    'REACTOME_ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION',
    'REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS',
    'REACTOME_SIGNALING_BY_ROBO_RECEPTORS',
    'REACTOME_NGF_STIMULATED_TRANSCRIPTION',
    'BIOCARTA_VIP_PATHWAY',
    'REACTOME_TOLL_LIKE_RECEPTOR_CASCADES',
    'REACTOME_NERVOUS_SYSTEM_DEVELOPMENT'
]
rule albany_up_cannonical_heatmaps:
    output:
        [
            '{out}/ICU_COVID_VS_Englart_ARDS_DE_heatmaps/Up_Albany/{g_set}.heatmap.pdf'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in UP_CANON_SETS
        ]
    run:
        cmds = [
            'python heatmap_from_de.py ../gsvapy/gene_sets/c2.cp.v7.1.symbols.gmt {g_set} ebseq_ARDS.COVID.v.NO_HSCT/Up.Genes.pp99.Normed.tsv ebseq_ARDS.COVID.v.NO_HSCT/Up.Genes.pp99.Normed.tsv AHNMJYDMXX_rsem/DCD_v2.tsv englert_meta.tsv AHNMJYDMXX_rsem/genes.tpm.no_hg.no_C054.tab ARDS/g.tpm.tab -o {out}/ICU_COVID_VS_Englart_ARDS_DE_heatmaps/Up_Albany -u'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in UP_CANON_SETS
        ]
        for c in cmds:
            shell('echo "{}"'.format(c))
            shell(c)

UP_HALLMARK_SETS = [
    'HALLMARK_TNFA_SIGNALING_VIA_NFKB'
]
rule albany_up_hallmarks_heatmaps:
    output:
        [
            '{out}/ICU_COVID_VS_Englart_ARDS_DE_heatmaps/Up_Albany/{g_set}.heatmap.pdf'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in UP_HALLMARK_SETS
        ]
    run:
        cmds = [
            'python heatmap_from_de.py ../gsvapy/gene_sets/h.all.v7.1.symbols.gmt {g_set} ebseq_ARDS.COVID.v.NO_HSCT/Up.Genes.pp99.Normed.tsv ebseq_ARDS.COVID.v.NO_HSCT/Up.Genes.pp99.Normed.tsv AHNMJYDMXX_rsem/DCD_v2.tsv englert_meta.tsv AHNMJYDMXX_rsem/genes.tpm.no_hg.no_C054.tab ARDS/g.tpm.tab -o {out}/ICU_COVID_VS_Englart_ARDS_DE_heatmaps/Up_Albany -u'.format(
                out=config['output_location'],
                g_set=g_set
            )
            for g_set in UP_HALLMARK_SETS
        ]
        for c in cmds:
            shell('echo "{}"'.format(c))
            shell(c)


