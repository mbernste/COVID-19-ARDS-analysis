##

configfile: 'config.json'

rule all:
    input:
        '{}/COVID_ICU_DE_GSEA.json'.format(config['output_location']),
        '{}/NONCOVID_ICU_DE_GSEA.json'.format(config['output_location']),
        '{}/no_hg.no_C054.only_COVID.PCA_hospital_free.pdf'.format(
            config['output_location']
        ),
        '{}/no_hg.no_C054.only_COVID.tSNE_hospital_free.pdf'.format(
            config['output_location']
        ),
        '{}/no_hg.no_C054.only_COVID.LDA_hospital_free.pdf'.format(
            config['output_location']
        )
        
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

rule run_gsva_GO_biological_function:
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
        meta='{}/DCD.tab'.format(
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


# Run LASSO analysis on GSVA output for only COVID samples
rule run_LASSO_regression_GSVA_COVID:
    input:
        gsva='{}/gsva.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        meta='{}/DCD.tab'.format(
            config['raw_data_location']
        )
    output:
        '{}/LASSO_GSVA_GO_biological_processes.only_COVID.kept_features.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python regression_hospital_free.py {{input.gsva}} {{input.meta}} -s COVID -o {out}/LASSO_GSVA_GO_biological_processes.only_COVID'.format(
            out=config['output_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)


rule run_LASSO_regression_GSVA_KEGG_COVID:
    input:
        gsva='{}/gsva_kegg.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        meta='{}/DCD.tab'.format(
            config['raw_data_location']
        )
    output:
        '{}/LASSO_GSVA_KEGG.only_COVID.kept_features.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python regression_hospital_free.py {{input.gsva}} {{input.meta}} -s COVID -o {out}/LASSO_GSVA_KEGG.only_COVID'.format(
            out=config['output_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)


rule run_LASSO_regression_expression_COVID:
    input:
        tpm='{}/genes.tpm.no_hg.tab'.format(
            config['raw_data_location']
        ),
        meta='{}/DCD.tab'.format(
            config['raw_data_location']
        )
    output:
        '{}/LASSO_genes.only_COVID.kept_features.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python regression_hospital_free.py {{input.tpm}} {{input.meta}} -l -s COVID -o {out}/LASSO_genes.only_COVID'.format(
            out=config['output_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)


# Run LASSO on the GSVA output for all samples
rule run_LASSO_regression_GSVA_all:
    input:
        gsva='{}/gsva.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        meta='{}/DCD.tab'.format(
            config['raw_data_location']
        )
    output:
        '{}/LASSO_GSVA_GO_biological_processes.kept_features.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python regression_hospital_free.py {{input.gsva}} {{input.meta}} -o {out}/LASSO_GSVA_GO_biological_processes'.format(
            out=config['output_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule run_ElasticNet_regression_GSVA_COVID:
    input:
        gsva='{}/gsva.no_hg.no_C054.tsv'.format(
            config['output_location']
        ),
        meta='{}/DCD.tab'.format(
            config['raw_data_location']
        )
    output:
        '{}/EN_GSVA_GO_biological_processes.only_COVID.kept_features.tsv'.format(
            config['output_location']
        )
    run:
        cmd='python elasticnet_hospital_free.py {{input.gsva}} {{input.meta}} -l -s COVID -o {out}/EN_GSVA_GO_biological_processes.only_COVID'.format(
            out=config['output_location']
        )
        shell('echo "{}"'.format(cmd))
        shell(cmd)


rule run_ElasticNet_regression_expression_COVID:
    input:
        tpm='{}/genes.tpm.no_hg.tab'.format(
            config['raw_data_location']
        ),
        meta='{}/DCD.tab'.format(
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
        '{}/EN_genes.only_COVID.kept_features.tsv'.format(
            config['output_location']
        )
    output:
        '{}/EN_genes.only_COVID.GSEA.json'.format(
            config['output_location']
        )
    run:
        cmd='python gsea.py {input} -o {output}'
        shell('echo "{}"'.format(cmd))
        shell(cmd)

#######################################################
#   GSEA on DE genes
#######################################################

# COVID ICU vs. non-ICU
rule covid_ICU_DE:
    output:
        '{}/COVID_ICU_DE_GSEA.json'.format(config['output_location'])
    run:
        cmd='python gsea.py COVID_ICU_ebseq/Up.Genes.pp95.txt,COVID_ICU_ebseq/Down.Genes.pp95.txt -o {output}'

rule noncovid_ICU_DE:
    output:
        '{}/NONCOVID_ICU_DE_GSEA.json'.format(config['output_location'])
    run:
        cmd='python gsea.py NONCOVID_ICU_ebseq/Up.Genes.pp95.txt,NONCOVID_ICU_ebseq/Down.Genes.pp95.txt -o {output}'





