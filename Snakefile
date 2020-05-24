##

configfile: 'config.json'

rule remove_duplicate_patient:
    input:
        '{}/genes.ec.no_hg.tab'.format(
            config['raw_data_location']
        )
    output:
        '{}/genes.ec.no_hg.no_C054.tab'.format(
            config['intermediate_data_location']
        )
    run:
        cmd='python remove_duplicate.py {input} {output}'
        shell('echo "{}"'.format(cmd))
        shell(cmd)

rule run_gsva:
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
