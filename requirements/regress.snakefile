rule all_regression:
    input:
        
expand('{phenotype}/deeprvat/repeat_{repeat}/results/burden_associations_testing.parquet',
               phenotype=phenotypes, type=['deeprvat'], 
repeat=range(n_repeats)),
        
expand('{phenotype}/deeprvat/repeat_{repeat}/results/burden_associations_replication.parquet',
               phenotype=phenotypes, type=['deeprvat'], 
repeat=range(n_repeats))

rule combine_regression_chunks:
    input:
        testing = 
expand('{{phenotype}}/{{type}}/repeat_{{repeat}}/results/burden_associations_testing_{chunk}.parquet', 
chunk=range(n_regression_chunks)),
        replication = 
expand('{{phenotype}}/{{type}}/repeat_{{repeat}}/results/burden_associations_replication_{chunk}.parquet', 
chunk=range(n_regression_chunks))
    output:
        testing = 
'{phenotype}/{type}/repeat_{repeat}/results/burden_associations_testing.parquet',
        replication = 
'{phenotype}/{type}/repeat_{repeat}/results/burden_associations_replication.parquet'
    threads: 1
    resources:
        mem_mb = 2048,
        load = 2000
    shell:
        ' && '.join([
            conda_check,
            py + 'univariate_burden_regression.py 
combine-regression-results '
            '--model-name repeat_{wildcards.repeat} '
            '{input.testing} '
            '{output.testing}',
            py + 'univariate_burden_regression.py 
combine-regression-results '
            '--model-name repeat_{wildcards.repeat} '
            '{input.replication} '
            '{output.replication}'
        ])

rule regress:
    input:
        # '{phenotype}/config.yaml',
        # '{phenotype}/burdens/genes.npy',
        # '{phenotype}/burdens/burdens.zarr',
        # '{phenotype}/burdens/y.zarr',
        # '{phenotype}/burdens/x.zarr',
        # expand('{{phenotype}}/{{type}}/burdens/chunk{chunk}.finished',
        #        chunk=range(n_chunks))
        # # ,
        # # phenotype=phenotypes
        config = "{phenotype}/deeprvat/hpopt_config.yaml",
        chunks = lambda wildcards: expand(
            ('{{phenotype}}/{{type}}/burdens/chunk{chunk}.' +
             ("finished" if wildcards.phenotype == phenotypes[0] else 
"linked")),
            chunk=range(n_burden_chunks)
        ),
        phenotype_0_chunks =  expand(
            phenotypes[0] + '/{{type}}/burdens/chunk{chunk}.finished',
            chunk=range(n_burden_chunks)
        ),
    output:
        
temp('{phenotype}/{type}/repeat_{repeat}/results/burden_associations_testing_{chunk}.parquet'),
        
temp('{phenotype}/{type}/repeat_{repeat}/results/burden_associations_replication_{chunk}.parquet')
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: 24572 + (attempt - 1) * 4098,
        # mem_mb = 16000,
        load = lambda wildcards, attempt: 24000 + (attempt - 1) * 4000
    shell:
        ' && '.join([
            conda_check,
            py + 'univariate_burden_regression.py regress '
            + debug +
            '--chunk {wildcards.chunk} '
            '--n-chunks ' + str(n_regression_chunks) + ' '
            '--use-bias '
            '--repeat {wildcards.repeat} '
            + do_scoretest +
            '{input.config} '
            '{wildcards.phenotype}/{wildcards.type}/burdens ' #TODO make 
this w/o repeats
            
'{wildcards.phenotype}/{wildcards.type}/repeat_{wildcards.repeat}/results'
        ])

rule all_burdens:
    input:
        [
            (f'{p}/deeprvat/burdens/chunk{c}.' +
             ("finished" if p == phenotypes[0] else "linked"))
            for p in phenotypes
            for c in range(n_burden_chunks)
        ]

rule link_burdens:
    priority: 1
    input:
        # checkpoints = lambda wildcards: [
        #     
f'{phenotypes[]}/deeprvat/repeat_{repeat}/models/best/bag_{bag}.ckpt'
        #     for repeat in range(n_repeats) for bag in range(n_bags)
        # ],
        checkpoints = lambda wildcards: [
            f'models/repeat_{repeat}/best/bag_{bag}.ckpt'
            for repeat in range(n_repeats) for bag in range(n_bags)
        ],
        dataset = '{phenotype}/deeprvat/association_dataset.pkl',
        config = 'models/repeat_0/config.yaml' #TODO make this more 
generic
    output:
        '{phenotype}/deeprvat/burdens/chunk{chunk}.linked'
    # params:
    #     link_burdens = lambda wildcards: (
    #         '' if wildcards.phenotype == phenotypes[0]
    #         else f'--link-burdens 
../../../{phenotypes[0]}/deeprvat/burdens/burdens.zarr '
    #     )
    threads: 8
    resources:
        mem_mb = lambda wildcards: 16384,
        load = lambda wildcards: 16000,
    shell:
        ' && '.join([
            conda_check,
            cuda_visible_devices,
            (py + 'univariate_burden_regression.py compute-burdens '
             + debug +
             ' --n-chunks '+ str(n_burden_chunks) + ' '
             f'--link-burdens 
../../../{phenotypes[0]}/deeprvat/burdens/burdens.zarr '
             '--chunk {wildcards.chunk} '
             '--dataset-file {input.dataset} '
             '{input.config} '
             '{input.checkpoints} '
             '{wildcards.phenotype}/deeprvat/burdens'),
            'touch {output}'
        ])

rule compute_burdens:
    priority: 10
    input:
        # checkpoints = lambda wildcards: [
        #     
f'{phenotypes[]}/deeprvat/repeat_{repeat}/models/best/bag_{bag}.ckpt'
        #     for repeat in range(n_repeats) for bag in range(n_bags)
        # ],
        checkpoints = lambda wildcards: [
            f'models/repeat_{repeat}/best/bag_{bag}.ckpt'
            for repeat in range(n_repeats) for bag in range(n_bags)
        ],
        dataset = '{phenotype}/deeprvat/association_dataset.pkl',
        config = 'models/repeat_0/config.yaml' #TODO make this more 
generic
    output:
        '{phenotype}/deeprvat/burdens/chunk{chunk}.finished'
        # lambda wildcards: temp(
        #     '{phenotype}/deeprvat/burdens/chunk{chunk}.' +
        #     ("finished" if wildcards.phenotype == phenotype[0] else 
"linked")
        # )
    # params:
    #     link_burdens = lambda wildcards: (
    #         '' if wildcards.phenotype == phenotypes[0]
    #         else f'--link-burdens 
../../../{phenotypes[0]}/deeprvat/burdens/burdens.zarr '
    #     )
    threads: 8
    resources:
        mem_mb = 2000000,        # Using this value will tell our modified 
lsf.profile not to set a memory resource
        load = 8000,
        gpus = 1
        # mem_mb = lambda wildcards: (32000 if wildcards.phenotype == 
phenotypes[0]
        #                             else 16384),
        # load = lambda wildcards: (32000 if wildcards.phenotype == 
phenotypes[0]
        #                           else 16000),
    shell:
        ' && '.join([
            conda_check,
            cuda_visible_devices,
            (py + 'univariate_burden_regression.py compute-burdens '
             + debug +
             ' --n-chunks '+ str(n_burden_chunks) + ' '
             '--chunk {wildcards.chunk} '
             '--dataset-file {input.dataset} '
             '{input.config} '
             '{input.checkpoints} '
             '{wildcards.phenotype}/deeprvat/burdens'),
            'touch {output}'
        ])

rule all_association_dataset:
    input:
        expand('{phenotype}/deeprvat/association_dataset.pkl',
               phenotype=phenotypes)

rule association_dataset:
    input:
        config = '{phenotype}/deeprvat/hpopt_config.yaml'
    output:
        '{phenotype}/deeprvat/association_dataset.pkl'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: 32000 * (attempt + 1),
        load = 64000
    shell:
        ' && '.join([
            conda_check,
            (py + 'univariate_burden_regression.py make-dataset '
             + debug +
             '{input.config} '
             '{output}')
        ])

# rule reverse_models:
#     input:
#         checkpoints = lambda wildcards: [
#             
f'{training_phenotypes[wildcards.phenotype]}/deeprvat/repeat_{repeat}/models/best/bag_{bag}.ckpt'
#             for repeat in range(n_repeats) for bag in range(n_bags)
#         ],
#         config = '{phenotype}/deeprvat/repeat_0/config.yaml', #TODO make 
this more generic
#         # (the config from any repeat can be used cause they are all the 
same)
#     output:
#         "{phenotype}/deeprvat/reverse_finished.tmp"
#     threads: 4
#     resources:
#         mem_mb = 20480,
#         load = 20480
#     shell:
#         " && ".join([
#             conda_check,
#             (py + "univariate_burden_regression.py reverse-models "
#              "{input.config} "
#              "{input.checkpoints}"),
#             "touch {output}"
#         ])

rule all_bagging:
    input:
        expand('models/repeat_{repeat}/best/bag_{bag}.ckpt',
               bag=range(n_bags), repeat=range(n_repeats)),
        expand('models/repeat_{repeat}/config.yaml',
               repeat=range(n_repeats))

rule best_bagging_run:
    input:
        expand('models/repeat_{{repeat}}/trial{trial_number}/config.yaml',
               trial_number=range(n_trials)),
        # 
expand('models/repeat_{{repeat}}/trial{trial_number}/finished.tmp',
        #        trial_number=range(n_trials)),
    output:
        checkpoints = 
expand('models/repeat_{{repeat}}/best/bag_{bag}.ckpt',
                             bag=range(n_bags)),
        config = 'models/repeat_{repeat}/config.yaml'
    threads: 1
    resources:
        mem_mb = 2048,
        load = 2000
    shell:
        ' && '.join([
            conda_check,
            py + 'train.py best-bagging-run '
            + debug +
            'models/repeat_{wildcards.repeat} '
            'models/repeat_{wildcards.repeat}/best '
            
'models/repeat_{wildcards.repeat}/hyperparameter_optimization.db '
            '{output.config}'
        ])
