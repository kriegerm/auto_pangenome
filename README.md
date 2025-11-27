# auto_pangenome

This pipeline downloads reference genomes from NCBI and whole genome shotgun sequencing reads from SRA. It will then assemble contigs from the WGS reads, build a pangenome with both the reference genomes and the assembled contigs using anvio, and construct a single-copy core gene tree.

You can set parameters in the `config/project_<project_number>.yaml` file. You can make a new project YAML file and for each separate genome build you'd like to do. An example would be `config/project_001.yaml`, so the project number/identifier for this will be 001. 

All you have to do to run the whole pipeline is:
`./analysis/pipeline_commands.sh <project_number>`

So if you want to run project 001, you can do:
`./analysis/pipeline_commands.sh 001`


Relevant environments are stashed in the `env` folder. You need a couple specific conda environments for this to work, which are `ncbi_datasets`, `wgs_spades`, and `anvio-8`. The script uses a funciton that is defined at the top called `activate_and_export` to activate conda environments and export environment YAMLs to the env folder if they don't already exist. 

The final step is using `analyze_tree.Rmd` to use `ggtree` to build and save your final tree. 