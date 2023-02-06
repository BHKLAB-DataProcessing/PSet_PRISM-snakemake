from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider(
#     access_key_id=config["key"],
#     secret_access_key=config["secret"],
#     host=config["host"],
#     stay_on_remote=False
# )

prefix = config["prefix"]

rule get_prism_pset:
    input:
        prefix + "download/profile.sensitivity.PRISM.rds",
        prefix + "download/raw.sensitivity.prismii_v3.rds",
        prefix + "download/secondary-screen-replicate-collapsed-treatment-info.csv",
        prefix + "download/secondary-screen-cell-line-info.csv",
        prefix + "download/secondary-screen-replicate-collapsed-logfold-change.csv",
        prefix + "download/secondary-screen-dose-response-curve-parameters.csv",
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv"
    output:
        prefix + "PRISM.rds"
    shell:
        """
        Rscript {prefix}scripts/getPRISM.R {prefix}
        """

rule download_annotation:
    output:
        prefix + "download/drugs_with_ids.csv",
        prefix + "download/cell_annotation_all.csv"
    shell:
        """
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/drugs_with_ids.csv' \
            -O {prefix}download/drugs_with_ids.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/cell_annotation_all.csv' \
            -O {prefix}download/cell_annotation_all.csv
        """

rule download_prism_data:
    output:
        prefix + "download/profile.sensitivity.PRISM.rds",
        prefix + "download/raw.sensitivity.prismii_v3.rds",
        prefix + "download/secondary-screen-replicate-collapsed-treatment-info.csv",
        prefix + "download/secondary-screen-cell-line-info.csv",
        prefix + "download/secondary-screen-replicate-collapsed-logfold-change.csv",
        prefix + "download/secondary-screen-dose-response-curve-parameters.csv"
    shell:
        """
        Rscript {prefix}scripts/downloadPRISMData.R {prefix}
        """
