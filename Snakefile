from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(
    access_key_id=config["key"],
    secret_access_key=config["secret"],
    host=config["host"],
    stay_on_remote=False
)

prefix = config["prefix"]

rule get_prism_pset:
    input:
        S3.remote(prefix + "download/profile.sensitivity.PRISM.rds"),
        S3.remote(prefix + "download/raw.sensitivity.prismii_v3.rds"),
        S3.remote(
            prefix + "download/secondary-screen-replicate-collapsed-treatment-info.csv"),
        S3.remote(prefix + "download/secondary-screen-cell-line-info.csv"),
        S3.remote(
            prefix + "download/secondary-screen-replicate-collapsed-logfold-change.csv"),
        S3.remote(
            prefix + "download/secondary-screen-dose-response-curve-parameters.csv"),
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/cell_annotation_all.csv")
    output:
        prefix + "PSet_PRISM.rds"
    shell:
        """
        Rscript scripts/getPRISM.R {prefix}
        """

rule download_annotation:
    output:
        S3.remote(prefix + "download/drugs_with_ids.csv"),
        S3.remote(prefix + "download/cell_annotation_all.csv")
    shell:
        """
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/drugs_with_ids.csv' \
            -O {prefix}download/drugs_with_ids.csv
        wget 'https://github.com/BHKLAB-DataProcessing/Annotations/raw/master/cell_annotation_all.csv' \
            -O {prefix}download/cell_annotation_all.csv
        """

rule download_prism_data:
    output:
        S3.remote(prefix + "download/profile.sensitivity.PRISM.rds"),
        S3.remote(prefix + "download/raw.sensitivity.prismii_v3.rds"),
        S3.remote(
            prefix + "download/secondary-screen-replicate-collapsed-treatment-info.csv"),
        S3.remote(prefix + "download/secondary-screen-cell-line-info.csv"),
        S3.remote(
            prefix + "download/secondary-screen-replicate-collapsed-logfold-change.csv"),
        S3.remote(
            prefix + "download/secondary-screen-dose-response-curve-parameters.csv")
    shell:
        """
        Rscript scripts/downloadPRISMData.R {prefix}
        """
