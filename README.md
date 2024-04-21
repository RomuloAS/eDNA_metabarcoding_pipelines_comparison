# eDNA Metabarcoding pipelines comparison

1. First of all, download the raw FASTQ files from the NCBI Sequence Read Archive under accession number PRJNA611963.
2. Then, before the pipelines execution, adapters in the 3’ end of the read need to be removed using Cutadapt version 1.18.
    - python Remove_Adapter.py FASTQ_files_folder/ Barcode_Adapters_Information/