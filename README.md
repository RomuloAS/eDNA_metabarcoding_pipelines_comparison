# eDNA Metabarcoding pipelines comparison

1. First of all, download the raw FASTQ files from the NCBI Sequence Read Archive under accession number PRJNA611963.
2. Then, before the pipelines execution, adapters in the 3’ end of the read need to be removed using Cutadapt version 1.18.
    ```python
        python Codes/Remove_Adapter.py FASTQ_files_folder/ Barcode_Adapters_Information/
    ```
3. For each pipeline, convert the genbank file to fasta file by executing the command below for each one of them:
   ```python
       python Codes/genbank2Fasta.py Reference_Database/reference_database.gb --pipeline anacapa --rank superkingdom
   ``` 
5. Execute each one of the pipelines:
    - [Anacapa](https://github.com/limey-bean/Anacapa)
      * Follow the official documentation to install any required package.
      * Use the configuration files that can found inside the /Config_Files folder
      * The previously converted reference database in FASTA format and a taxonomic table generated in the conversion process are used to convert the database to Bowtie format by using the command below ([For details on how to create it](https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools/blob/master/Manual_addition_of_reads_to_CRUX.txt)):
          ```console
              bowtie2-build -f 12S_.fasta ../12S_bowtie2_database/12S_bowtie2_index
          ```
      * Then, execute the commands below after substituting the input and output folder from the commands:
          ```console
            bash anacapa_QC_dada2.sh -i path_to_input_data_folder -o path_to_output_data_folder -a nextera -t MiSeq -l -f forward_primers.txt -r reverse_primers.txt -q 20 -m 90 -x 0 -y 0 -e minimum_length_for_the_overlapping_region.txt
            bash anacapa_classifier.sh -o path_to_output_data_folder -d path_to_Anacapa_db -l -b 1 -p 0.85 -n 1000
          ```
    - [Barque](https://github.com/enormandeau/barque)
      * Execute it according to the official documentation and use the configuration that can be found inside the /Config_Files folder.
    - [metaBEAT](https://github.com/HullUni-bioinformatics/metaBEAT)
      * Follow the execution of the jupyter notebook /Codes/metaBEAT_workflow.ipynb that can be found inside the codes folder.
    - [MiFish](https://doi.org/10.5061/dryad.54v2q)
      * Download the scripts from the repository
      * The header of the uc_size_fas_integrator.pl script needs to be modified from /$OTUname/ to /\Q$OTUname\E/ and the uc_size_processor.pl script needs to be modified from /$otuname/ to /\Q$otuname\E/, both to deal with Illumina FASTQ files header.
      * Execute it according to the documentation.
    - SEQme
      * First of all, train the classfifier by executing the command below. Use the Config_Files/SEQme_rRNAClassifier.properties file that can be found inside the /Config_Files folder :
        ```console
            classifier train -o Classifier -s reference\_database.fasta -t reference\_database\ _taxid.txt
        ```
      * Then, execute the /Codes/SEQme_workflow.py file that can be found inside the /Codes folder.

6. Execute the /Codes/Filter_Samples_by_Threshold.R file, which can be found inside the /Codes folder, to remove false positive species assignment where the number of reads assigned fell below 0.1 % considering the sample total of reads.
