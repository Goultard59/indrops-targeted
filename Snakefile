SAMPLES = ["006"]

rule all:
    input:
        expand("/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/sam_fold/AR{sample}_R1.sam", sample=SAMPLES),
        expand('/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/csv_construct/B{sample}_8N_tot.csv', sample=SAMPLES),
        expand('/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/hist/reads_per_BC_{sample}.csv', sample=SAMPLES),
        expand('/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/hist/BC_{sample}.jpg', sample=SAMPLES),
        expand('/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/hist/reads_per_BC_{sample}.csv', sample=SAMPLES),
        expand('/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/matrice/mean_analysis_umi_Fbc_B{sample}_nooutlier.csv', sample=SAMPLES)

rule bowtie:
    input:
        "/home/adrien.dufour/NeuroDev_ADD/00_DATA/SingleCell/Targeted/Novaseq_R2_240720/AR{sample}_R1.fastq"
    params:
        "/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/index/target_280720"
    output:
        "/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/sam_fold/AR{sample}_R1.sam"
    shell:
        "bowtie2 --local -x {params} {input} -S {output} --no-hd"

rule csv_construct:
    input:
        sam="/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/sam_fold/AR{sample}_R1.sam",
        rd2="/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/read_2/AR{sample}_R2_dpBC.csv",
        fast="/home/adrien.dufour/NeuroDev_ADD/00_DATA/SingleCell/Targeted/Novaseq_R2_240720/AR{sample}_R2.fastq"
    output:
        "/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/csv_construct/B{sample}_8N_tot.csv"
    shell:
        "python BC_gene_umi_csv_construction_SAstyle.py -s {input.sam} -r {input.rd2} -o {output} -f {input.fast}"

rule hist:
    input:
        csv="/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/csv_construct/B{sample}_8N_tot.csv",
        barcode='/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/read_2/AR{sample}_R2_dpBC.csv'
    output:
        csv="/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/hist/reads_per_BC_{sample}.csv",
        fig="/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/hist/BC_{sample}.jpg"
    shell:
        "python reads_per_BC.py -c {input.csv} -b {input.barcode} -o {output.csv} -f {output.fig}"

rule matrice_generation:
    input:
        csv="/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/csv_construct/B{sample}_8N_tot.csv",
        poly="/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/csv_construct/B{sample}_8N_tot.csv"
    output:
        csv="/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/matrice/mean_analysis_umi_Fbc_B{sample}_nooutlier.csv",
        matrice="/home/adrien.dufour/NeuroDev_ADD/SingleCell/demultiplex/matrice/Matrice_{sample}F.csv"
    params:
        min=1000,
        max=20000
    shell:
        "python Matrice_generation.py -c {input.csv} -p {input.poly} -o {output.csv} -m {output.matrice} --minimum {params.min} --maximum {params.max}"
