# **DeepSAP Manuscript Shared Data and Code Repository**

This repository contains the data, figures, and scripts used in the **DeepSAP manuscript**. The repository is organized into multiple folders, each corresponding to different figures and analyses.

## **Usage Instructions**
1. **Fine-tuning & evaluation**: Refer to `Figure_2_dataset_for_fine-tuning_evaluation/tune_transformers.html`.
2. **Comapring fine-tuned models**: Refer to `Figure_2_dataset_for_fine-tuning_evaluation/visulize_models.html`.
3. **RNA-seq aligners benchmarking Nextflow**: See `Figure_3_benchmarking_data_results/Nextflow.md`.
4. **RNA-seq aligners benchmarking results**: See `Figure_3_benchmarking_data_results/results.tsv`.
5. **Splice junction analysis**:
   - Baruzzo datasets: `Figure_4_sj_scatter_plots_from_counts_baruzzo_datasets/workflow_generating_scatter_plots_baruzzo.html`
   - Real datasets: `Figure_5_sj_scatter_plots_from_counts_real_datasets/workflow_generating_scatter_plots_real.html`
6. **InDels analysis in Baruzzo datasets**: Refer to `Figure_6_InDels_plots_from_baruzzo_datasets/analysze_INDELs_baruzzo_datasets.html`
7. **Source data for figures**:  Refer to `source_data_for_figures.xlsx`

## **Repository Structure**
Below is an overview of the repository structure and its contents:

### **DeepSAP Figures**
📂 `DeepSAP Manuscript and Figures`
- **DeepSAP.pdf** - DeepSAP manuscript.
- **DeepSAP-main-figures.pdf** - Main figures from the DeepSAP manuscript.
- **DeepSAP-supplemental-figures-tables.pdf** - Supplemental figures and tables included in the manuscript.
- **DeepSAP-main-figures-wot-legends.pdf** - Main figures in a simplified format.
- **DeepSAP-supplemental-figures-wot-legends.pdf** - Supplemental figures in a simplified format.

### **Figure 2 - Dataset for Fine-Tuning & Evaluation**
📂 `Figure_2_dataset_for_fine-tuning_evaluation/`
- Contains datasets and scripts for **fine-tuning and evaluating** DeepSAP transformer models shown in Fig. 2 and Fig. 2S.
- **Key folders and files**:
  - `tune_transformers.html` - HTML Notebook for fine-tuning transformers.
  - `visulize_models.html` - HTML notebook for visualization of models performance.
  - `evaluation_of_transformer_models` - Evaluation results, contains the evaluation.pdf and evaluation.csv of evaluating transformer models plot in Fig. 2a. 
  - `fine_tuning_transformer_models/` - Contains evaluation metrics (metric).pdf and (metric).csv files during fine-tuning per window size used in Fig. 2b
  - `comparing_transformer_models/` - Contains comparison between transfromer models DNABERT1__MS150_vs_DNABERT1__PF150_PF150, DNABERT1__MS150 vs DNABERT1__RefSeq150_RefSeq150 and SpliceAI400 vs DNABERT1__RefSeq400_RefSeq400. Each folder includes 
    - `acceptor_donor_score_mean_by_intron_length/` - compare_scores_by_intron_length.svg file besides acceptor_score_mean_by_intron_bin_(model).csv and donor_score_mean_by_exon_bin_(model).csv which contains the mean of acceptor and donor score per exon bin for each model. 
    - `acceptor_donor_score_mean_by_exon_length/` - compare_scores_by_exon_length.svg file besides acceptor_score_mean_by_exon_bin_(model).csv and donor_score_mean_by_intron_bin_(model).csv which contains the mean of acceptor and donor score per intron bin for each model. 
    - `histogram_of_junction_score_by_biotype/` - Contains histogram .psd and .csv files of splice junction scores per biotype for each model. 
    - `violin_comparison/` - Contains comparison_violin_plots_all.csv.zip and pdf file which represents the violin plots of acceptor and donor score accros biotypes and splicing signals, besides the orgingial df in .csv.zip format used to generate the comparison dataframe. 

### **Figure 3 - Benchmarking Data & Results**
📂 `Figure_3_benchmarking_data_results/`
- Contains **benchmarking results** for Baruzzo and SimBA datasets shown in Fig. 3 and Fig. 3S.
- **Key folders and files**:
  - `All/`, `Human/`, `Malaria/` - Benchmarking .png and .csv files for Human only datasets, Malaria only datasets and Humand & Malaria datasets.
  - `Nextflow.md` - Documentation for Nextflow pipeline.
  - `results.tsv` - Table with benchmarking results produced by the Nextflow pipeline.

### **Figure 4 - Splice Junction Scatter Plots (Baruzzo Datasets)**
📂 `Figure_4_sj_scatter_plots_from_counts_baruzzo_datasets/`
- Contains splice junction scatter plots for **Comparing DeepSAP with other aligners in Baruzzo datasets** shown in Fig. 4 and Fig. 4S.
- **Key folders and files**:
  - `workflow_generating_scatter_plots_baruzzo.html` - Steps for generating scatter plots for Baruzzo Human T1, T2 and T3 datasets.
  - `counts/` - Splice junctions counts from featureCount per Baruzzo dataset.
  - `processed_counts/` - Processed splice junctions counts files.
  - `figures_DeepSAP_vs_other_aligners_sj/` - Figures comparing DeepSAP vs. other aligners.
  - `figures_aligners_vs_ground_truth_sj/` - Figures comparing all aligners vs. ground truth.

### **Figure 5 - Splice Junction Scatter Plots (Real Datasets)**
📂 `Figure_5_sj_scatter_plots_from_counts_real_datasets/`
- Contains splice junction scatter plots for **Comparing DeepSAP with other aligners in real datasets** shown in Fig. 5 and Fig. 5S.
- **Key folders and files**:
  - `workflow_generating_scatter_plots_real.html` - Workflow documentation to generate the plots. 
  - `SRR5280319/` - FeatureCounts results for the dataset SRR5280319 used for analysis. Splice junctions counts are stored in zipped file pHGG_SRR5280319_counts_all_MAPQ0.jcounts.zip (MAPQ0 means here that all reads are allowed to be counted depsite the MAPQ)
  - `SRR6781181/` - FeatureCounts results for the dataset SRR6781181 used for analysis. Splice junctions counts are stored in zipped file sQTL_SRR6781181_counts_all_MAPQ0.jcounts.zip (MAPQ0 means here that all reads are allowed to be counted depsite the MAPQ)

### **Figure 6 - InDels Analysis (Baruzzo Datasets)**
📂 `Figure_6_InDels_plots_from_baruzzo_datasets/`
- Contains analysis of insertions and deletions (InDels).
- **Key folders and filess**:
  - `analyze_INDELs_baruzzo_datasets.html` - HTML report on InDels analysis.
  - `analyze_INDELs_baruzzo/` - Processed results folder contains subfolders per Baruzzo dataset, each contains: 
    - `(dataset)_(aligner)_benchmark.txt` - Benchmarking result file generated by Baruzzo benchmarking script `compare2truth.rb`
    - `(dataset)_(aligner)_comparison.txt.gz` - Comparison of indels between the aligner and ground truth. 

### **Individual Figures**
📂 `Individual_Figures/`
- High-quality **PDF versions** of the manuscript figures.
- **Key files**:
  - `Figure_1.pdf`, ..., `Figure_6.pdf` - Main figures.
  - `Figure_1S.pdf`, ..., `Figure_5S.pdf` - Supplemental figures.
  - `Plain_Main_Figures.pdf` - All main figures in a single file.
  - `Plain_Supplemental_Figures.pdf` - All supplemental figures in a single file.