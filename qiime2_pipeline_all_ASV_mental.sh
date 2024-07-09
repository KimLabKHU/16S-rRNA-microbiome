#!/bin/bash 

if [ $# -ne 4 ]; then
	echo "USAGE: ./qiime2_pipeline_all_ASV.sh DATA_dir, Menifest_file, ClinicData Reference (Silva, Greengene)"
	exit 1
else

data_dir=$1
menifest=$2
metadata=$3
reference=$4

eval "$(conda shell.bash hook)"

conda activate qiime2-new

mkdir $data_dir

echo --------------------------------------------
echo "import"
echo --------------------------------------------

qiime tools import \
	--type 'SampleData[JoinedSequencesWithQuality]' \
	--input-path $menifest \
	--output-path $data_dir/1_data.qza \
	--input-format SingleEndFastqManifestPhred33

# Summarize counts per sample.
qiime demux summarize \
--i-data $data_dir/1_data.qza \
--o-visualization $data_dir/1_data_summary.qzv

echo "==================================="

echo "Trimming index of reads and Denoising with Deblur..."

echo "==================================="

echo "Trimming Forward read length : 20"
echo "Truncated read length : 430"

qiime deblur denoise-16S \
  --i-demultiplexed-seqs $data_dir/1_data.qza \
  --p-left-trim-len 20 \
  --p-trim-length 430 \
  --p-sample-stats \
  --o-representative-sequences $data_dir/2_denoised_rep-seqs.qza \
  --o-table $data_dir/2_denoised_table.qza \
  --o-stats $data_dir/2_denoised_denoising-stats.qza

# summary
qiime feature-table summarize \
	--i-table $data_dir/2_denoised_table.qza \
	--o-visualization $data_dir/2_denoised_table.qzv

qiime feature-table tabulate-seqs \
	--i-data $data_dir/2_denoised_rep-seqs.qza \
	--o-visualization $data_dir/2_denoised_rep-seqs.qzv

qiime deblur visualize-stats \
--i-deblur-stats $data_dir/2_denoised_denoising-stats.qza \
--o-visualization $data_dir/2_denoised_denoising-stats.qzv

echo "================================="

echo "Taxonomic classification..."

echo "================================="

if [ $reference = "Silva" ]; then

REF_DIR=/kimlab_wd/rhdfyd/silva
mkdir $data_dir/Silva

qiime feature-classifier extract-reads \
  --i-sequences ${REF_DIR}/silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-trunc-len 270 \
  --o-reads $data_dir/Silva/silva-138-99-ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $data_dir/Silva/silva-138-99-ref-seqs.qza \
  --i-reference-taxonomy ${REF_DIR}/silva-138-99-tax.qza \
  --o-classifier $data_dir/Silva/silva-138-99-classifier.qza

##test
qiime feature-classifier classify-sklearn \
  --i-reads ${data_dir}/2_denoised_rep-seqs.qza \
  --i-classifier $data_dir/Silva/silva-138-99-classifier.qza \
  --o-classification $data_dir/Silva/3_taxonomy.qza

#visualization
qiime metadata tabulate \
  --m-input-file $data_dir/Silva/3_taxonomy.qza \
  --o-visualization $data_dir/Silva/3_taxonomy.qzv

qiime tools export \
--input-path ${data_dir}/Silva/3_taxonomy.qzv \
--output-path ${data_dir}/Silva

rm -rf ${data_dir}/Silva/css* ${data_dir}/Silva/js ${data_dir}/Silva/q2templateassets ${data_dir}/Silva/index*

for i in 2 3 4 5 6 7
do

mkdir $data_dir/Silva/level-${i}
level_DIR=$data_dir/Silva/level-${i}
#taxonomy collapse
qiime taxa collapse \
	--i-table ${data_dir}/2_denoised_table.qza \
	--i-taxonomy ${data_dir}/Silva/3_taxonomy.qza \
	--p-level ${i} \
	--o-collapsed-table ${level_DIR}/3_1_table_collapse_level_${i}.qza \

#table heatmap
qiime feature-table heatmap \
	--i-table ${level_DIR}/3_1_table_collapse_level_${i}.qza \
	--p-metric braycurtis \
	--p-method centroid \
	--p-cluster both \
	--o-visualization ${level_DIR}/3_2_table_heatmap_collapse_level_${i}.qzv

#barplot
qiime taxa barplot \
  --i-table ${data_dir}/2_denoised_table.qza \
  --i-taxonomy ${data_dir}/Silva/3_taxonomy.qza \
  --m-metadata-file $metadata \
  --o-visualization ${level_DIR}/3_3_taxa_bar_plots_level_${i}.qzv

done

taxonomy_DIR=${data_dir}/Silva

elif [ $reference = "Greengene" ]; then

REF_DIR=/home/data/References/Metagenomics/gg_13_8_otus
mkdir $data_dir/Greengene

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path ${REF_DIR}/rep_set/97_otus.fasta \
  --output-path $data_dir/Greengene/97_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path ${REF_DIR}/taxonomy/97_otu_taxonomy.txt \
  --output-path $data_dir/Greengene/97_ref_taxonomy.qza

qiime feature-classifier extract-reads \
  --i-sequences $data_dir/Greengene/97_otus.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-trunc-len 270 \
  --o-reads $data_dir/Greengene/97_ref_seqs.qza

##train
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads $data_dir/Greengene/97_ref_seqs.qza \
  --i-reference-taxonomy $data_dir/Greengene/97_ref_taxonomy.qza \
  --o-classifier $data_dir/Greengene/Greengene_97_classifier.qza \

##test
qiime feature-classifier classify-sklearn \
  --i-reads ${data_dir}/2_denoised_rep-seqs.qza \
  --i-classifier $data_dir/Greengene/Greengene_97_classifier.qza \
  --o-classification $data_dir/Greengene/3_taxonomy.qza \
  --p-n-jobs 2

#visualization
qiime metadata tabulate \
  --m-input-file $data_dir/Greengene/3_taxonomy.qza \
  --o-visualization $data_dir/Greengene/3_taxonomy.qzv

qiime tools export \
--input-path ${data_dir}/Greengene/3_taxonomy.qzv \
--output-path ${data_dir}/Greengene

rm -rf ${data_dir}/Greengene/css* ${data_dir}/Greengene/js ${data_dir}/Greengene/q2templateassets ${data_dir}/Greengene/index*

for i in 2 3 4 5 6 7
do

mkdir $data_dir/Greengene/level-${i}
level_DIR=$data_dir/Greengene/level-${i}
#taxonomy collapse
qiime taxa collapse \
	--i-table ${data_dir}/2_denoised_table.qza \
	--i-taxonomy ${data_dir}/Greengene/3_taxonomy.qza \
	--p-level ${i} \
	--o-collapsed-table ${level_DIR}/3_1_table_collapse_level_${i}.qza \

#table heatmap
qiime feature-table heatmap \
	--i-table ${level_DIR}/3_1_table_collapse_level_${i}.qza \
	--p-metric braycurtis \
	--p-method centroid \
	--p-cluster both \
	--o-visualization ${level_DIR}/3_2_table_heatmap_collapse_level_${i}.qzv

#barplot
qiime taxa barplot \
  --i-table ${data_dir}/2_denoised_table.qza \
  --i-taxonomy ${data_dir}/Greengene/3_taxonomy.qza \
  --m-metadata-file $metadata \
  --o-visualization ${level_DIR}/3_3_taxa_bar_plots_level_${i}.qzv

done
fi

taxonomy_DIR=${data_dir}/$reference
# taxonomy_DIR=${data_dir}/Greengene

for i in 2 3 4 5 6 7
do
level_dir=$taxonomy_DIR/level-${i}
RESULTS_DIR=$level_dir/metadata
mkdir $RESULTS_DIR

qiime tools export \
--input-path ${level_dir}/3_3_taxa_bar_plots_level_${i}.qzv \
--output-path ${level_dir}

mv ${level_dir}/level-${i}.csv ${RESULTS_DIR}/level-${i}-w-clinic.txt
rm -rf ${level_dir}/level-*.csv ${level_dir}/dist ${level_dir}/q2templateassets ${level_dir}/index.html ${level_dir}/level-*.jsonp

sed -i 's/index/#SampleID/g' $RESULTS_DIR/level-${i}-w-clinic.txt
sed -i 's/,/\t/g' $RESULTS_DIR/level-${i}-w-clinic.txt
done


echo "================================="

echo "Alignment"

echo "================================="

##aligment
qiime alignment mafft \
--i-sequences ${data_dir}/2_denoised_rep-seqs.qza \
--o-alignment ${data_dir}/4_1_aligned_seqs.qza \

#mask (positional consevation and gap filtering)
qiime alignment mask \
--i-alignment ${data_dir}/4_1_aligned_seqs.qza \
--p-max-gap-frequency 1 \
--p-min-conservation 0.4 \
--o-masked-alignment ${data_dir}/4_2_aligned_masked_seq.qza

#phylogeny
qiime phylogeny fasttree \
--i-alignment ${data_dir}/4_2_aligned_masked_seq.qza \
--o-tree ${data_dir}/4_3_tree_unrooted.qza

#rooted tree
qiime phylogeny midpoint-root \
--i-tree ${data_dir}/4_3_tree_unrooted.qza \
--o-rooted-tree ${data_dir}/4_4_tree_rooted.qza

# filter features not in tree
qiime phylogeny filter-table \
--i-table ${data_dir}/2_denoised_table.qza \
--i-tree ${data_dir}/4_4_tree_rooted.qza \
--o-filtered-table ${data_dir}/4_5_filtered_table.qza

#visualization
qiime feature-table summarize \
--i-table ${data_dir}/4_5_filtered_table.qza \
--o-visualization ${data_dir}/4_5_filtered_table.qzv

#export to .biom file
qiime tools export \
--input-path ${data_dir}/4_5_filtered_table.qza \
--output-path ${data_dir}

qiime feature-table relative-frequency \
--i-table ${data_dir}/4_5_filtered_table.qza \
--o-relative-frequency-table ${data_dir}/4_6_relative-frequency

qiime feature-table summarize \
--i-table ${data_dir}/4_6_relative-frequency.qza \
--o-visualization ${data_dir}/4_6_relative-frequency.qzv \
--m-sample-metadata-file $metadata

# Clinic DIR
clinic_DIR="CLINC_DIR"

qiime feature-table filter-samples \
--i-table ${data_dir}/4_5_filtered_table.qza \
--m-metadata-file ${clinic_DIR}/UC_HADS_clinic_metadata.txt \
--o-filtered-table ${data_dir}/4_5_filtered_table_UC.qza

qiime feature-table filter-samples \
--i-table ${data_dir}/4_5_filtered_table.qza \
--m-metadata-file ${clinic_DIR}/CD_HADS_clinic_metadata.txt \
--o-filtered-table ${data_dir}/4_5_filtered_table_CD.qza

qiime feature-table relative-frequency \
--i-table ${data_dir}/4_5_filtered_table_UC.qza \
--o-relative-frequency-table ${data_dir}/4_6_relative-frequency_UC

qiime feature-table relative-frequency \
--i-table ${data_dir}/4_5_filtered_table_CD.qza \
--o-relative-frequency-table ${data_dir}/4_6_relative-frequency_CD

qiime feature-table summarize \
--i-table ${data_dir}/4_6_relative-frequency_UC.qza \
--o-visualization ${data_dir}/4_6_relative-frequency_UC.qzv \
--m-sample-metadata-file ${clinic_DIR}/UC_HADS_clinic_metadata.txt

qiime feature-table summarize \
--i-table ${data_dir}/4_6_relative-frequency_CD.qza \
--o-visualization ${data_dir}/4_6_relative-frequency_CD.qzv \
--m-sample-metadata-file ${clinic_DIR}/CD_HADS_clinic_metadata.txt


echo "================================="

echo "Diversity"

echo "================================="

# Minimum sampling depth
mkdir ${data_dir}/tmp
tmp_dir=${data_dir}/tmp
cp ${data_dir}/4_5_filtered_table.qzv $data_dir/table.zip
unzip -q $data_dir/table.zip -d $tmp_dir
file_path=`realpath $tmp_dir/*`
awk -F, '{print $NF}' $file_path/data/sample-frequency-detail.csv | tail -n1 > $tmp_dir/sample
depth=$(awk -F. '{print $1}' $tmp_dir/sample)
echo -e 'Sampling depth :' $depth
echo -e 'Sampling depth :' $depth > ${data_dir}/Sampling_depth
rm -rf $tmp_dir $data_dir/table.zip

# #
taxonomy_DIR=${data_dir}/$reference
# taxonomy_DIR=${data_dir}/Greengene

for i in 2 3 4 5 6 7
do

RESULTS=$taxonomy_DIR/level-${i}
metadata_DIR=$RESULTS/metadata
#alpha diversity
qiime diversity core-metrics-phylogenetic \
--i-table ${data_dir}/4_5_filtered_table.qza \
--i-phylogeny ${data_dir}/4_4_tree_rooted.qza \
--p-n-jobs-or-threads 2 \
--p-sampling-depth ${depth} \
--m-metadata-file ${metadata} \
--o-rarefied-table ${RESULTS}/5_1_rarefied_table.qza \
--o-faith-pd-vector ${RESULTS}/5_2_alpha_faith_pd.qza \
--o-shannon-vector ${RESULTS}/5_2_alpha_shannon.qza \
--o-evenness-vector ${RESULTS}/5_2_alpha_evenness.qza \
--o-observed-features-vector ${RESULTS}/5_2_alpha_obs_features.qza \
--o-unweighted-unifrac-distance-matrix ${RESULTS}/5_2_beta_un_unifrac.qza \
--o-weighted-unifrac-distance-matrix ${RESULTS}/5_2_beta_weighted_unifrac.qza \
--o-jaccard-distance-matrix ${RESULTS}/5_2_beta_jaccard.qza \
--o-bray-curtis-distance-matrix ${RESULTS}/5_2_beta_braycurtis.qza \
--o-unweighted-unifrac-pcoa-results ${RESULTS}/5_3_pcoa_un_unifrac.qza \
--o-weighted-unifrac-pcoa-results ${RESULTS}/5_3_pcoa_weighted_unifrac.qza \
--o-jaccard-pcoa-results ${RESULTS}/5_3_pcoa_jaccard.qza \
--o-bray-curtis-pcoa-results ${RESULTS}/5_3_pcoa_braycurtes.qza \
--o-unweighted-unifrac-emperor ${RESULTS}/5_4_plot_un_unifrac.qzv \
--o-weighted-unifrac-emperor  ${RESULTS}/5_4_plot_weighted_unifrac.qzv \
--o-jaccard-emperor ${RESULTS}/5_4_plot_jaccard.qzv \
--o-bray-curtis-emperor ${RESULTS}/5_4_plot_braycurtis.qzv

qiime diversity alpha \
--i-table ${data_dir}/4_5_filtered_table.qza \
--p-metric chao1 \
--o-alpha-diversity ${RESULTS}/5_2_alpha_chao1.qza

qiime diversity alpha \
--i-table ${data_dir}/4_5_filtered_table.qza \
--p-metric simpson \
--o-alpha-diversity ${RESULTS}/5_2_alpha_simpson.qza

#qiime diversity alpha-correlation (numeric sample metadata ~ alpha diversity)
qiime diversity alpha-correlation \
--i-alpha-diversity ${RESULTS}/5_2_alpha_chao1.qza \
--m-metadata-file ${metadata} \
--p-method spearman \
--o-visualization ${RESULTS}/5_3_alpha_chao1_correlation.qzv

qiime diversity alpha-correlation \
--i-alpha-diversity ${RESULTS}/5_2_alpha_simpson.qza \
--m-metadata-file ${metadata} \
--p-method spearman \
--o-visualization ${RESULTS}/5_3_alpha_simpson_correlation.qzv

qiime diversity alpha-correlation \
--i-alpha-diversity ${RESULTS}/5_2_alpha_faith_pd.qza \
--m-metadata-file ${metadata} \
--p-method spearman \
--o-visualization ${RESULTS}/5_3_alpha_faith_pd_correlation.qzv

qiime diversity alpha-correlation \
--i-alpha-diversity ${RESULTS}/5_2_alpha_shannon.qza \
--m-metadata-file ${metadata} \
--p-method spearman \
--o-visualization ${RESULTS}/5_3_alpha_shannon_correlation.qzv

qiime diversity alpha-correlation \
--i-alpha-diversity ${RESULTS}/5_2_alpha_evenness.qza \
--m-metadata-file ${metadata} \
--p-method spearman \
--o-visualization ${RESULTS}/5_3_alpha_evenness_correlation.qzv

qiime diversity alpha-correlation \
--i-alpha-diversity ${RESULTS}/5_2_alpha_obs_features.qza \
--m-metadata-file ${metadata} \
--p-method spearman \
--o-visualization ${RESULTS}/5_3_alpha_obs_features_correlation.qzv

# alpha_shannon
qiime tools export \
--input-path ${RESULTS}/5_2_alpha_faith_pd.qza \
--output-path ${RESULTS}/alpha_faith_pd

qiime tools export \
--input-path ${RESULTS}/5_2_alpha_obs_features.qza \
--output-path ${RESULTS}/alpha_obs_features

qiime tools export \
--input-path ${RESULTS}/5_2_alpha_chao1.qza \
--output-path ${RESULTS}/alpha_chao1

qiime tools export \
--input-path ${RESULTS}/5_2_alpha_simpson.qza \
--output-path ${RESULTS}/alpha_simpson

qiime tools export \
--input-path ${RESULTS}/5_2_alpha_shannon.qza \
--output-path ${RESULTS}/alpha_shannon

done

taxonomy_DIR=${data_dir}/$reference

bash alpha_sig_image.sh $taxonomy_DIR
bash alpha_sig_image_group.sh $taxonomy_DIR ControlSevere
bash alpha_sig_image_group_UC+CD.sh $taxonomy_DIR 

echo "================================="

echo "Alpha-rarefaction"

echo "================================="

eval "$(conda shell.bash hook)"

conda activate qiime2

taxonomy_DIR=${data_dir}/$reference

# Minimum sampling depth
qiime tools export \
--input-path ${data_dir}/4_5_filtered_table.qzv \
--output-path ${data_dir}/rarefaction

awk -F, '{print $NF}' ${data_dir}/rarefaction/sample-frequency-detail.csv | tail -n10 | head -n1 > $data_dir/sample_rare
depth=$(awk -F. '{print $1}' $data_dir/sample_rare)
echo -e 'Sampling depth :' $depth
echo -e 'Sampling depth :' $depth > ${data_dir}/Sampling_depth_Rarefaction

# # # # #
for i in 2 3 4 5 6 7
do

RESULTS=$taxonomy_DIR/level-${i}
metadata_DIR=$RESULTS/metadata

qiime diversity alpha-rarefaction \
--i-table ${data_dir}/4_5_filtered_table.qza \
--i-phylogeny ${data_dir}/4_4_tree_rooted.qza \
--p-max-depth ${depth} \
--m-metadata-file ${metadata} \
--o-visualization ${RESULTS}/5_alpha-rarefaction.qzv
done

fi