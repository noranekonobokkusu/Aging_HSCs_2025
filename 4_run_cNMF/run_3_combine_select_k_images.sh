#!/bin/bash
#$ -l h_vmem=8G
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=2:00:00
#$ -R y
#$ -cwd

base_dir=/broad/vangalenlab/safina/projects/aging_hsc/prepare_for_publication/4_run_cNMF

for s in `cat ${base_dir}/sample_list.txt`
do
  cd ${base_dir}/cNMF_out/${s}/sample

  # Create a placeholder image for missing components
  convert -size 600x600 xc:white placeholder.png
  convert -size 1500x1200 xc:white heatmap_placeholder.png

  # Find all unique values of k
  ks=$(ls sample.clustering.k_*.png | sed 's/.*.k_//g' | sed 's/.dt.*//g' | sort -n | uniq)

  # Prepare a directory for temporary images if needed
  rm -r temp_images
  mkdir -p temp_images

  # Loop over each unique value of k
  for k in {3..10}; do
    for dt in 0_05 0_1 0_15; do
      if ! [ -f "sample.clustering.k_${k}.dt_${dt}.png" ]; then
        cp placeholder.png "temp_images/sample.clustering.k_${k}.dt_${dt}.png"
        cp heatmap_placeholder.png temp_images/k_${k}.dt_${dt}.heatmap.png

      else
        cp "sample.clustering.k_${k}.dt_${dt}.png" temp_images/
        cp postplots/k_${k}.dt_${dt}.heatmap.png temp_images/
      fi
    done
  done

  # Run montage on the organized files in the temp_images directory
  montage -label "$s - %t" temp_images/sample.clustering.k_*.png ${base_dir}/cNMF_out/${s}/sample/sample.k_selection.png -tile 3x9 -geometry 600x+2+2 -pointsize 100 -pointsize 20 ${base_dir}/cNMF_out/result.${s}.sample.png

  montage -label "$s - %t" temp_images/*.heatmap.png ${base_dir}/cNMF_out/${s}/sample/sample.k_selection.png -tile 3x9 -geometry 600x+2+2 -pointsize 100 -pointsize 20 ${base_dir}/cNMF_out/result.${s}.heatmap.png

  montage ${base_dir}/cNMF_out/result.${s}.*.png -tile 2x1 -geometry x3000+2+2 -pointsize 100 -pointsize 20 ${base_dir}/cNMF_out/${s}.final.png

  # Clean up
  rm ${base_dir}/cNMF_out/result.${s}.*png
  rm -r temp_images
  rm *placeholder.png
done

cd ${base_dir}/cNMF_out/
montage -label "%t" *.final.png -tile 22x1 -geometry +5+5 -pointsize 20 total_result.sample.png
convert total_result.sample.png -resize 50% total_result.sample.resized.png
rm *.final.png
# rm total_result.sample.png
