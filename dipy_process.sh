#!/bin/sh

dir1=/media/juan2634/ROCKET-NANO/MRI_Analysis/ELP/BIDS_ELP
dir2=/media/juan2634/ROCKET-NANO/MRI_Analysis/ELP/BIDS_ELP/derivatives/DWI/DIPY/ELA-ELP
dir3=/media/juan2634/ROCKET-NANO/MRI_Analysis/ELP/BIDS_ELP/derivatives/DWI/DIPY/ELA-ELP/TRACT

for x in `cat ids`
do
    echo ${x}
    mkdir ${dir2}/tmp
    mkdir ${dir2}/ELP
    mkdir ${dir2}/tmp/${x}
    mkdir ${dir2}/tmp/${x}/vol
    mkdir ${dir2}/tmp/${x}/b0_correction_input
    mkdir ${dir2}/tmp/${x}/b0_correction_output
    mkdir ${dir2}/tmp/${x}/parameters_dwi
    mkdir ${dir2}/tmp/${x}/tractography
    
    mkdir ${dir3}/subjects
    mkdir ${dir3}/subjects/${x}
    mkdir ${dir3}/subjects/${x}/anatomical_measures
    mkdir ${dir3}/subjects/${x}/rec_bundles
    mkdir ${dir3}/subjects/${x}/org_bundles
    
    cp ${dir1}/${x}/anat/${x}_run-1_T1w.nii.gz ${dir2}/tmp/${x}/${x}_run-1_T1w.nii.gz
    cp ${dir1}/${x}/dwi/${x}_run-1_dwi.nii.gz ${dir2}/tmp/${x}/${x}_run-1_dwi.nii.gz
    cp ${dir1}/${x}/dwi/${x}_run-1_dwi.bval ${dir2}/tmp/${x}/${x}_run-1_dwi.bval
    cp ${dir1}/${x}/dwi/${x}_run-1_dwi.bvec ${dir2}/tmp/${x}/${x}_run-1_dwi.bvec

    dipy_denoise_nlmeans \
        ${dir2}/tmp/${x}/${x}_run-1_T1w.nii.gz \
        --sigma 2 \
        --patch_radius 2 \
        --out_denoised \
        ${dir2}/tmp/${x}/${x}_run-1_T1w_denoise.nii.gz

    # dipy_denoise_lpca \
    #    ${dir2}/tmp/${x}/${x}_run-1_dwi.nii.gz \
    #    ${dir2}/tmp/${x}/${x}_run-1_dwi.bval \
    #    ${dir2}/tmp/${x}/${x}_run-1_dwi.bvec \
    #    --patch_radius 3 \
    #    --out_denoised \
    #    ${dir2}/tmp/${x}/${x}_run-1_dwi_denoise.nii.gz

    dipy_gibbs_ringing \
        ${dir2}/tmp/${x}/${x}_run-1_dwi_denoise.nii.gz \
        --num_processes -1 \
        --out_unring ${dir2}/tmp/${x}/${x}_run-1_dwi_denoise_unr.nii.gz
        
    cp ${dir2}/tmp/${x}/${x}_run-1_dwi_denoise_unr.nii.gz ${dir2}/tmp/${x}/vol/
    cd ${dir2}/tmp/${x}/vol/

    fslsplit ${dir2}/tmp/${x}/${x}_run-1_dwi_denoise_unr.nii.gz -t
    cd ${dir2}/tmp/${x}/
    cp ${dir2}/tmp/${x}/vol/vol0000.nii.gz ${dir2}/tmp/${x}/vol0000.nii.gz
    rm -fr ${dir2}/tmp/${x}/vol
    cp /media/juan2634/ROCKET-NANO/MRI_Analysis/ELP/BIDS_ELP/derivatives/Morfometria/DARTEL/CTR-ELP/${x}/c2${x}.nii ${dir2}/tmp/${x}/c2${x}_run-1_T1w.nii

    dipy_align_affine \
        ${dir2}/tmp/${x}/vol0000.nii.gz \
        ${dir2}/tmp/${x}/c2${x}_run-1_T1w.nii\
        --transform "affine" \
        --out_moved ${dir2}/tmp/${x}/c2${x}_c2_to_dwi.nii

    dipy_median_otsu ${dir2}/tmp/${x}/${x}_run-1_dwi_denoise_unr.nii.gz  --out_mask ${dir2}/tmp/${x}/${x}_dwi_mask.nii.gz --vol_idx 3 4 5 6

    #Procesado para corregir las distoriones de eddy mendiante machine learning usando el contenedor docker
    
    	#Preparamos las carpetas para el procesado
    cp ${dir2}/tmp/${x}/vol0000.nii.gz ${dir2}/tmp/${x}/b0_correction_input/b0.nii.gz
    #cp ${dir2}/tmp/${x}/${x}_run-1_T1w_denoise.nii.gz ${dir2}/tmp/${x}/b0_correction_input/T1.nii.gz
    cp ${dir2}/acqparams.txt ${dir2}/tmp/${x}/b0_correction_input/acqparams.txt
    
    docker run --rm -v ${dir2}/tmp/${x}/b0_correction_input/:/INPUTS/ -v ${dir2}/tmp/${x}/b0_correction_output/:/OUTPUTS/ \
    	-v $FREESURFER_HOME/license.txt:/extra/freesurfer/license.txt leonyichencai/synb0-disco --stripped
    	
    eddy --imain=${dir2}/tmp/${x}/${x}_run-1_dwi_denoise_unr.nii.gz --mask=${dir2}/tmp/${x}/${x}_dwi_mask.nii.gz \
    	--acqp=${dir2}/tmp/${x}/b0_correction_input/acqparams.txt --index=${dir2}/index.txt \
    	--bvecs=${dir2}/tmp/${x}/${x}_run-1_dwi.bvec --bvals=${dir2}/tmp/${x}/${x}_run-1_dwi.bval \
    	--topup=${dir2}/tmp/${x}/b0_correction_output/topup --out=${dir2}/tmp/${x}/eddy_unwarped_images
    	
    #cambio de nombre
    cp ${dir2}/tmp/${x}/eddy_unwarped_images.nii.gz ${dir2}/tmp/${x}/${x}_run-1_dwi_denoise_unr_ecorr.nii.gz

    dipy_fit_csd \
        ${dir2}/tmp/${x}/${x}_run-1_dwi_denoise_unr_ecorr.nii.gz \
        ${dir2}/tmp/${x}/${x}_run-1_dwi.bval \
        ${dir2}/tmp/${x}/eddy_unwarped_images.eddy_rotated_bvecs \
        ${dir2}/tmp/${x}/${x}_dwi_mask.nii.gz \
        --fa_thr 0.7 \
        --sh_order 8 \
        --parallel \
        --out_dir ${dir2}/tmp/${x}/parameters_dwi \
        --parallel \
        --num_processes -1 \
        --out_pam csd_peaks.pam5

    dipy_fit_dti \
        ${dir2}/tmp/${x}/${x}_run-1_dwi_denoise_unr_ecorr.nii.gz  \
        ${dir2}/tmp/${x}/${x}_run-1_dwi.bval \
        ${dir2}/tmp/${x}/eddy_unwarped_images.eddy_rotated_bvecs \
        ${dir2}/tmp/${x}/${x}_dwi_mask.nii.gz \
        --save_metrics "md" "fa" "ad" "rd" \
        --out_dir ${dir2}/tmp/${x}/parameters_dwi

    dipy_fit_dki \
        ${dir2}/tmp/${x}/${x}_run-1_dwi_denoise_unr_ecorr.nii.gz  \
        ${dir2}/tmp/${x}/${x}_run-1_dwi.bval \
        ${dir2}/tmp/${x}/eddy_unwarped_images.eddy_rotated_bvecs \
        ${dir2}/tmp/${x}/${x}_dwi_mask.nii.gz \
        --b0_threshold 70.0 \
        --out_dir ${dir2}/tmp/${x}/parameters_dwi \
        --out_mk mk.nii.gz \
        --out_ak ak.nii.gz \
        --out_rk rk.nii.gz \
        --out_fa fa_dwi.nii.gz \
        --out_ad ad_dwi.nii.gz \
        --out_rd rk_dwi.nii.gz \
        --out_md mk_dwi.nii.gz \
        --force

    dipy_fit_csa \
       ${dir2}/tmp/${x}/${x}_run-1_dwi_denoise_unr_ecorr.nii.gz  \
       ${dir2}/tmp/${x}/${x}_run-1_dwi.bval \
       ${dir2}/tmp/${x}/eddy_unwarped_images.eddy_rotated_bvecs \
       ${dir2}/tmp/${x}/${x}_dwi_mask.nii.gz \
       --extract_pam_values \
       --out_dir ${dir2}/tmp/${x}/parameters_dwi  \
       --out_pam csa_peaks.pam5
        
    cp ${dir2}/tmp/${x}/parameters_dwi/ad.nii.gz ${dir3}/subjects/${x}/anatomical_measures/ad.nii.gz
    cp ${dir2}/tmp/${x}/parameters_dwi/fa.nii.gz ${dir3}/subjects/${x}/anatomical_measures/fa.nii.gz
    cp ${dir2}/tmp/${x}/parameters_dwi/md.nii.gz ${dir3}/subjects/${x}/anatomical_measures/md.nii.gz
    cp ${dir2}/tmp/${x}/parameters_dwi/rd.nii.gz ${dir3}/subjects/${x}/anatomical_measures/rd.nii.gz
    cp ${dir2}/tmp/${x}/parameters_dwi/csd_peaks.pam5 ${dir3}/subjects/${x}/anatomical_measures/csd_peaks.pam5
    
    dipy_mask \
        ${dir2}/tmp/${x}/parameters_dwi/gfa.nii.gz \
        0.25 \
        --out_dir ${dir2}/tmp/${x}/tractography

    dipy_track \
        ${dir3}/subjects/${x}/anatomical_measures/csd_peaks.pam5 \
        ${dir2}/tmp/${x}/tractography/mask.nii.gz \
        ${dir2}/tmp/${x}/c2${x}_c2_to_dwi.nii \
        --seed_density 2 \
        --tracking_method "prob" \
        --out_dir ${dir2}/tmp/${x}/tractography\
        --out_tractogram streamlines.trk

    dipy_slr \
        ${dir2}/Atlas_30_Bundles/whole_brain/whole_brain_MNI.trk \
        ${dir2}/tmp/${x}/tractography/streamlines.trk \
        --force \
        --out_dir ${dir2}/tmp/${x}/tractography

    for y in `cat ${dir2}/bundles_filesname`
    do
        echo ${y}
        dipy_recobundles \
            ${dir2}/tmp/${x}/tractography/moved.trk \
            ${dir2}/Atlas_30_Bundles/bundles/${y}.trk \
            --force \
            --mix_names \
            --out_dir ${dir3}/subjects/${x}/rec_bundles

        dipy_labelsbundles \
            ${dir2}/tmp/${x}/tractography/streamlines.trk \
            ${dir3}/subjects/${x}/rec_bundles/moved_${y}__labels.npy \
            --mix_names \
            --out_dir ${dir3}/subjects/${x}/org_bundles
    done
done 
