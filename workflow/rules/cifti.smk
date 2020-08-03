## Combine volume and surface data into CIFTI .dtseries.nii file

rule transform_rois:
    input:
        parc = join(config['fs_dir'],'{subject}/mri/aparc+aseg.mgz'),
        warp = rules.fnirt.output.fout
    output:
        parc_anat = 'output/parc_anat/{subject}/aparc+aseg.nii.gz',
        parc_mni = 'output/parc_mni/{subject}/aparc+aseg.mni.nii.gz'
    params:
        mni = config['MNI'],
        fs_setup = config['fs_setup']    
    container: config['fmriprep']        
    shell:
        "{params.fs_setup} mri_convert {input.parc} {output.parc_anat} && "
        "applywarp -i {output.parc_anat} -r {params.mni} -w {input.warp} -o {output.parc_mni} -d int --interp=nn"

rule extract_rois:
    input: rules.transform_rois.output.parc_mni
    output:
        wm = 'output/rois_mni/{subject}/wm.nii.gz',
        csf = 'output/rois_mni/{subject}/csf.nii.gz',
        rois = 'output/rois_mni/{subject}/rois.nii.gz'
    params:
        fs_setup = config['fs_setup']         
    container: config['fmriprep']
    shell:
        "{params.fs_setup} mri_binarize --i {input} --match 2 41 --erode 2 --o {output.wm} && "
        "mri_binarize --i {input} --match 4 43 --erode 1 --o {output.csf} && "
        "mri_binarize --i {input} --match 8 10 11 12 13 16 17 18 26 28 47 49 50 51 52 53 54 58 60 --o {output.rois} && "
        "fslmaths {output.rois} -mul {input} {output.rois}"

rule generate_gii_label:
    input:
        rois = rules.extract_rois.output.rois,
        labels = 'resources/Atlas_ROIs.1p6mm.txt'
    output: 'output/rois_labels/{subject}/rois.nii.gz'
    container: config['connectome_workbench']
    shell:
        "wb_command -volume-label-import {input.rois} {input.labels} {output}"

rule create_dtseries:
    input:
        lh = 'output/epi_surf_fsLR/{subject}/lh.rsfMRI.32k_fs_LR.func.gii',
        rh = 'output/epi_surf_fsLR/{subject}/rh.rsfMRI.32k_fs_LR.func.gii',
        vol = rules.combine_warps_and_apply.output.warped_mni,
        rois = rules.generate_gii_label.output
    output: 'output/cifti/{subject}/rsfMRI.32k_fs_LR.dtseries.nii'
    container: config['connectome_workbench'] 
    threads: 8
    resources:
        mem_mb = 32000
    shell:
        "wb_command -cifti-create-dense-timeseries {output} -volume {input.vol} {input.rois} -left-metric {input.lh} -right-metric {input.rh} -timestep 2.0"
        