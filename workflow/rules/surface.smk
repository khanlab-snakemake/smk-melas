## Generate midthickness and inflated surfaces

rule convert_surfaces:
    input: join(config['fs_dir'],'{subject}/surf/{hemi}.{surf}')
    output: 'output/midthickness_surface/{subject}/{hemi}.{surf}.surf.gii'
    params:
        fs_setup = config['fs_setup']
    container: config['fmriprep']
    shell:
        "{params.fs_setup} mris_convert {input} {output}"

rule generate_midthickness:
    input:
        white = 'output/midthickness_surface/{subject}/{hemi}.white.surf.gii',
        pial = 'output/midthickness_surface/{subject}/{hemi}.pial.surf.gii'
    output: 'output/midthickness_surface/{subject}/{hemi}.midthickness.surf.gii'
    container: config['connectome_workbench']      
    shell:
        "wb_command -surface-cortex-layer {input.white} {input.pial} 0.5 {output}"      

rule get_tkr2scanner:
    input: join(config['fs_dir'],'{subject}/mri/orig.mgz')
    output: 'output/tkr2scanner/{subject}/tkr2scanner.xfm'
    params:
        fs_setup = config['fs_setup']
    container: config['fmriprep']
    shell: 
        "{params.fs_setup} mri_info {input} --tkr2scanner > {output}"

rule apply_surf_tkr2scanner:
    input: 
        surf = rules.generate_midthickness.output,
        tkr2scanner = rules.get_tkr2scanner.output
    output: 'output/midthickness_surface_anat/{subject}/{hemi}.midthickness.native.surf.gii'
    threads: 8
    container: config['connectome_workbench']
    shell: 
        "wb_command -surface-apply-affine {input.surf} {input.tkr2scanner} {output}"    

rule inflate_surface_native:
    input: rules.apply_surf_tkr2scanner.output
    output:
        inflated = 'output/midthickness_surface_anat/{subject}/{hemi}.inflated.native.surf.gii',
        very_inflated = 'output/midthickness_surface_anat/{subject}/{hemi}.very_inflated.native.surf.gii',
    params: "-iterations-scale 2.5"
    threads: 8
    container: config['connectome_workbench']
    shell:     
        "wb_command -surface-generate-inflated {input} {output.inflated} {output.very_inflated} {params}"   

## Bring midthickness surface in EPI space and sample data

rule get_scanner2epi:
    input: 
        bbr = rules.inverse_bbr.output.inverse_mat,
        t1w = join(config['fs_dir'],'{subject}/mri/orig.nii.gz'),
        epi = rules.apply_topup.output.firstvol_gdc
    output: 'output/surf2epi/{subject}/surf2epi.xfm'
    threads: 8
    container: config['connectome_workbench']    
    shell:
        "wb_command -convert-affine -from-flirt {input.bbr} {input.t1w} {input.epi} -to-world {output}"

rule apply_surf_scanner2epi:
    input:
        surf = rules.apply_surf_tkr2scanner.output,
        scanner2epi = rules.get_scanner2epi.output
    output: 'output/midthickness_surface_epi/{subject}/{hemi}.midthickness.epi.surf.gii'
    threads: 8
    container: config['connectome_workbench']
    shell:     
        "wb_command -surface-apply-affine {input.surf} {input.scanner2epi} {output}"

rule inflate_surface_epi:
    input: rules.apply_surf_scanner2epi.output
    output:
        inflated = 'output/midthickness_surface_epi/{subject}/{hemi}.inflated.epi.surf.gii',
        very_inflated = 'output/midthickness_surface_epi/{subject}/{hemi}.very_inflated.epi.surf.gii',
    params: "-iterations-scale 2.5"
    threads: 8
    container: config['connectome_workbench']
    shell:     
        "wb_command -surface-generate-inflated {input} {output.inflated} {output.very_inflated} {params}"   

rule volume_to_surface:
    input: 
        epi = rules.apply_topup.output.gdc,
        surf = rules.apply_surf_scanner2epi.output
    output: 'output/epi_surf/{subject}/{hemi}.rsfMRI.func.gii'
    threads: 8
    container: config['connectome_workbench']    
    shell: 
        "wb_command -volume-to-surface-mapping {input.epi} {input.surf} {output} -trilinear"
        
## Resample surface data to 32k_fs_LR space

rule resample_surface:
    input: 
        surf = rules.apply_surf_tkr2scanner.output,
        sphere_old = 'output/midthickness_surface/{subject}/{hemi}.sphere.reg.surf.gii',
        sphere_new = 'resources/fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii'
    output: 'output/midthickness_surface_fsLR/{subject}/{hemi}.midthickness.32k_fs_LR.surf.gii'
    container: config['connectome_workbench']      
    shell:    
        "wb_command -surface-resample {input.surf} {input.sphere_old} {input.sphere_new} BARYCENTRIC {output}"

rule inflate_surface_fsLR:
    input: rules.apply_surf_scanner2epi.output
    output:
        inflated = 'output/midthickness_surface_fsLR/{subject}/{hemi}.inflated.32k_fs_LR.surf.gii',
        very_inflated = 'output/midthickness_surface_fsLR/{subject}/{hemi}.very_inflated.32k_fs_LR.surf.gii',
    params: "-iterations-scale 1.0"
    threads: 8
    container: config['connectome_workbench']
    shell:     
        "wb_command -surface-generate-inflated {input} {output.inflated} {output.very_inflated} {params}"  

rule resample_surface_data:
    input:
        metric = rules.volume_to_surface.output,
        sphere_old = 'output/midthickness_surface/{subject}/{hemi}.sphere.reg.surf.gii',
        sphere_new = 'resources/fs_LR-deformed_to-fsaverage.{hemi}.sphere.32k_fs_LR.surf.gii',
        area_old = rules.generate_midthickness.output,
        area_new = rules.resample_surface.output
    output: 'output/epi_surf_fsLR/{subject}/{hemi}.rsfMRI.32k_fs_LR.func.gii'
    threads: 8
    container: config['connectome_workbench']    
    shell:     
        "wb_command -metric-resample {input.metric} {input.sphere_old} {input.sphere_new} ADAP_BARY_AREA {output} -area-surfs {input.area_old} {input.area_new}"
        
## Warp midthickness surface to MNI space

rule warp_surface:
    input:
        surf = rules.resample_surface.output,
        inverse_warp = rules.inverse_fnirt.output,
        warp = rules.fnirt.output.fout
    output: 'output/midthickness_surface_mni/{subject}/{hemi}.midthickness.32k_fs_LR.surf.gii'
    threads: 8
    container: config['connectome_workbench']   
    shell:
        "wb_command -surface-apply-warpfield {input.surf} {input.inverse_warp} {output} -fnirt {input.warp}"

rule inflate_surface_mni:
    input: rules.warp_surface.output
    output:
        inflated = 'output/midthickness_surface_mni/{subject}/{hemi}.inflated.mni.32k_fs_LR.surf.gii',
        very_inflated = 'output/midthickness_surface_mni/{subject}/{hemi}.very_inflated.mni.32k_fs_LR.surf.gii',
    params: "-iterations-scale 1.0"
    threads: 8
    container: config['connectome_workbench']
    shell:     
        "wb_command -surface-generate-inflated {input} {output.inflated} {output.very_inflated} {params}"  