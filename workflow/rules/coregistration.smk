## Do boundary-based registration to anatomy 

rule fs_bbr:
    input:
        epi = rules.apply_topup.output.firstvol,
        t1w = join(config['fs_dir'],'{subject}/mri/orig.mgz')
    output:
        reg = 'output/fs_bbr/{subject}/epi_rest_bbr.dat',
        mat = 'output/fs_bbr/{subject}/epi_rest_bbr.mat',
        vol = 'output/fs_bbr/{subject}/epi_rest_bbr.nii.gz'
    params:
        fs_setup = config['fs_setup'],      
        optional = '--bold --init-fsl'
    group: 'coregistration'
    log: 'logs/fs_bbr/{subject}.log'  
    singularity: config['fmriprep']        
    shell:
        "{params.fs_setup} bbregister --s {wildcards.subject} --mov {input.epi} --reg {output.reg} --o {output.vol} {params.optional} && "
        "tkregister2 --noedit --reg {output.reg} --mov {input.epi} --targ {input.t1w} --fslregout {output.mat} &> {log}"      

rule inverse_bbr:
    input:
        epi = rules.apply_topup.output.gdc,
        t1w = join(config['fs_dir'],'{subject}/mri/orig.nii.gz'),
        mat = rules.fs_bbr.output.mat
    output:
        inverse_mat = 'output/fs_bbr/{subject}/epi_rest_bbr_inverse.mat',
        t1w_coreg = 'output/fs_bbr/{subject}/t1w_bbr_inverse.nii.gz'
    group: 'coregistration'
    shell:
        "convert_xfm -omat {output.inverse_mat} -inverse {input.mat} && "
        "flirt -interp spline -in {input.t1w} -ref {input.epi} -applyxfm -init {output.inverse_mat} -out {output.t1w_coreg}"

rule apply_bbr:
    input:
        epi = rules.apply_topup.output.gdc,
        bbr = rules.fs_bbr.output.reg,
        targ = join(config['fs_dir'],'{subject}/mri/orig.mgz')
    output: 'output/apply_bbr/{subject}/epi_rest_gdc_bbr.nii.gz'
    group: 'coregistration'
    params:
        fs_setup = config['fs_setup']
    container: config['fmriprep']
    shell:
        "{params.fs_setup} mri_vol2vol --mov {input.epi} --targ {input.targ} --o {output} --reg {input.bbr} --no-resample"

## Do volume-based registration to MNI space 

rule flirt:
    input: join(config['fs_dir'],'{subject}/mri/orig.nii.gz')
    output: 'output/coreg_anat/{subject}/linear/Anat2MNILinear.mat'
    params:
        prefix = 'output/coreg_anat/{subject}/linear/Anat2MNILinear',
        MNI = config['MNI'],
        optional = '-interp spline -dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -v'
    group: 'coregistration'
    log: 'logs/flirt/{subject}.log'    
    shell:
        "flirt {params.optional} -in {input} -ref {params.MNI} -omat {output} -out {params.prefix} &> {log}"

rule fnirt:
    input:
        t1w = join(config['fs_dir'],'{subject}/mri/orig.nii.gz'),
        affine = rules.flirt.output
    output:
        fout = 'output/coreg_anat/{subject}/nonlinear/NonlinearRegWarp.nii.gz',
        jout = 'output/coreg_anat/{subject}/nonlinear/NonlinearRegJacobians.nii.gz',
        refout = 'output/coreg_anat/{subject}/nonlinear/IntensityModulatedT1.nii.gz',
        iout = 'output/coreg_anat/{subject}/nonlinear/NonlinearWarped.nii.gz',
        logout = 'output/coreg_anat/{subject}/nonlinear/NonlinearReg.txt',
        intout = 'output/coreg_anat/{subject}/nonlinear/NonlinearIntensities.nii.gz',
        cout = 'output/coreg_anat/{subject}/nonlinear/NonlinearReg.nii.gz'
    params:
        MNI = config['MNI']
    group: 'coregistration'
    log: 'logs/fnirt/{subject}.log'
    threads: 8
    resources:
        time = 120,
        mem_mb = 32000
    shell:
        "fnirt --in={input.t1w} --ref={params.MNI} --aff={input.affine} --fout={output.fout} --jout={output.jout} "
        "--refout={output.refout} --iout={output.iout} --logout={output.logout} --intout={output.intout} --cout={output.cout} -v &> {log}"

rule inverse_fnirt:
    input:
        t1w = join(config['fs_dir'],'{subject}/mri/orig.nii.gz'),
        warp = rules.fnirt.output.fout
    output: 'output/coreg_anat/{subject}/nonlinear/NonlinearRegInverseWarp.nii.gz'
    group: 'coregistration'
    shell:
        "invwarp -w {input.warp} -o {output} -r {input.t1w} --noconstraint"

rule combine_warp:
    input:
        topup = rules.run_topup.output,
        bbr = rules.fs_bbr.output.mat,
        fnirt = rules.fnirt.output.cout
    output:
        scaled = 'output/combine_warp/{subject}/epi_fout_scaled.nii.gz',
        combined = 'output/combine_warp/{subject}/epi_rest_gdc_bbr_mni.nii.gz'
    params:
        MNI = config['MNI']
    group: 'coregistration' 
    shell:
        "fslmaths {input.topup}/epi_fout.nii.gz -mul 0.0343097 {output.scaled} && "
        "convertwarp --relout --rel -s {output.scaled} --premat={input.bbr} --warp1={input.fnirt} --ref={params.MNI} --out={output.combined}"

## Register EPI data to MNI space 

rule apply_combine_warp:    
    input: 
        warp = rules.combine_warp.output.combined,
        epi = rules.apply_topup.output.firstvol,
        ref = config['MNI']
    output: 'output/final/{subject}/epi_rest_final.nii.gz'
    group: 'coregistration' 
    shell:
        "applywarp -i {input.epi} -r {input.ref} -o {output} -w {input.warp} --rel --interp=spline"
