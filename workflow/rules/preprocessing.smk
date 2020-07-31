rule slicetimer:
    input: join(config['input_dir'],'{subject}/rsfmri/{scan}.nii.gz')
    output: 'output/slicetimer/{subject}/st_{scan}.nii.gz'
    params:
        slicetiming = config['slicetiming']
    group: 'preprocessing'     
    singularity: config['fmriprep']
    log: 'logs/slicetimer/{subject}_{scan}.log'
    shell:
        "TR=`fslval {input[0]} pixdim4` && "
        "3dTshift -prefix {output} -tpattern @{params.slicetiming} -TR $TR -tzero 0.0 {input} &> {log}"

rule gradient_unwarp:
    input: 'output/slicetimer/{subject}/st_rest.nii.gz'
    output:
        firstvol = 'output/gradcorrect/{subject}/epi_rest_firstvol.nii.gz',
        warp = 'output/gradcorrect/{subject}/epi_gdc_warp.nii.gz',
    params:
        coeff = config['coeff'],
        scratch = directory('output/gradcorrect/{subject}/scratch')
    group: 'preprocessing'
    log: 'logs/gradient_unwarp/{subject}.log'             
    singularity: config['gradcorrect']
    shell:
        "fslroi {input} {output.firstvol} 0 1 && "
        "procGradCorrect -i {output.firstvol} -g {params.coeff} -w {output.warp} -s {params.scratch} &> {log}"

rule apply_gradient_unwarp:
    input:
        epi = rules.slicetimer.output,
        warp = rules.gradient_unwarp.output.warp
    output:
        unwarped = 'output/gradcorrect/{subject}/epi_{scan}_unwarped.nii.gz',
        jacobian = 'output/gradcorrect/{subject}/epi_{scan}_jacobian.nii.gz',
        corrected = 'output/gradcorrect/{subject}/epi_{scan}_intcor.nii.gz'
    group: 'preprocessing'          
    singularity: config['gradcorrect']         
    shell:
        "applywarp -i {input.epi} -o {output.unwarped} -r {input.epi} -w {input.warp} --abs --interp=spline && "
        "fslroi {output.unwarped} {output.jacobian} 0 1 && "
        "reg_jacobian -ref {output.jacobian} -def {input.warp} -jac {output.jacobian} && "
        "fslmaths {output.jacobian} -mul -1 -abs {output.jacobian} && "
        "fslmaths {output.unwarped} -mul {output.jacobian} {output.corrected} && "
        "fslcpgeom {input.epi} {output.corrected}"

rule prep_topup:
    input: expand('output/gradcorrect/{{subject}}/epi_{scan}_intcor.nii.gz', scan=config['scans'])
    output:
        first_vols = 'output/topup/{subject}/prep_topup/epi_rest_firstvols.nii.gz',
        concat_out = 'output/topup/{subject}/prep_topup/epi_rest-topup_concat.nii.gz'
    group: 'preprocessing'     
    shell:
        "fslroi {input[0]} {output.first_vols} 0 5 && "
        "fslmerge -tr {output.concat_out} {output.first_vols} {input[1]} 2"

rule run_topup:
    input: rules.prep_topup.output.concat_out
    output: directory('output/topup/{subject}/run_topup')
    params: 
        topup = config['topup_config'],
        acquisition = config['acqparams']
    group: 'preprocessing'        
    log: 'logs/run_topup/{subject}.log'
    threads: 8
    resources:
        time = 360,
        mem_mb = 32000
    shell:
        "mkdir -p {output} && "
        "topup --imain={input} --datain={params.acquisition} --config={params.topup} "
        "--out={output}/epi --fout={output}/epi_fout --iout={output}/epi_iout -v &> {log}"

rule mc_tseries:
    input: 'output/gradcorrect/{subject}/epi_rest_intcor.nii.gz'
    output: 'output/motioncor/{subject}/epi_rest_mc.nii.gz'
    params:
        optional = '-refvol 0 -sinc_final -plots -mats -verbose 2'
    group: 'preprocessing'
    log: 'logs/motioncor/{subject}.log' 
    shell:
        "mcflirt -in {input} -out {output} {params.optional} &> {log}"

rule apply_topup:
    input:
        epi = rules.mc_tseries.output,
        topup = rules.run_topup.output
    output:
        gdc = 'output/topup/{subject}/apply_topup/epi_rest_gdc.nii.gz',
        firstvol = 'output/topup/{subject}/apply_topup/epi_rest_firstvol.nii.gz',
        firstvol_gdc = 'output/topup/{subject}/apply_topup/epi_rest_firstvol_gdc.nii.gz'
    params:
        acquisition = config['acqparams']
    group: 'preprocessing'
    log: 'logs/apply_topup/{subject}.log'       
    shell:
        "fslroi {input.epi} {output.firstvol} 0 1 && "
        "applytopup --imain={output.firstvol} --inindex=1 --datain={params.acquisition} --topup={input.topup}/epi --out={output.firstvol_gdc} --method=jac -v &> {log} && "
        "applytopup --imain={input.epi} --inindex=1 --datain={params.acquisition} --topup={input.topup}/epi --out={output.gdc} --method=jac -v &> {log}"
