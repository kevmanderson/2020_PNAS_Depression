function plotParcelData(dscalar_template, write_val_file, save_path)

    addpath('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/matlab');
    addpath('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/matlab/gifti-1.6');

    wb_path = '/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_command'
    cii_dscalar = ciftiopen(dscalar_template, wb_path);

    num_idxs = length(cii_dscalar.cdata)/2;

    rh_hemi = cii_dscalar.cdata((num_idxs+1):length(cii_dscalar.cdata));
    rh_hemi(rh_hemi ~= 0) = rh_hemi(rh_hemi ~= 0)+35;

    cii_dscalar.cdata((num_idxs+1):length(cii_dscalar.cdata)) = rh_hemi;

    write_vals = csvread(write_val_file);
    for row = 1:length(write_vals);
        cii_dscalar.cdata(cii_dscalar.cdata == row) = write_vals(row);
    end
    ciftisave(cii_dscalar, save_path, wb_path)

end
