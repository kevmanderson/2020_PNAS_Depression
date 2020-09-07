function plotParcelData(dscalar_template, parcel_info_file, write_val_file, save_path)

    addpath('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/matlab');
    addpath('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/matlab/gifti-1.6');

    wb_path = '/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_command'
    cii_dscalar = ciftiopen(dscalar_template, wb_path);


    fileID    = fopen(parcel_info_file, 'r');
    parcel_in = textscan(fileID, '%s', 'Delimiter', '\n');
    net_idxs  = cellfun(@isempty, strfind(parcel_in{1}, 'Network')) == 0;
    parcel_list = parcel_in{1}(net_idxs);

    write_vals = csvread(write_val_file);
    for row = 1:length(write_vals);
        cii_dscalar.cdata(cii_dscalar.cdata == row) = write_vals(row);
    end
    ciftisave(cii_dscalar, save_path, wb_path)

end

