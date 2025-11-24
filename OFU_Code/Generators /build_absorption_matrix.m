function A = build_absorption_matrix(min_w, max_w, species_bool, num_points, data_folder)
% BUILD_ABSORPTION_MATRIX
%   Builds an absorption matrix A by loading real absorption spectra
%   from a user-specified folder, cropping each spectrum to the desired
%   wavelength range, and interpolating them to a fixed number of points.
%
% INPUTS:
%   min_w, max_w   : wavelength range to crop
%   species_bool   : logical vector selecting species to include
%                    [Hb, HbO, Lipid, ICG, MethyleneBlue]
%   num_points     : number of interpolated wavelength samples
%   data_folder    : directory containing the spectral CSV files
%
% OUTPUT:
%   A              : num_species × num_points matrix
%
% EXAMPLE:
%   folder = '/Users/calvin/Data/Spectra/';
%   A = build_absorption_matrix(700, 900, [1 1 1 1 0], 200, folder);

    % Ensure the folder exists
    if ~isfolder(data_folder)
        error('Data folder does not exist: %s', data_folder);
    end

    %% -----------------------------------------------------------
    % Build file paths using fullfile (OS-safe)
    % -----------------------------------------------------------
    path_Hb      = fullfile(data_folder, 'HbR_Spectrum.csv');
    path_HbO      = fullfile(data_folder, 'Hb0_Spectrum.csv');
    path_lipid      = fullfile(data_folder, 'Lipid_Spectrum.csv');
    path_ICG        = fullfile(data_folder, 'ICG_Spectrum.csv');
    path_methylene  = fullfile(data_folder, 'MethyleneBlue_spectrum.csv');


    %% -----------------------------------------------------------
    % Load Hb + HbO spectra
    % -----------------------------------------------------------
    load_Hb = readmatrix(path_Hb);
    Hb_data     = load_Hb(3:end, [1,2]);   % wavelength + Hb
  

    load_HbO = readmatrix(path_HbO);
    HbO_data     = load_HbO(3:end, [1,2]);   % wavelength + Hb

    %% -----------------------------------------------------------
    % Load lipid spectrum
    % -----------------------------------------------------------
    load_lipid = readmatrix(path_lipid);
    lipid_data = load_lipid;
    lipid_data(:,2) = lipid_data(:,2) * 1e3;   % scale to match others


    %% -----------------------------------------------------------
    % Load ICG spectrum
    % -----------------------------------------------------------
    load_ICG = readmatrix(path_ICG);
    ICG_data = load_ICG(3:end, [1,2]);
    ICG_data(:,2) = 0.05 * ICG_data(:,2);      % amplitude normalization


    %% -----------------------------------------------------------
    % Load Methylene Blue
    % -----------------------------------------------------------
    load_MB = readmatrix(path_methylene);
    methyl_blue = load_MB(3:end, :);
    methyl_blue(:,2) = 0.1 * methyl_blue(:,2);


    %% -----------------------------------------------------------
    % Package species in consistent order
    % -----------------------------------------------------------
    species_data = {Hb_data, HbO_data, lipid_data, ICG_data, methyl_blue};

    % How many species are included?
    num_species = sum(species_bool);


    %% -----------------------------------------------------------
    % Build the absorption matrix
    % -----------------------------------------------------------
    A = zeros(num_species, num_points);

    count = 1;
    for n = 1:length(species_bool)
        if species_bool(n)
            curve = species_data{n};

            % Crop + interpolate
            selected_curve_interp = pick_bandwidth(curve, min_w, max_w, num_points);

            A(count, :) = selected_curve_interp;
            count = count + 1;
        end
    end
end
