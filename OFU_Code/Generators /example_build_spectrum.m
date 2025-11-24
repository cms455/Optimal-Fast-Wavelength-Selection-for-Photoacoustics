data_folder = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data';

species_bool = [1 1 1 1 0];   % Include Hb, HbO, Lipid, ICG only

A = build_absorption_matrix(700, 900, species_bool, 200, data_folder);
