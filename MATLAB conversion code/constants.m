


%% Constants

orig_secs_per_frame = 37;  % in seconds
orig_mins_per_frame = orig_secs_per_frame/60;  % in seconds
fluo_per_polymerase = 3.5e4/100; % in Stucken.	Am I sure about this number?
forced_start_nc_time_13 = 35;
forced_start_nc_time_14 = 50;
% nuc_cyc = 14;
time_shift_limit = 5;  % in mins
k = 26*60;     % Bulk jump rate in bp/min
dk = 2 * 60;   % Error in k, bp/min
L = 5400;      % Gene length in base pairs
l = 45;        % Polymerase footprint in bp
dl = 5;         % Error in l
l_min = 40;     % Min polymerase footprint in bp
l_max = 60;     % Max polymerase footprint in bp
tau_exp = 5;    % in Seconds, for theoretical predictions
dtau_exp = 1;  % Error in tau, in seconds
L_array = [6444, 6444, 6444];       % Gene lengths in base pairs for Hb, Sn and Kn.
dL_array = [0, 0, 0];         % Error in gene length for Hb, Sn and Kn.
bootstrap_samples_count = 2000;			% How many times to get bootstrapped statistics
bootstrap_only_bin_value = 0.37;        % Only used for plots. Bootstrap calculated for all bins
bootstrap_only_each_frame = 5;  % 5, 8, 12
% N_filter_threshold = 15;
integral_threshold_value = 200;     % in pol * mins through one nuclear cycle. Filters out traces that do not get above a certain value. To reconsider
% filter_time_interval_mins = [48, 60];
init_slope_length = 4.0;		% Interval of measurement of the initial slope, in mins. Change to 4 min for the new gene length of 6444 bp. 
								% I am using two passes to sync correclty. Was 2.5 min
init_slope_length_frames = 5;								
raw_data_plot_every_which_point = 100;
data_folder = '/media/aserov/DATA/Experimental_Data/Transcription. New data from Madhav (2016_07)/';
output_figures_folder = './figures_for_article/';
bl_save_figures = true;
% bl_stop_early = true;
bl_stop_early = false;
half_width_max_rgn_ind = 2; % Number of points to consider around maximum when calculating the value of
                            % the plateau of concentration
input_folder = './processed_data/';
gene_names_array = {'HunchBack', 'SNAIL', 'Knirps'};
gene_short_names_array = {'Hb', 'Kn', 'Sn'};
dataset_names_array = {'bac', 'no_primary', 'no_shadow'};
dataset_short_names_array = {'bac', 'no_pr', 'no_sh'};
nuc_cyc_array = [13, 14];

hist_color = [117, 145, 41] / 255;

% Refined time mesh
desired_time_step = 5/60;		% in mins
new_time_start = 0;				% in mins
new_time_end = 75;				% in mins



%% Manually set min and max AP limits
manual_bins=true;
if manual_bins
    min_ms2_AP = 0.15;
    max_ms2_AP = 0.75;
end;



% Defining colors
% my_color_sequence = get(groot,'defaultAxesColorOrder');
my_bins_color_sequence = [...
    0 0.251 0;...       %	Dark green       
    0 0 1;...           %	Blue
    0.9412 0.4706 0;... %   Orange
    0.502 0.251 0;...   %	Brown                    
    0 0.502 0.502;...   %	Turquoise
    1 0 0;...           %	Bright red
    1 0.502 0.502;...   %	Peach
    0 1 1;...           %	Cyan
    0.502 0.502 1;...   %	Light purple
    0 1 0;...           %	Bright green
    0.502 0 0;...       %	Burgundy 
    1 0 1;...           %	Pink
    0.251 0 0.502;...   %	Purple
    1 1 0;...           %	Yellow                    
    0 0 0;...           %	Black
    0 0.251 0;...       %	Dark green       
    0 0 1;...           %	Blue
    0.9412 0.4706 0;... %   Orange
    0.502 0.251 0;...   %	Brown                    
    0 0.502 0.502;...   %	Turquoise
    1 0 0;...           %	Bright red
    1 0.502 0.502;...   %	Peach
    0 1 1;...           %	Cyan
    0.502 0.502 1;...   %	Light purple
    0 1 0;...           %	Bright green
    0.502 0 0;...       %	Burgundy 
    1 0 1;...           %	Pink
    0.251 0 0.502;...   %	Purple
    1 1 0;...           %	Yellow                    
    0 0 0];...          %	Black
                    
my_constructs_color_sequence = [...
    0.502 0 0;...       %	Burgundy 
    1 0 1;...           %	Pink
    0.251 0 0.502;...   %	Purple
    1 1 0];             %	Yellow

% % Setting default colors
% set(groot,'defaultAxesColorOrder', my_bins_color_sequence);



