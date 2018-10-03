%% The function converts all dat from numerous subfolders to csv to be exported to Python

constants;

data_folder = 'D:\\Experimental_Data\Transcription. New data from Madhav (2016_07)\';
csv_output_path = 'D:\\Experimental_Data\Transcription. New data from Madhav (2016_07)\all_data_new.csv';

gene_folders = {'HunchBack', 'Knirps', 'SNAIL'};
gene_folders_short = {'hb', 'kn', 'sn'};

% Initialize an empty structure
data = struct([]);
error_datasets = {};
AP_missing_datasets = {};
required_fields_list = {'OriginalParticle','Frame','Index','xPos','yPos','APpos', 'APPos', 'APposParticle', 'NuclearAP','MeanAP','MedianAP','DVpos','MeanDV','MedianDV','FirstFrame','Approved','Fluo','Off','Off2','Fluo2','FluoOld','FluoRaw','FluoError','FitType','nc','Nucleus','PParticle','DParticle','EParticle','NucStart','NucEnd','TotalmRNA','TotalmRNAError','dataset','construct','gene'};


% For a given gene get a subfolder list
g_inds = length(gene_folders);
trace_id = 0;
dataset_id = 0;
for g_ind = 1:g_inds
    gene_folder = fullfile(gene_folders{g_ind}, '_approved');
    datasets = dir(fullfile(data_folder, gene_folder));
    isub = [datasets(:).isdir]; %# returns logical vector
    datasets = {datasets(isub).name}';
    datasets(ismember(datasets,{'.','..'})) = [];
    sub_inds = length(datasets);
    % datasets

    % For a given subfolder load the ms2 data file
    data_sub = struct([]);
    for sub_ind = 1:sub_inds
        dataset = datasets{sub_ind}; 
        disp('Processing:');
        sub_path = fullfile(data_folder, gene_folder, dataset)
        
        try
            compiled_particles = load(fullfile(sub_path, 'CompiledParticles.mat'), 'CompiledParticles');
        catch
            error_datasets{end + 1} = fullfile(gene_folder, dataset);
            continue;
        end
        
        compiled_particles = compiled_particles.('CompiledParticles');
        
%         % Check if AP location is present
%         if 

        % Identify the dataset type
        if contains(dataset, 'shad', 'IgnoreCase', true)
            construct = 'no_sh';
            construct_id = 2;
        elseif contains(dataset, 'prim', 'IgnoreCase', true)
            construct = 'no_pr';
            construct_id = 1;
        else
            construct = 'bac';
            construct_id = 0;
        end

        % Get particles number
        particles = length(compiled_particles);
        
        data_part = struct([]);
        for particle = 1:particles
            fprintf('Processing gene %i/%i, dataset %i/%i, particle %i/%i\n', g_ind, g_inds, sub_ind, sub_inds, particle, particles);

            % Get frames number
            frames = length(compiled_particles(particle).Frame);

            data_fr = struct([]);
            for frame=1:frames
                % Expand frames to a flat structure
                line_struct = compiled_particles(particle);
                % Drop fields
                line_struct = rmfield(line_struct, {'optFit1', 'SlopeTrace', 'SDSlopeTrace'});

                field_names = fieldnames(line_struct);
                for field_ind = 1:length(field_names)
                    field=field_names{field_ind};
                    length(line_struct.(field));
                    if ~ischar(line_struct.(field)) && length(line_struct.(field)) > 1 
                        line_struct.(field) = line_struct.(field)(frame);
                    end
                end
                % Add & drop fields
                line_struct.dataset = dataset;
                line_struct.dataset_id = dataset_id;
                line_struct.construct = construct;
                line_struct.construct_id = construct_id;
                line_struct.gene = gene_folders_short{g_ind};
                line_struct.gene_id = g_ind - 1;
                line_struct.trace_id = trace_id;
                
        %         line_struct.fluo_per_polymerase = fluo_per_polymerase;
        %         line_struct.minutes_per_frame = orig_mins_per_frame;
                
                % Check for & add missing fields
                for field_ind = 1:length(required_fields_list)
                    field = required_fields_list{field_ind};
                    if ~isfield(line_struct, field)
                        line_struct.(field) = NaN;
                    end
                end


                % Calculate new fields
                line_struct.polymerases = line_struct.Fluo / fluo_per_polymerase;
                line_struct.time_min = line_struct.Frame * orig_mins_per_frame;

                % Add to the big structure
                data_fr = [data_fr, line_struct];
            end
            trace_id = trace_id + 1;
            data_part = [data_part, data_fr];
        end
        data_sub = [data_sub, data_part];
    dataset_id = dataset_id + 1;
    end
    data = [data, data_sub];
end

data;

% Ouput to file
data_table = struct2table(data);
% fid = fopen(csv_output_path,'w'); 
% fprintf(fid, data_table);
% fclose(fid);
writetable(data_table, csv_output_path, 'Delimiter', ';');

disp('Folders combined!');
count = length(error_datasets);
if count > 0
    disp('Encountered errors in the subfolders:');
    for i =1:count
        fprintf('%s\n', error_datasets{i});
    end
else
    disp('No errors while importing!')
end


















