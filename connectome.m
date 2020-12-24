% CONNECTOME static class with tools to input functional niftis, a mask/atlas, and
% output 2D connectomes (per nifti) in .mat format
%
% Methods
%   'batch'     - create job file to split across HPC tasks
%   'combine'   - combine outputs from 'batch'
%   'run'       - run analysis locally
%
% SurveyBott, 2018, info@surveybott.com
classdef connectome < handle
   properties
      data
      analysis
   end
   % load
   methods 
       function obj = connectome
           % setup properties
           fields.data = {'info' 'edges' 'coverage' 'inputs'};
           f = fieldnames(fields);
           for i=1:numel(fields)
              for j=1:numel(fields.(f{i}))
                 fields.(f{i}).(fields.(f{i}){j}) = []; 
              end
           end
       end
       function load(obj,input,varargin)
           % LOAD (input)
           %  Required 'Input' can be individual .mat or path to *_n*.mat batch outputs
           p = inputParser;
           p.addRequired('input',@(x) ischar(x) && (exist(x,'file') || isfolder(x)));
           p.parse(input,varargin{:});
           inputs = p.Results;
           if isfolder(input)
              out = obj.combine(input);
           else
              % single .mat
              out = load(input);
           end
           f = fieldnames(out);
           for i=1:numel(f)
               switch f{i}
                   case {'info' 'coverage'}
                       obj.data.(f{i}) = out.(f{i});
                   case {'edges','corr_clusters'}
                       obj.data.edges = out.(f{i});
               end
           end
           edges = out.edges;
           coverage = out.coverage;
           info = out.info;
       end
   end
   % connectome creation
   methods (Static)
       function jid = batch(funcPath,funcStr,savePath,saveStr,atlasFile,varargin)
           % BATCH(funcPath,funcStr,savePath,saveStr,atlasFile,...)
           %
           % INPUTS
           % Required
           %   'funcPath'      directory to run find command on
           %   'funcStr'       search string using find -wholename. e.g. *.nii
           %   'savePath'      location to save .mat output
           %   'saveStr'       string to save .mat as (_n# gets added)
           %   'atlasFile'      3D or 4D .nii or .nii.gz mask - assumes intensity values are 1:number of clusters
           %
           % Optional parameters
           %   'submit'        default: true
           %   'partition'     default: 'general'
           %   'cores'         default: 100
           %   'cpus'          default: 10
           %   'hours'         default: 6
           %   'mem'           default: 6000
           %
           %   'include'       cellstring to include (or file containing it)
           %       default:    {}
           %   'exclude'       cellstring to exclude (or file containing it)
           %       default:    {}
           %   'atlasVolume'    4th dimension volume to select. 1-indexed. Supports multiple - empty selects all
           %       default:    1
           %   'fishers'       Fisher's z transform correlation coeffs
           %       default:    true
           %   'crop':         two-elment vector of start point (1-indexed) and duration to crop time dim to
           %       default:    [1 -1]
           %   'scriptPath'    location of nii tools scripts
           %       default:    /ysm-gpfs/home/je45/scripts
           %   'JobStorageLocation'     tmp folder for parpool (to avoid collisions) 
           %       default":            '' (empty)
           
           % inputs
           p = inputParser;
           p.KeepUnmatched = true;
           p.addParameter('submit',true,@(x) islogical(x) || isnumeric(x));
           p.addParameter('partition','general',@ischar);
           p.addParameter('cores',100,@isnumeric);
           p.addParameter('cpus',10,@isnumeric);
           p.addParameter('hours',6,@isnumeric);
           p.addParameter('mem',6000,@isnumeric);
           p.addParameter('include',{},@(x) (ischar(x) && exist(x,'file')) || iscellstr(x));
           p.addParameter('exclude',{},@(x) (ischar(x) && exist(x,'file')) || iscellstr(x));
           p.parse(varargin{:});
           inputs = p.Results;
           
           % get scriptPath
           scriptPath = fileparts(mfilename('fullpath'));
           
           % find funcs
           [~,func] = unix(sprintf('find %s -wholename "%s"',funcPath,funcStr));
           func = strsplit(func,'\n');
           func(cellfun(@isempty,func)) = [];
           func = unique(func);
           if ~isempty(inputs.include)
               func_tmp = {};
               for i=1:numel(inputs.include)
                   func_tmp = [func_tmp func(~cellfun(@isempty,strfind(func,inputs.include{i})))];
               end
               func = func_tmp;
           end
           if ~isempty(inputs.exclude)
               for i=1:numel(inputs.exclude)
                   func(cellfun(@(x) ~isempty(regexp(x,inputs.exclude{i},'once')),func)) = [];
               end
           end
           
           % split into 'include' batches
           nodes = floor(inputs.cores/inputs.cpus);
           fprintf('%d runs split into %d array jobs (%d cores each)\n',numel(func),nodes,inputs.cpus);
           include = cell(1,nodes);
           idx = 1;
           for i=1:nodes
               % number of runs per batch
              if i > mod(numel(func),nodes)
                 n = floor(numel(func)/nodes);
              else
                 n = ceil(numel(func)/nodes); 
              end
              include{i} = func(idx:idx+n-1);
              idx = idx + n;
           end
           
           % setup and save inputs
           key = fieldnames(p.Unmatched);
           varargin = {};
           for i=1:numel(key)
               varargin{end+1} = key{i};
               varargin{end+1} = p.Unmatched.(key{i});
           end
           if ~isdir(savePath)
               mkdir(savePath);
           end
           inputFile = fullfile(savePath,sprintf('%s_inputs.mat',saveStr));
           save(inputFile,'funcPath','funcStr','savePath','saveStr','atlasFile','include','varargin');
           
           % write slurm file
           slurmFile = fullfile(savePath,sprintf('%s_batch.slurm',saveStr));
           fid = fopen(slurmFile,'w');
           fprintf(fid,'#!/bin/bash\n');
           fprintf(fid,'#SBATCH --job-name=connectome\n');
           fprintf(fid,'#SBATCH --partition=%s\n',inputs.partition);
           fprintf(fid,'#SBATCH --mem-per-cpu=%d\n',inputs.mem);
           fprintf(fid,'#SBATCH --time=%d:00:00\n',inputs.hours);
           fprintf(fid,'#SBATCH --cpus-per-task=%d\n',inputs.cpus);
           fprintf(fid,'#SBATCH --array=1-%d\n\n',nodes);
           
           fprintf(fid,'source $HOME/.bashrc\n');
           fprintf(fid,['matlab -nodisplay -r "addpath(''%s'');load(''%s'');idx=$SLURM_ARRAY_TASK_ID;include = include{idx};' ...
                        'connectome.run(funcPath,funcStr,savePath,sprintf(''%%s_%%d'',saveStr,idx),atlasFile,''include'',include,' ...
                        '''JobStorageLocation'',fullfile(savePath,''tmp'',num2str(idx)),varargin{:});exit"' ...
                        ],scriptPath,inputFile);
            fclose(fid);
           
            % submit
            jid = [];
            if inputs.submit
                owd = pwd;
                cd(savePath);
                [~,jid] = system(['sbatch ' slurmFile]);
                jid = strsplit(strtrim(jid));
                jid = jid{end};
                cd(owd);
            end
       end
       function out = combine(varargin)
           % COMBINE multiple .mat outputs from a batch run
           
           % inputs
           p = inputParser;
           p.addParameter('matPath',pwd,@isfolder);
           p.addParameter('matStr','*_n*.mat',@ischar);
           p.addParameter('volume',1,@isnumeric);
           p.addParameter('JobStorageLocation','',@ischar);
           p.addParameter('verbose',false);
           p.parse(varargin{:});
           inputs = p.Results;
           
           % outputs
           out.edges = [];
           out.coverage = [];
           out.info = [];
           f = fieldnames(out);
           
           % load mat batch files
           connectome.start_parpool('JobStorageLocation',inputs.JobStorageLocation);
           mat = dir(fullfile(inputs.matPath,inputs.matStr));
           edges = cell(size(mat));
           coverage = cell(size(mat));
           info = cell(size(mat));
           if ~isempty(mat)
               parfor i=1:numel(mat)
                   tmp = load(fullfile(mat(i).folder,mat(i).name));
                   for j=1:numel(f)
                       % select volume
                       if iscell(tmp.(f{j})) && numel(tmp.(f{j})) >= inputs.volume
                           if numel(tmp.(f{j})) >= inputs.volume
                               tmp.(f{j}) = tmp.(f{j}){inputs.volume};
                           else
                               error('Only %d output volumes and %d requested',numel(tmp.(f{j})),inputs.volume);
                           end
                       end
                   end
                   edges{i} = tmp.edges;
                   coverage{i} = tmp.coverage;
                   info{i} = tmp.info;
                   if inputs.verbose
                       fprintf('%d\n',i);
                   end
               end
           end
           
           % combine variables
           start = 1;
           for i=1:numel(mat)
               stop = start + size(edges{i},3) - 1;
               
               out.edges(:,:,start:stop) = edges{i};
               out.coverage(:,start:stop) = coverage{i};
               if i==1
                   out.info = info{i};
               else
                   out.info(start:stop) = info{i};
               end
               
               start = stop + 1;
           end
           
           % calculate corr_scans
           dims = size(out.edges);
           tmp_pos = tril(ones(dims(1),dims(2)),-1)==1; % lower triangle
           tmp_corr = reshape(out.edges,[dims(1)*dims(2) dims(3)]);
           tmp_corr = tmp_corr(tmp_pos,:);
           out.corr_scans = corr(tmp_corr,'rows','complete');
       end       
       function run(funcPath,funcStr,savePath,saveStr,atlasFile,varargin)
           % RUN(funcPath,funcStr,savePath,saveStr,atlasFile,...)
           %
           % INPUTS
           % Required
           %   'funcPath'      directory to run find command on
           %   'funcStr'       search string using find -wholename. e.g. *.nii
           %   'savePath'      location to save .mat output
           %   'saveStr'       string to save .mat as (_n# gets added)
           %   'atlasFile'      3D or 4D .nii or .nii.gz mask - assumes intensity values are 1:number of clusters
           %
           % Optional parameters
           %   'include'       cellstring to include (or file containing it)
           %       default:    {}
           %   'exclude'       cellstring to exclude (or file containing it)
           %       default:    {}
           %   'atlasVolume'    4th dimension volume to select. 1-indexed. Supports multiple - empty selects all
           %       default:    1
           %   'fishers'       Fisher's z transform correlation coeffs
           %       default:    true
           %   'crop':         two-column matrix of start point (1-indexed) and duration to crop time dim to. can be multiple rows
           %       default:    [1 -1]
           %   'scriptPath'    location of nii tools scripts
           %       default:    /ysm-gpfs/home/je45/scripts
           %   'JobStorageLocation'     tmp folder for parpool (to avoid collisions) 
           %       default":            '' (empty)
           % SurveyBott, 2016, info@surveybott.com
           
           % parse inputs
           p = inputParser;
           p.addParameter('atlasVolume',1,@(x) isnumeric(x));
           p.addParameter('crop',[1 -1],@(x) isnumeric(x) && size(x,2) == 2);
           p.addParameter('stat','mean',@(x) ismember(x,{'mean','pca'}));
           p.addParameter('scriptPath',fullfile(filesep,'ysm-gpfs','home','je45','scripts'),@isdir);
           p.addParameter('include',{},@(x) (ischar(x) && exist(x,'file')) || iscellstr(x));
           p.addParameter('exclude',{},@(x) (ischar(x) && exist(x,'file')) || iscellstr(x));
           p.addParameter('fishers',true,@(x) isnumeric(x) || islogical(x));
           p.addParameter('JobStorageLocation','',@ischar);
           p.parse(varargin{:});
           inputs = p.Results;
           
           atlasVolume = inputs.atlasVolume;
           crop = inputs.crop;
           stat = inputs.stat;
           fishers = inputs.fishers;
           scriptPath = inputs.scriptPath;
           include = inputs.include;
           exclude = inputs.exclude;
           clear p varargin
           
           if ~exist('funcStr','var') || ~ischar(funcStr)
               error('Please provide a ''funcStr'' to pass to find -wholename');
           elseif ~strcmp(funcStr(1),'*')
               funcStr = ['*' funcStr];
           end
           if ~exist('saveStr','var') || ~ischar(saveStr)
               saveStr = '';
           elseif strcmp(saveStr(1),'_')
               saveStr = saveStr(2:end);
           end
           if ischar(include) && exist(include,'file')
               tmp = load(include);
               if isfield(tmp,'include') && isstruct(tmp.include) && iscellstr(tmp.include)
                   include = tmp.include;
               else
                   error('%s needs to contain cellstr ''include''',include);
               end
               clear tmp
           elseif ~iscellstr(include)
               error('''include'' needs to be a cellstr or a file containing it');
           end
           if ischar(exclude) && exist(exclude,'file')
               tmp = load(exclude);
               if isfield(tmp,'exclude') && isstruct(tmp.exclude) && iscellstr(tmp.exclude)
                   exclude = tmp.exclude;
               else
                   error('%s needs to contain cellstr ''exclude''',exclude);
               end
               clear tmp
           elseif ~iscellstr(exclude)
               error('''exclude'' needs to be a cellstr or a file containing it');
           end
           if ~exist('atlasFile','var') || ~ischar(atlasFile) || ~exist(atlasFile,'file') || ...
                   (~strcmp(atlasFile(end-6:end),'.nii.gz') && ~strcmp(atlasFile(end-3:end),'.nii'))
               error('Please enter a valid .nii.gz or .nii ''atlasFile''');
           end
           if any(crop(:,1) < 1) || any(crop(:,2) ~= -1 & crop(:,2) < 1)
               error('One or more volume crops invalid');
           end
           
           % addpaths for scripts
           paths{1} = fullfile(scriptPath);
           paths{2} = fullfile(paths{1},'mri','NIFTI_tools');
           for i=1:numel(paths)
               addpath(paths{i});
           end
           clear scriptPath paths
           
           % load mask from 'atlasFile' - already validated to exist and have .nii.gz or .nii ending
           if strcmp(atlasFile(end-2:end),'.gz')
               gunzip(atlasFile);
               newatlasFile = regexprep(atlasFile,'\.gz','');
               mask = load_nii(newatlasFile);
               delete(newatlasFile);
               clear newatlasFile
           else
               mask = load_nii(atlasFile);
           end
           % select spatial map volume, reshape to 1D voxel index
           mask = mask.img;
           if ndims(mask) == 4
               if isempty(atlasVolume)
                   atlasVolume = 1:size(mask,4);
               end
               mask = mask(:,:,:,atlasVolume);
           elseif ndims(mask) < 3 || ndims(mask)  > 4
               error('Mask has %d dims. Should be 3 or 4.',ndims(mask));
           end
           dims = size(mask);
           mask = reshape(mask,[dims(1)*dims(2)*dims(3),numel(atlasVolume)]);
           clear atlasFile
           
           % make sure we're not asking for multiple atlasVolumes AND crops
           if numel(atlasVolume) > 1 && size(crop,1) > 1
               error('%d atlasVolumes and %d crops requested - only one can be multiple',numel(atlasVolume),size(crop,1));
           end
           
           % find funcs
           [~,func] = unix(sprintf('find %s -wholename "%s"',funcPath,funcStr));
           func = strsplit(func,'\n'); % parse unix find output
           func(cellfun(@isempty,func)) = []; % remove empties
           func = unique(func); % remove duplicates
           % include matching funcs
           if ~isempty(include)
               func_tmp = {};
               for i=1:numel(include)
                   func_tmp = [func_tmp func(~cellfun(@isempty,strfind(func,include{i})))];
               end
               func = func_tmp;
           end
           % exclude matching funcs
           if ~isempty(exclude)
               for i=1:numel(exclude)
                   func(cellfun(@(x) ~isempty(regexp(x,exclude{i},'once')),func)) = [];
               end
           end
           if isempty(func)
               error('No functionals found');
           else % create 'info' struct with (infoect)'id','scan', and 'cond' fields
               info = [];
               for i=1:numel(func)
                   info(i).file = func{i};
               end
           end
           clear path filename ext func_* info_* r
           
           % setup parpool to open with right num workers
           connectome.start_parpool('JobStorageLocation',inputs.JobStorageLocation);
           
           % funcs - gunzip, load, parcellate, corr matrix
           N = numel(func);
           fprintf('Processing %d runs\n',N);
           edges = cell(1,N); %zeros(max(mask),max(mask),N);
           coverage = cell(1,N); %zeros(max(mask),N);
           tmp_n = max([numel(atlasVolume) size(crop,1)]);
           parfor i=1:numel(func)
               fprintf('Processing run %2d...',i);
               
               % gunzip if necessary
               [tmp_path,tmp_filename,tmp_ext] = fileparts(func{i});
               if strcmp(tmp_ext,'.gz')
                   gunzip(func{i});
                   func{i} = fullfile(tmp_path,tmp_filename);
               end
               
               % load nii
               tmp_nii = load_nii(func{i});
               img = tmp_nii.img;
               dims = size(img);
               if numel(dims) == 3 % account for 1 vol in 4th dim
                   dims(4) = 1;
               end
               img = reshape(img,[dims(1)*dims(2)*dims(3),dims(4)]);
               
               % delete *.nii
               if strcmp(tmp_ext,'.gz')
                   delete(func{i});
               end
               
               % loop over crops/maskVols
               edges{i} = cell(1,numel(tmp_n));
               coverage{i} = cell(1,numel(tmp_n));
               for j=1:size(crop,1)
                   if crop(j,1) <= dims(4)
                       tmp_start = crop(1);
                   else
                       tmp_start = 1;
                       warning('Starting crop volume %d is greater than number of vols %d. Setting to %d.',crop(j,1),dims(4),tmp_start);
                   end
                   if crop(j,2) == -1
                       tmp_stop = dims(4);
                   elseif crop(j,1)+crop(j,2)-1 <= dims(4)
                       tmp_stop = crop(j,1)+crop(j,2)-1;
                   else
                       tmp_stop = dims(4);
                       warning('Crop duration %d from starting vol %d is greater than number of vols %d. Setting to max.',crop(j,2),crop(j,1),dims(4));
                   end
                   tmp_img = img(:,tmp_start:tmp_stop);
                   tmp_dims = dims;
                   tmp_dims(4) = size(tmp_img,2);
                   
                   % voxel-wise coverage
                   tmp_coverage_img = tmp_img;
                   tmp_coverage_img(isnan(tmp_coverage_img)) = 0;
                   tmp_coverage_img = range(tmp_coverage_img,2) ~= 0;
                   
                   % loop over parcellations, compute correlations and coverage
                   for k=1:numel(atlasVolume)
                       tmp_mask = mask(:,k);
                       % get mean timeseries for each cluster
                       tmp_clusters = zeros(max(tmp_mask),tmp_dims(4));
                       tmp_coverage = zeros(max(tmp_mask),1);
                       for m=1:max(tmp_mask)
                           switch stat
                               case 'mean'
                                   tmp_clusters(m,:) = mean(tmp_img(tmp_mask==m,:),1);
                               case 'pca'
                                   pc = pca(tmp_img(tmp_mask==m,:));
                                   tmp_clusters(m,:) = pc(:,1)';
                           end
                           % get percentage coverage per cluster
                           tmp_coverage(m) = sum(tmp_coverage_img(tmp_mask==m)) / sum(tmp_mask==m);
                       end
                       if numel(atlasVolume) == 1
                           tmp_idx = j;
                       else
                           tmp_idx = k;
                       end
                       coverage{i}{tmp_idx} = tmp_coverage;
                       edges{i}{tmp_idx} = corr(tmp_clusters');
                       if fishers
                           edges{i}{tmp_idx} = atanh(edges{i}{tmp_idx});
                       end
                   end
               end
               
               fprintf('done\n');
           end
           
           % reshape edges from {scans}{atlasVolume/crop}(cluster,cluster) to
           %   {atlasVolume}(cluster_pair,scans)
           tmp_corr = edges;
           tmp_coverage = coverage;
           edges = cell(1,tmp_n);
           coverage = cell(1,tmp_n);
           
           % setup corr_scans (scans,scans,atlasVolume/crop)
           corr_scans = zeros(numel(tmp_corr),numel(tmp_corr),tmp_n);
           for i=1:numel(edges)
               edges{i} = zeros(size(tmp_corr{1}{i},1)*size(tmp_corr{1}{i},2),numel(tmp_corr));
               coverage{i} = zeros(numel(tmp_coverage{1}{i}),numel(tmp_coverage));
               for j=1:size(edges{i},2)
                   edges{i}(:,j) = reshape(tmp_corr{j}{i},[size(edges{i},1) 1]);
                   coverage{i}(:,j) = tmp_coverage{j}{i};
               end
               tmp_pos = tril(ones(size(tmp_corr{1}{i},1),size(tmp_corr{1}{i},2)),-1)==1; % lower triangle
               tmp_pos = reshape(tmp_pos,[size(tmp_pos,1)*size(tmp_pos,2) 1]); % reshape for clarity
               corr_scans(:,:,i) = corr(edges{i}(tmp_pos,:),'rows','complete');
               edges{i} = reshape(edges{i},sqrt(size(edges{i},1)),sqrt(size(edges{i},1)),size(edges{i},2));
           end
           if numel(edges) == 1
               edges = edges{1};
               coverage = coverage{1};
           end
           
           if ~isdir(savePath)
               mkdir(savePath);
           end
           try
               save(fullfile(savePath,[saveStr '_n' num2str(N) '.mat']),'corr_scans','edges','coverage','info','inputs');
           catch
               save(fullfile(savePath,[saveStr '_n' num2str(N) '.mat']),'corr_scans','edges','coverage','info','inputs','-v7.3');
           end
       end
       function cleanup(outPath,varargin)
           % CLEANUP slurm output files generated during a batch run
           if ~exist('outPath','var')
              outPath = pwd; 
           end
           if ~isfolder(outPath)
              error('''outPath'' not a directory');
           end
           delete(fullfile(outPath,'slurm-*.out'));
           rmdir(fullfile(outPath,'tmp'),'s');
       end
   end
   % connectome manipulation
   methods
       % corr_scans
       
       % r to z
       
       % crad yeo
       
   end
   % connectome analysis
   methods
       % edge-wise
       
   end
   % compute util
   methods (Static)
       function pool = start_parpool(varargin)
           p = inputParser;
           p.addParameter('JobStorageLocation','',@ischar);
           p.parse(varargin{:});
           inputs = p.Results;
          
           cluster = parcluster;
           cluster.NumWorkers = feature('numCores');
           if ~isempty(inputs.JobStorageLocation)
              if ~isfolder(inputs.JobStorageLocation)
                  mkdir(inputs.JobStorageLocation)
              end
              cluster.JobStorageLocation = inputs.JobStorageLocation;
           end
           pool = gcp('nocreate');
           if ~isempty(pool) && (pool.NumWorkers < cluster.NumWorkers || ~strcmp(cluster.JobStorageLocation,pool.Cluster.JobStorageLocation))
               delete(pool);
               pool = gcp('nocreate');
           end
           if isempty(pool)
               pool = parpool(cluster,cluster.NumWorkers);
           end
       end 
   end
end