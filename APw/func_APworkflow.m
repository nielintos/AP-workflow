function func_APworkflow(files,chnames,savedir,dotransf,namestransfsel2,...
    alg,n_red,class_str,nTimes,pthR)
% func_APworkflow(files,chnames,savedir,dotransf,namestransfsel2,alg,n_red,
%                   class_str,nTimes,k_neigh,distance,pthR)
% 
% Run Automated Phenotyping workflow for flow cytometry data.
% Comprises 5 stages: 
%     1.Pre-processing of FCM data to transform fluorescence signals into 
%       a linear scale
%     2.Subsampling of FCM events to facilitate suitable processing of 
%       full datasets and multiple samples
%     3.Dimensionality reduction of high-dimensional fluorescence signals 
%       for visualization purposes (tSNE)
%     4.Automatic clustering for unbiased detection of subpopulations
%       (Phenograph)
%     5.Upsampling of FCM events to extend results to the whole set of 
%       original events and files export
%
% The following parameters need to be specified as arguments:
%
%   Argument    Value
%   ----------  ----------------------------------------------------------
%   files       Array cell of Strings to specify location of fcs files to
%               analyze.
%
%   chnames     Array cell of Strings with names of channels (parameters)
%               included in the fcs files that will be used for analysis
%
%   savedir     String specifying the output folder to save final and
%               intermediate results
%
%   dotransf    Flag to specify is biexponential transformation will be
%               used (0: no transformation; 1: transform all files together
%               by using the same tranformation; 2: transform each file 
%               individually
%
%   namestransfsel2
%               Array cell where each cell is an array cell with strings
%               defining channel/parameter names that will be transformed
%               together with a single logicle/biexponential funtion
%
%   alg         Algorithm to use for downsampling ('LocalD': local density
%               subsampling; 'Random': random subsampling)
%
%   n_red       Number of events/cells to keep after downsampling
%
%   class_str   Strings with name of channel/parameter containing a
%               reference standard (e. g. manual gating id numbers). It can
%               be empty
%
%   nTimes      Number of runs for checking repeteability/variability of
%               results
%
%   pthR        String specifying the Path to R binaries
%
%
% OUTPUT        fcs files corresponding to the original ones (inputs) with
%               the inclusion of new parameters (tSNE coordinates, 
%               clustering results, and subsampling assignment).
%               Name for outputs: out_[filename]_i[N].fcs, where
%               [filename] is the original name of the fcs file, and [N]
%               is the number of repetition.
%
% 
% Version: January 2018 - CNIC
% Author:  Daniel Jimenez-Carretero
% Email:   djimenez at cnic dot es
%
% Copyright (c) 2018 CNIC
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%


%-----------------------
% PARAMETERS
%-----------------------
if ~exist('pthR','var') || isempty(pthR) || exist([pthR filesep 'R.exe'], 'file')==0
    disp('ERROR: invalid R path');
    return;
end
disp(['PATH TO R: ' pthR]);

pth=pwd;
if ~exist('files','var') || isempty(files) || any(cellfun(@(x) exist(x,'file')==0, files))
    disp('ERROR: invalid fcs files');
    return;
end
disp('FILES');
disp(files(:));
[pathnames,filenames,fileexts]=cellfun(@fileparts,files,'UniformOutput',0);

% Read first file to get parameter/channel names
[data,h,~,~]=fca_readfcs(files{1});
t=array2table(data, 'VariableNames', matlab.lang.makeValidName({h.par.name}));
t.File=repmat(1,size(t{:,1}));

if ~exist('chnames','var') || isempty(chnames) || ~all(ismember(chnames,t.Properties.VariableNames))
    disp('ERROR: invalid channels names');
    return;
end
disp('CHANNELS');
disp(chnames(:));

if ~exist('savedir','var') || isempty(savedir)
    disp('ERROR: invalid channels names');
    return;
elseif ~exist(savedir,'dir')
    mkdir(savedir);
end
disp(['OUTPUTDIR: ' savedir]);

if ~exist('dotransf','var') || isempty(dotransf) || length(dotransf)>1 || ~ismember(dotransf,0:2)
    disp('ERROR: invalid channels names');
    return;
end

if dotransf==0
    disp('TRANSFORMATION: NONE');
else
    if dotransf==1
        disp('TRANSFORMATION: All tubes together');
    else
        disp('TRANSFORMATION: Tubes independently');
    end
        
    savedir_biexp=[savedir filesep 'transf' filesep];
    % Check if transformed files already exist
    transf_fcs=strcat(savedir_biexp,'biexp_', filenames, fileexts);    
    if all(cellfun(@(x) exist(x,'file'),transf_fcs))
        disp('Files already transformed');
        namestransfsel2={};
    else
        % Set of channels to transform
        if ~exist('namestransfsel2','var') || isempty(namestransfsel2) || ~all(cellfun(@(x) all(ismember(x,t.Properties.VariableNames)),namestransfsel2)) 
            disp('ERROR: invalid sets of channels to transform');
            return;
        end
        disp('CHANNELS TO TRANSFORM:');
        cellfun(@(x) disp(x(:)), namestransfsel2)
    end
end

algnames={'Random','LocalD'};
if ~exist('alg','var') || isempty(alg) || ~ismember(alg, algnames)
    disp('ERROR: invalid subsampling algorithm');
    return;
end

if ~exist('n_red','var') || isempty(n_red) || length(n_red)>1 || ~isnumeric(n_red) || n_red<0
    disp('ERROR: invalid n_red argument');
    return;
end
k_clusters=n_red;
disp(['SUBSAMPLING ALGORITHM (N=' num2str(k_clusters) '): ' alg]);

if ~exist('class_str','var') || isempty(class_str)
    % No reference standard --> continue
    class_str=[];
    disp('No reference standard');
elseif ~isstr(class_str) || ~ismember(class_str, t.Properties.VariableNames)
    disp('ERROR: invalid name for CLASS parameter');
    return;
end
disp(['Reference Standard parameter: ' class_str]);

nT0=0;
if ~exist('nTimes','var') || isempty(nTimes) || length(nTimes)>1 || ~isnumeric(nTimes)
    disp('ERROR: invalid nTimes argument');
    return;
end
disp(['Repetitions: ' num2str(nTimes)]);

% Get coincident part of channels
 n=1;
 KK=chnames;
 while all(cellfun(@(x) strncmp(KK{1},x,n),KK(2:end)))
    n=n+1;
 end
 short_chnames=[];
 if n>1
    short_chnames=KK{1}(1:n-1);
 end

%-------------------------------
% TRANSFORMATION
%-------------------------------
disp('PREPROCESSING FILES...');
% FOLDER 4 TRANSF
if ~exist([savedir filesep 'transf'],'dir');
    mkdir([savedir filesep 'transf']);
end
biexp_files=strcat([savedir filesep 'transf' filesep 'biexp_'], strcat(filenames, fileexts));

% If they are not already created
if all(cellfun(@exist, biexp_files))
    disp('Biexponential transformation already computed');
else
    % If all together
    if dotransf==1
        for i=2:length(files)
            [data,h,~,~]=fca_readfcs(files{i});
            knam=matlab.lang.makeValidName({h.par.name});
            t_i=array2table(data, 'VariableNames', knam);
            t_i.File=repmat(i,size(t_i{:,1}));
            t=[t; t_i];
        end
        tout=transformdatabiexp(t,namestransfsel2);
        G=tout.File;
        tout.File=[];
        cellfun(@(x,f) fca_writefcs(f,tout{G==x,:},tout.Properties.VariableNames,tout.Properties.VariableNames), num2cell(1:length(files)), biexp_files);
    % Individually
    elseif dotransf==2
        for i=1:length(files)
            [data,h,~,~]=fca_readfcs(files{i});
            knam=matlab.lang.makeValidName({h.par.name});
            t_i=array2table(data, 'VariableNames', knam);

            disp(['Transforming file: ' filenames{i}]);
            tout=transformdatabiexp(t_i,namestransfsel2);
            fca_writefcs(biexp_files{i},tout{:,:},tout.Properties.VariableNames,tout.Properties.VariableNames);
        end
    elseif dotransf==0
        biexp_files=files;
    end
end

%-------------------------------
% Create Ref STD
%-------------------------------

RS={};
if ~isempty(class_str) & ~exist([savedir filesep 'RS.mat'],'file')
    disp('CREATING REFERENCE STANDARD...');
    for i=1:length(files)
        [data,h,~,~]=fca_readfcs(files{i});
        knam=matlab.lang.makeValidName({h.par.name});
        t_i=array2table(data, 'VariableNames', knam);
        RS{i}=t_i(:,class_str);
    end
    save([savedir filesep 'RS.mat'], 'RS');
end

[~,bi_filenames,bi_fileexts]=cellfun(@fileparts,biexp_files,'UniformOutput',0);
clear t t_i tout;

%-------------------------------
% REPETITIONS
%-------------------------------

for itN=nT0:(nT0+nTimes-1)
    disp(['REPETITION NUMBER ' num2str(itN)]);
    subname_red_file=[alg num2str(k_clusters) '_' short_chnames '_i' num2str(itN)];
    transf_biexp_fcs=strcat(savedir, filesep, bi_filenames, bi_fileexts);
    
    %-------------------------------
    % SUBSAMPLING
    %-------------------------------
    disp('     SUBSAMPLING...');
    savedir_red=[savedir, filesep, 'red'];
    if ~exist(savedir_red,'dir')
        mkdir(savedir_red);
    end
    transf_red_fcs=strcat(savedir_red, filesep, bi_filenames,'_',subname_red_file, bi_fileexts);
    idx_need_subsamp=find(cellfun(@(x) ~exist(x,'file'),transf_red_fcs));
    if isempty(idx_need_subsamp) & exist([savedir filesep 'IDX_' subname_red_file '.mat'],'file')
        disp('     All reduced files are created');
        load([savedir filesep 'IDX_' subname_red_file '.mat']);
    else
        % Read all of them in a different table
        t_i=[];
        t_o=[];
        for ii=idx_need_subsamp(:)'
            [data,h,~,~]=fca_readfcs(biexp_files{ii});
            knam=matlab.lang.makeValidName({h.par.name});
            tt_i=array2table(data, 'VariableNames', knam);
            clear data;
            t_i{ii}=tt_i{:,chnames};
            clear tt_i;
        end
        
        if exist([savedir filesep 'IDX_' subname_red_file '.mat'],'file')
            disp('     IDX files are loaded');
            load([savedir filesep 'IDX_' subname_red_file '.mat']);
        else
            IDX={};
            local_density={};
            % Check if it is a parallel pool
            poolobj = gcp('nocreate');
            if ~isempty(poolobj)
                delete(poolobj);
            end
            % Create best parpool
            maxjobs=length(idx_need_subsamp(:));
            if maxjobs<=8
                njobs=maxjobs;
            else
                num_jobs_2try=2:8;
                modules=mod(maxjobs,num_jobs_2try);
                njobs=num_jobs_2try(find(modules==min(modules),1,'last'));
            end
            parpool(njobs);
            
            parfor ii=idx_need_subsamp(:)'
                data=t_i{ii};
                if strcmp(alg,'LocalD')>0
                    % Codified from SPADE 2
                    % REFERENCE
                    %       Peng Qiu, Erin F. Simonds, Sean C. Bendall, 
                    %       Kenneth D. Gibbs Jr., Robert V. Bruggner, 
                    %       Michael D. Linderman, Karen Sachs, 
                    %       Garry P. Nolan, Sylvia K. Plevritis, 
                    %       "Extracting a Cellular Hierarchy from 
                    %       High-dimensional Cytometry Data with SPADE", 
                    %       Nature Biotechnology, 29(10):886-891, 2011.
                    kernel_width_factor=5;
                    density_estimation_optimization_factor=1.5;
                    rng('shuffle');
                    new_data = data(randi(size(data,1),2000,1),:);
                    % since new_data is part of data, this function does not compute min L1 dist. It computes the min non-zero distance, because there will be one 0. However, this creates a problem, what if the data simply contains a lot of identical entries, and result in multiple 0's? This code ignores them all, because this being zero does no good to the subsequent definination of neighborhood 
                    [min_dist,NN_ind] = compute_min_dist_downsample(new_data',data');
                    median_min_dist = median(min_dist);
                    kernel_width = median_min_dist*kernel_width_factor;
                    optimizaiton_para = median_min_dist*density_estimation_optimization_factor;
                    local_density{ii} = compute_local_density(data', kernel_width, optimizaiton_para);
                    target_density = downsample_to_certain_num_cells(data', local_density{ii}, k_clusters);
                    keep_prob = min(1,(target_density./local_density{ii}));
                    is_keep = rand(1,length(local_density{ii}))<keep_prob;
                    data2=data(is_keep,:);
                % Random
                elseif strcmp(alg,'Random')>0
                    rng('shuffle');
                    local_density{ii}={};
                    data2=data(randperm(size(data,1),k_clusters),:);
                else
                    disp('     WARNING: Subsampling algorithm is not valid... Using random...');
                    %alg='Random';
                    rng('shuffle');
                    local_density{ii}={};
                    data2=data(randperm(size(data,1),k_clusters),:);
                end
                % Now find the correspondences
                NN_index = zeros(1,size(data,1));
                block_size = 1000;
                ns = createns(data2,'nsmethod','kdtree','Distance','cityblock');
                for k=1:block_size:size(data,1)
                    ind_tmp = k:min([k+block_size-1,size(data,1)]);
                    [a] = knnsearch(ns,data(ind_tmp,:),'k',1);
                    NN_index(ind_tmp) = a';
                end 
                IDX{ii}=NN_index;
            end
            clear new data data2 is_keep keep_prob min_dist NN_index;
            save([savedir filesep 'LD_' subname_red_file '.mat'], 'local_density');
            save([savedir filesep 'IDX_' subname_red_file '.mat'], 'IDX');
        end

        disp('     SAVING DOWNSAMPLED FILES...');
        for ii=idx_need_subsamp(:)'
            [data,h,~,~]=fca_readfcs(biexp_files{ii});
            knam=matlab.lang.makeValidName({h.par.name});
            tt_i=array2table(data, 'VariableNames', knam);
            
            tt_i.(subname_red_file)=IDX{ii}(:);
            tt_i_Red=grpstats(tt_i,subname_red_file,{@(x) nanmean(x,1)});
            tt_i_Red.GroupCount=[];
            tt_i_Red=[tt_i_Red(:,2:end), tt_i_Red(:,1)];
            tt_i_Red.Properties.VariableNames=tt_i.Properties.VariableNames;
            
            fca_writefcs(transf_red_fcs{ii},tt_i_Red{:,:},tt_i_Red.Properties.VariableNames,tt_i_Red.Properties.VariableNames);
        end
        clear data tt_i tt_i_Red;
        
    end 
    
    disp('     TSNE AND CLUSTERING...');
    
    % LOAD ALL REDUCED FILES
    tt_i_Red_All=[];
    for ii=1:length(transf_red_fcs)
        [data,h,~,~]=fca_readfcs(transf_red_fcs{ii});
        knam=matlab.lang.makeValidName({h.par.name});
        tt_i_Red=array2table(data, 'VariableNames', knam);
        tt_i_Red.File=repmat(ii,size(tt_i_Red{:,1}));
        tt_i_Red_All=[tt_i_Red_All; tt_i_Red];
    end
    clear tt_i_Red data;
    
    % PHENOGRAPH
    if ~exist([savedir filesep 'tmp'],'dir')
        mkdir([savedir filesep 'tmp']);
    end
    
    % R clustering
    if exist([savedir filesep 'tmp' filesep 'out' '_i' num2str(itN) '.fcs'],'file')
        disp('     tSNE and clustering already computed');
        [datares,hres,~,~]=fca_readfcs([savedir filesep 'tmp' filesep 'out' '_i' num2str(itN) '.fcs']);
        knamres=matlab.lang.makeValidName({hres.par.name});
        t_res=array2table(datares, 'VariableNames', knamres);
    else
        wdir=pwd;
        cd([savedir filesep 'tmp']);
        
        % Overcome cluster-based problems derived from repetitions of the same event
        [~,ia,~]=unique(tt_i_Red_All(:,chnames),'stable');
        repeatedrows=setdiff(1:size(tt_i_Red_All(:,chnames),1), ia);
        if ~isempty(repeatedrows)
            tt_i_Red_All{repeatedrows,chnames}=(1e-6*rand([length(repeatedrows), size(tt_i_Red_All(:,chnames),2)]))+tt_i_Red_All{repeatedrows,chnames};
        end
        
        % R: preparation of inputs, execution and parsing of outputs
        fca_writefcs(['in_it' num2str(itN) '.fcs'],tt_i_Red_All{:,chnames},chnames,chnames);
        save('it.mat','itN');

        system(['"' pthR filesep 'R.exe" CMD BATCH ' wdir filesep 'APworkflow_Phenograph_tSNE.R log' '_i' num2str(itN) '.txt']);
        delete('it.mat');
        
        clfiles=dir(['*it' num2str(itN) '.mat']);
        cl_names=strrep({clfiles.name}, ['_it' num2str(itN) '.mat'],'');
        t_res=[];
        for kkl=1:length(clfiles)
            kkkk=load([clfiles(kkl).folder filesep clfiles(kkl).name]);
            kkkknames=fieldnames(kkkk);
            %kkkkdata=double(struct2array(kkkk));
			c = struct2cell(kkkk);
			kkkkdata = double([c{:}]);
            if strncmpi(cl_names{kkl},'tsne',4)
                t_res=[array2table(kkkkdata,'VariableNames',kkkknames),t_res];
            else
                t_res=[t_res,array2table(kkkkdata,'VariableNames',kkkknames)];
            end
        end
        
        % Find tsne
        knamres2=t_res.Properties.VariableNames;
        idx_notsne=find(cellfun(@(x) isempty(x), strfind(knamres2, 'tsne')));
        knamres2(idx_notsne)=strcat('CL_', knamres2(idx_notsne));
        knamres2=strcat(knamres2, ['_' short_chnames '_i' num2str(itN)]);
        
        t_res.Properties.VariableNames=knamres2;
        
        % Save
        fca_writefcs([savedir filesep 'tmp' filesep 'out' '_i' num2str(itN) '.fcs'],t_res{:,:},knamres2,knamres2);
        delete('*.RData');
        cd(wdir);
    end
    
    savedir_red_cl=[savedir filesep 'red_cl'];
    transf_red_cl_fcs=strrep(transf_red_fcs, savedir_red, savedir_red_cl);
    if exist([savedir_red_cl filesep 'CLS_i' num2str(itN) '.mat'],'file')
        disp('     tSNE and clusterings loaded');
        load([savedir_red_cl filesep 'CLS_i' num2str(itN) '.mat']);
    else
        % Save results from subsampled files
        cols_ok=~varfun(@(x) all(isnan(x)), t_res, 'OutputFormat','uniform');
        t_res2=t_res(:,cols_ok);

        tt_i_Red_All=[tt_i_Red_All t_res2];
        
        if ~exist(savedir_red_cl,'dir')
            mkdir(savedir_red_cl);
        end
        disp('     SAVING DOWNSAMPLED FILES WITH TSNE AND CLUSTERING RESULTS...');
        
        transf_red_cl_fcs=strrep(transf_red_fcs, savedir_red, savedir_red_cl);
        cls_aux=arrayfun(@(x) tt_i_Red_All(tt_i_Red_All.File==x,setdiff(tt_i_Red_All.Properties.VariableNames,'File','stable')),...
            1:length(transf_red_fcs), 'UniformOutput', false);
        cls=cellfun(@(x) x(:,[subname_red_file t_res2.Properties.VariableNames]), cls_aux, 'UniformOutput', false);
        save([savedir_red_cl filesep 'CLS_i' num2str(itN) '.mat'], 'cls');

        cellfun(@(x,f) fca_writefcs(f, x{:,:}, x.Properties.VariableNames, x.Properties.VariableNames), ...
            cls_aux(:), transf_red_cl_fcs(:));
        clear cls_aux;
    end
    if exist([savedir filesep 'CLS_' subname_red_file '.mat'],'file')
        disp('     Upsampled clustering results loaded');
        load([savedir filesep 'CLS_' subname_red_file '.mat']);
    else
        disp('     UPSAMPLING CLUSTERING RESULTS...');
        resampled_cls=cellfun(@(c,i) c(i,:), cls, IDX, 'UniformOutput', false);
        save([savedir filesep 'CLS_' subname_red_file '.mat'],'resampled_cls');
    end
    
    if exist([savedir filesep 'RS.mat'],'file') & ~exist([savedir filesep 'f1_' subname_red_file '.csv'],'file')
        disp('     COMPUTING F1 SCORES...');
        load([savedir filesep 'RS.mat']);
        RS_all=vertcat(RS{:});
        resampled_cls_all=vertcat(resampled_cls{:});
        idx_clch=strncmp(resampled_cls_all.Properties.VariableNames,'CL_',3);
        [t_grp_assign,t_allclustP,t_allclustR,t_allclust]=calculate_f1(RS_all,resampled_cls_all(:,idx_clch));
        writetable(t_allclust,[savedir filesep 'f1_' subname_red_file '.csv']);
    end
    
    savedir_resampit=[savedir filesep 'resamp_it'];
    if ~exist(savedir_resampit,'dir')
        mkdir(savedir_resampit);
    end
    % Generar los grandes
    if exist([savedir filesep 'CLS_TSNE_' subname_red_file '.mat'],'file')
        disp('     Upampled tsne loaded');
        load([savedir filesep 'CLS_TSNE_' subname_red_file '.mat']);
    else
        disp('     UPSAMPLING TSNE...');
        resampled_cls_tsne={};
        for i=1:length(files)
            [data,h,~,~]=fca_readfcs(transf_red_cl_fcs{i});
            knam=matlab.lang.makeValidName({h.par.name});
            t_red=array2table(data, 'VariableNames', knam);
            t_red_dat{i}=t_red{:,chnames};
            idx_tsne=find(cellfun(@(x) ~isempty(x), strfind(knam,'tsne_')));
            names_tsne{i}=knam(idx_tsne);
            t_red_tsne{i}=t_red(:,idx_tsne);

            [data,h,~,~]=fca_readfcs(biexp_files{i});
            knam=matlab.lang.makeValidName({h.par.name});
            t_biexp=array2table(data, 'VariableNames', knam);
            t_biexp_data{i}=t_biexp{:,chnames};
            clear data t_biexp;
        end
        resampled_cls_tsne=cell(1,length(files));
        poolobj = gcp('nocreate');
        if ~isempty(poolobj)
            delete(poolobj);
        end
        % Create best parpool
        if ~exist('njobs','var')
            maxjobs=length(files);
            if maxjobs<=8
                njobs=maxjobs;
            else
                num_jobs_2try=2:8;
                modules=mod(maxjobs,num_jobs_2try);
                njobs=num_jobs_2try(find(modules==min(modules),1,'last'));
            end
        end
        parpool(njobs);
        
        parfor i=1:length(files)
            Mdl = KDTreeSearcher(t_red_dat{i});
            
            [Idx,D] = knnsearch(Mdl, t_biexp_data{i}, 'K', 5);
            kkk=varfun(@(x) (1-x./sum(x))./sum(1-x./sum(x)), array2table(D'));
            D=kkk{:,:}';

            idx_tsne_resamp=find(cellfun(@(x) ~isempty(x), strfind(resampled_cls{i}.Properties.VariableNames,'tsne_')));
            tsne_test=resampled_cls{i}{:,idx_tsne_resamp};
            tsne_noisy=nan(size(t_biexp_data{i},1),2);

            t_red_tsne_neigh1=varfun(@(x) t_red_tsne{i}{x,1}, array2table(Idx'));
            t_red_tsne_neigh1=t_red_tsne_neigh1{:,:}';
            t_red_tsne_neigh2=varfun(@(x) t_red_tsne{i}{x,2}, array2table(Idx'));
            t_red_tsne_neigh2=t_red_tsne_neigh2{:,:}';


            a=mat2cell(t_red_tsne_neigh1,ones(size(t_red_tsne_neigh1(:,1))));
            b=mat2cell(t_red_tsne_neigh2,ones(size(t_red_tsne_neigh2(:,1))));
            c=mat2cell(tsne_test,ones(size(tsne_test(:,1))));
            d=mat2cell(Idx, ones(size(Idx(:,1))));

            [D2,Idx2]=cellfun(@(x,y,z) pdist2([x(:) y(:)], z,'euclidean','Smallest',3), a,b,c,'UniformOutput',0);

            max_dist_tsne=5;
            Idx2=cellfun(@(x,y) y(x<=max(min(x(:)),max_dist_tsne)), D2,Idx2,'UniformOutput',0);
            D2=cellfun(@(x) x(x<=max(min(x(:)),max_dist_tsne)), D2,'UniformOutput',0);
            % Con esto de despues podemos quitar lo del eps
            D_f=cellfun(@(x) (eps+double(sum(x)-x))/(eps+double(sum(sum(x)-x))),D2,'UniformOutput',0);
            
            Idx_f=cellfun(@(x,y) x(y), d, Idx2,'UniformOutput',0);

            tsne_noisy(:,1)=cellfun(@(x,y) sum(x.*t_red_tsne{i}{y,1}), D_f, Idx_f);
            tsne_noisy(:,2)=cellfun(@(x,y) sum(x.*t_red_tsne{i}{y,2}), D_f, Idx_f);

            resampled_cls_tsne{i}=[resampled_cls{i}(:,1:3), array2table(tsne_noisy,'VariableNames',strcat(names_tsne{i},'_noisy')), resampled_cls{i}(:,4:end)];
        end
        clear kkk t_red_tsne_neigh1 t_red_tsne_neigh2 tsne_test Idx a b c D2 d Idx;
        clear t_red_dat;
        clear t_red_tsne;
        clear t_biexp_data;
        save([savedir filesep 'CLS_TSNE_' subname_red_file '.mat'],'resampled_cls_tsne');
    end
    
    disp('     SAVING FINAL RESULTS...');
    for i=1:length(files)
        [a,b,c]=fileparts(files{i});
        if exist([savedir filesep b '_i' num2str(itN) c],'file')
            continue;
        end
        [data,h,~,~]=fca_readfcs(biexp_files{i});
        knam=matlab.lang.makeValidName({h.par.name});
        t_biexp=array2table(data, 'VariableNames', knam);
        clear data;

        [a,b,c]=fileparts(biexp_files{i});
        t_biexp_out=[t_biexp resampled_cls_tsne{i}];
        fca_writefcs([savedir_resampit filesep b '_i' num2str(itN) c], t_biexp_out{:,:}, t_biexp_out.Properties.VariableNames, t_biexp_out.Properties.VariableNames);
        clear t_biexp_out t_biexp;

        [data,h,~,~]=fca_readfcs(files{i});
        knam=matlab.lang.makeValidName({h.par.name});
        t_orig=array2table(data, 'VariableNames', knam);
        [a,b,c]=fileparts(files{i});
        t_orig_out=[t_orig resampled_cls_tsne{i}];
        fca_writefcs([savedir_resampit filesep b '_i' num2str(itN) c], t_orig_out{:,:}, t_orig_out.Properties.VariableNames, t_orig_out.Properties.VariableNames);
        clear t_orig t_orig_out;
        copyfile([savedir_resampit filesep b '_i' num2str(itN) c], [savedir filesep 'out_' b '_i' num2str(itN) c]);
    end

end
% Remove temporal files
outfiles=regexp(genpath(savedir),['[^;]*'],'match');
cellfun(@(x) rmdir(x,'s'), outfiles(2:end));
delete([savedir filesep '*.mat']);
disp('FINISHED');

   

       

