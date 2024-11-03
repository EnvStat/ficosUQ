function [y] = evaluate_design(x,sim_f,post_f,tt_data,...
        useLogOfLogPosterior,outputPath, options)
   % EVALUATE_DESIGN Evaluates simulator for a set of design points.
   % y = evaluate_design(x, sim_f, post_f) Runs simulator at given design 
   % points x and returns the log posterior densities for each point. Function
   % handles sim_f and post_f define the functions for the simulator and the
   % negative log posterior densities, respectively.
   %
   % Used by various scripts in the calibration_main pipeline with additional
   % internal arguments.
   %
   % Copyright (c) 2017 - 2024 Karel Kaurila
   %
   arguments
       x {mustBeNumeric, mustBeReal}
       sim_f {mustBeA(sim_f, 'function_handle')}
       post_f {mustBeA(post_f, 'function_handle')}
       tt_data = [];
       useLogOfLogPosterior (1,1) {mustBeNumericOrLogical} = false;
       % (optional) where to save simulator outputs
       % if set to false, don't save outputs
       outputPath = false;
       % Additional options given as (name, value) pairs
       % --------------------------------------------------------
       % Whether to use old version of this function ('legacy') or the
       % current ('new') version
       options.mode {mustBeTextScalar, ...
           mustBeMember(options.mode, {'legacy', 'new'})} = 'new';
       % The remaining arguments are exclusive to the 'legacy' mode:
       % --- legacy arguments -----
       options.locations {mustBeMember(options.locations, {'intens'})} = ...
           'intens';
       options.log_q (:,1) {mustBeNumericOrLogical} = false;
       options.keepInputs (1,1) {mustBeNumericOrLogical} = false;
       options.max_cores (1,1) {mustBeInteger,mustBePositive} = 20;
       % Paths needed for attaching simulator files --------
       % path to folder containing simulator executable and .ini files
       options.simDir {mustBeFolder} = fullfile('bin/');
       % names of the simulator executable and .ini files
       options.simFile {mustBeTextScalar} = 'wqficos';
       options.simIniFile {mustBeTextScalar} = 'wqficos.ini';
       % path to folder containing .hdf5 input files
       options.hdf5Dir {mustBeFolder} = fullfile('data/');
       % names of the hdf5 input files
       options.hdf5loadName {mustBeTextScalar} = 'loading.hdf5';
       options.hdf5fileName {mustBeTextScalar} = 'hd_files.hdf5';
   end

   %% Input parsing

   log_q = options.log_q;
   predictions = false;
   posterior = false;
   keepInputs = options.keepInputs;
   max_cores = options.max_cores;

   dim_x = size(x,2);
   if any(outputPath)
       % Saving simulator output for predictive distribution
       predictions = true;
       fprintf('Saving predictions to %s\n', outputPath);
       wf_ID = [categorical(1900000274); ... % Brändö
                categorical(1900000089); ... % Utö
                categorical(1900000099)]; % Seili
       if(isnumeric(log_q) && ~isscalar(log_q))
           posterior = true;
       else
           posterior = false;
           log_q = zeros(size(x,1),1);
       end
   else
       wf_ID = [];
       log_q = zeros(size(x,1),1);
   end
   %% Evaluating simulator at design points

   switch options.mode
       case 'legacy'
           % Legacy implementation for compatibility
           filesToAttach = [fullfile(options.simDir, ...
               {options.simFile, options.simIniFile}), ...
               fullfile(options.hdf5Dir, ...
               {options.hdf5fileName, options.hdf5loadName})];
           try
               n_cores = feature('numcores');
           catch ME
               warning('Cannot use feature numcores.')
           end
           n_cores = min(n_cores, max_cores);
           poolobj = gcp('nocreate');
           if( isempty(poolobj))
               poolobj = parpool([2,n_cores]);
           elseif(poolobj.NumWorkers < n_cores)
               delete(poolobj)
               poolobj = parpool([2, n_cores]);
           end
           % Check if necessary files attached:
           ficos_attached = any(endsWith(poolobj.AttachedFiles,...
               options.simFile));
           hd_files_attached = any(endsWith(poolobj.AttachedFiles, ...
               options.hdf5fileName));
           loading_attached = any(endsWith(poolobj.AttachedFiles, ...
               options.hdf5loadName));
           ficos_ini_attached = any(endsWith(poolobj.AttachedFiles,...
               options.simIniFile));
           if( ~(ficos_attached && ...
                   hd_files_attached && ...
                   loading_attached && ...
                   ficos_ini_attached))
               addAttachedFiles(poolobj,filesToAttach);
           end
           tic
           n_design = size(x,1);
           Y = zeros(n_design,1);
           log_w = zeros(n_design,1);
           F_pred = {};
          
           parfor i=1:n_design
              w = getCurrentWorker;
              wid = w.ProcessId;
    
              f = sim_f(x(i,:), i, wid);
    
              Y(i) = post_f(x(i,:),f);
    
              if(predictions)
                  f.polyID = categorical(f.polyID);
                  G = groupfilter(f,'Aika',@(z)ismember(z,wf_ID),'polyID');
                  sampleID = i*ones(size(G,1),1);
                  F_pred{i} = addvars(G, sampleID);
                  if(posterior)
                      % Note the sign; log_q and Y are negative log densities, i.e.
                      % high value means smaller density
                      log_w(i) = log_q(i)-Y(i);
                  end
              end
              if(useLogOfLogPosterior)
                  Y(i) = log(Y(i));
              end
              fprintf('h(x_%d = [', i);
              fprintf('%1.3g, ', x(i,1:(dim_x-1)));   
              fprintf('%1.3g]) = %1.5g\n', x(i,dim_x), Y(i));
           end
           y = reshape(Y',n_design,1);
           toc
           if(predictions)
               save(outputPath,'F_pred','log_w');
           end
        
           if(~keepInputs)
               rmdir('input/w*', 's');
           end
           
           parfevalOnAll(@clearvars, 0);
       case 'new'
           % New implementation
           simOut = sim_f(x);
           y = NaN(size(x,1), 1);
           for i = 1:size(x,1)
               % cleanSimOutput ensures simulator outputs are in the format
               % expected by the posterior density function
               y(i) = post_f(x(i,:), cleanSimOutput(simOut{i}, 'old'));
           end
           if useLogOfLogPosterior
               y = log(y);
           end
       otherwise
       error('Invalid mode: %s', options.mode);
   end
end