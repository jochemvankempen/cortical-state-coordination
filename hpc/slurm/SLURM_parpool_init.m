function SLURM_parpool_init(setHPC, NumWorkers, IdleTimeout)
% initialise parpool environment on cluster (Rocket/SLURM).
%
% This function creates a new folder for each iteration of a parallel pool.
% This is necessary because when more than 1 Matlab session with parallel
% pool is called at the same time, there are communication issues (because
% Matlab assumes it has the machine (and filespace) to itself). There could
% thus be unwanted interaction between jobs that can cause errors, or worse,
% cross-talk between jobs. 
%
% Jochem van Kempen

%%% check whether parpool exists
pcheck = gcp('nocreate'); 
if ~isempty(pcheck)
    fprintf('parpool already running\n')
    return
end

%%% settings
defaultNumWorkers = 10;

%%% checks
if nargin<3 || isempty(IdleTimeout)
    IdleTimeout = 360; % set IdleTimeout to something large so that parallel pool doesn't get interupted and needs to start again.
end

if nargin<2 || isempty(NumWorkers)
    try
        NumWorkers = getenv('SLURM_CPUS_PER_TASK'); % try to read the number of workers set in sbatch script
        if isempty(NumWorkers)
            error('cannot read NumWorkers')
        end
        if ischar(NumWorkers)
            NumWorkers = str2double(NumWorkers);
        end
        fprintf('numWorkers (%d) read from getenv\n', NumWorkers)
    catch
        NumWorkers = defaultNumWorkers;
        fprintf('Could not read numWorkers, using default (%d) \n', NumWorkers)
    end
end

if nargin<1 
    global HPC
    setHPC = HPC;
end

%%% define and create directory for parpool
if setHPC
	%serverdir = getenv('TMPDIR') % scratch dir, somehow if I create directory there it doesn't show up
	%serverdir = [filesep 'scratch' filesep];
	serverdir   = [filesep '/mnt/nfs/home/']; % serverdir on rocket
    
    [arrayID, user, jobID] = SLURM_getEnv(1);

	if isempty(arrayID)
		arrayID = [filesep user filesep 'test' filesep]; % if run on interactive matlab, there's no array ID
		parpooldir = [serverdir arrayID]; % in this case, give the directory a different name
	else
		parpooldir = [serverdir user filesep 'tmpJobStorageLocation' filesep 'Job_' jobID '_' arrayID]; % set job specific folder
	end

	if ~exist(parpooldir, 'dir')
		[SUCCESS,~,~] = mkdir(parpooldir);
        if ~SUCCESS
            error('Cannot create parpooldir')
        end
	end

end

%%% start parpool
p = parcluster('local'); % open config for cluster
p.NumWorkers = NumWorkers;% set numworkers
if setHPC
    p.JobStorageLocation = parpooldir; % set the directory for this cluster
end

try
    pc = parpool(p); % open the parallel pool
catch
    
    try
        pc = parpool(p); % open the parallel pool
    catch
        fprintf('Could not load parpool \n')
    end
    
end

pc.IdleTimeout = IdleTimeout; % set idletimout to something big

