function [arrayID, user, jobID] = SLURM_getEnv(fid)
    
    if nargin==0
        fid=1;
    end

	arrayID     = str2double( getenv('SLURM_ARRAY_TASK_ID') ); % get current array id
	user        = getenv('USER'); % username
	jobID       = str2double( getenv('SLURM_JOBID') );% job ID
    
    if fid
        fprintf('Running SLURM array with arrayID:%d, user:%s and jobID:%d\n', arrayID, user, jobID)
    end
