function info=time_build_step(name,operation,outputDir)
%TIME_BUILD_STEP Time a build stage without skipping or changing its work.
%   Writes one JSON per stage to output/build_timings by default. Failed
%   stages record their timing and rethrow the ORIGINAL exception.
    validateattributes(name,{'char'},{'row','nonempty'});
    if isempty(regexp(name,'^[A-Za-z][A-Za-z0-9_-]*$','once'))
        error('time_build_step:Name','Use an alphanumeric stage name.');
    end
    if nargin<3
        root=fileparts(fileparts(mfilename('fullpath')));
        outputDir=fullfile(root,'output','build_timings');
    end
    info=struct('stage',name,'started',char(datetime('now')), ...
        'matlab',version,'seconds',NaN,'status','running','errorIdentifier','');
    fprintf('\nBUILD starting %s...\n',name); timer=tic;
    try
        operation(); info.status='passed';
    catch err
        info.status='failed'; info.errorIdentifier=err.identifier; record(); rethrow(err)
    end
    record();
    function record()
        info.seconds=toc(timer);
        fprintf('BUILD %s: %s in %.3f s.\n',name,info.status,info.seconds);
        % Timing is observational: a logging failure must neither prevent
        % the operation nor replace its original exception.
        try
            if ~isfolder(outputDir)
                [ok,message]=mkdir(outputDir);
                if ~ok, error('time_build_step:Directory','%s',message); end
            end
            fid=fopen(fullfile(outputDir,[name '.json']),'w');
            if fid<0, error('time_build_step:Write','Cannot open the timing report.'); end
            closer=onCleanup(@() fclose(fid)); fprintf(fid,'%s\n',jsonencode(info,'PrettyPrint',true));
        catch logError
            warning('time_build_step:Log','Could not write timing log for %s: %s',name,logError.message);
        end
    end
end
