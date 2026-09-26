function setup_contourroots(varargin)
%SETUP_CONTOURROOTS Add ContourRoots to the MATLAB path for this session.
%   SETUP_CONTOURROOTS adds the public ContourRoots folders (matlab/ and
%   matlab/delay/) to the MATLAB path of the current session and prints a
%   short confirmation. Nothing is saved: the path returns to its previous
%   state when MATLAB restarts. Use SAVEPATH afterwards if you want the
%   change to persist.
%
%   SETUP_CONTOURROOTS('-quiet') does the same without printing.
%
%   You can call it from any folder, for example:
%       run("ContourRoots/setup_contourroots.m")
%
%   See also CROOTS, CPOLES, CZEROS, CPZMAP, CONTOURROOTS_VERSION.

    quiet = any(strcmpi(varargin,'-quiet'));
    root = fileparts(mfilename('fullpath'));
    folders = {fullfile(root,'matlab'), fullfile(root,'matlab','delay')};
    for k = 1:numel(folders)
        assert(isfolder(folders{k}), 'setup_contourroots:Missing', ...
            'Folder not found: %s. Is the download complete?', folders{k});
    end
    addpath(folders{:});
    if ~quiet
        fprintf('ContourRoots %s is ready (added to the path for this session).\n', ...
            contourroots_version());
        fprintf('Try:  r = croots(@(s) s.^2 + s + 1 + exp(-s), [-8 2 -20 20])\n');
        fprintf('Help: help croots   |   doc: %s\n', fullfile(root,'docs','index.md'));
    end
end
