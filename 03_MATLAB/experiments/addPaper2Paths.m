function addPaper2Paths()
%Put the Paper 2 source tree on the MATLAB path.
%Called at the top of every experiment script so that a run does not depend on
%whatever path state the session happens to be in. Placed ahead of the Paper 1
%archive so that Paper 2 versions of shared function names win.
%
%This is version 1.0 (Last edited: 2026-09-07)
%License: GPLv2.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(genpath(fullfile(root,'src')), '-begin');
addpath(here, '-begin');
if isfolder(fullfile(root,'configs')), addpath(fullfile(root,'configs'), '-begin'); end
end