%% Download QAPLIB data
% This script downloads QAPLIB benchmark instances from 
% "https://coral.ise.lehigh.edu/data-sets/qaplib/"

QABLIB_url = 'http://coral.ise.lehigh.edu/wp-content/uploads/2014/07/qapdata.tar.gz';
QABLIB_soln_url = 'http://coral.ise.lehigh.edu/wp-content/uploads/2014/07/qapsoln.tar.gz';

if ~exist('data','dir')
    mkdir('data')
end

websave('./data/qapdata.tar.gz',QABLIB_url);
websave('./data/qapsoln.tar.gz',QABLIB_soln_url);
untar('./data/qapdata.tar.gz','./data/');
untar('./data/qapsoln.tar.gz','./data/');

% UPDATES: The library is slightly outdated, does not include some of the
% solutions or new results in the literature. We manually update the best 
% known solutions for: 
% esc32a, esc32b, esc32c, esc32d, esc32h, esc64a, kra32, tai100a
unzip('./data/qapsoln_update.zip','./data/');

fprintf('QAPLIB data Downloaded\n');