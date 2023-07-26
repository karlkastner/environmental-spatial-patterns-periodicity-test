% Mon 31 May 20:20:46 CEST 2021
% Karl Kästner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% metadata
%
function meta = pattern_periodicity_test_metadata()
	% set pflag to 1 to save figures as files
	meta.pflag = 0;

	% for quick testing, set skip to a large integer to process only
	% every skip-patern in the database
	meta.skip = 1;

	% confidence level for significant frequency components 
	meta.significance_level_a1 = 0.05;
	
	meta.url.repository             = 'https://github.com/karlkastner/';

	meta.filename.patterns_shp      = 'input/vegetation-patterns-%s-selected-regions.shp';
	meta.filename.patterns_sampling_interval = 'mat/vegetation-patterns-%s-sampling-interval';
	meta.filename.patterns_analyzed = 'output/vegetation-patterns-%s-analyzed';
	meta.filename.dependencies      = 'dependencies.csv';
	meta.filename.profile           = 'mat/profiling-information.mat';
	meta.filename.region_shp        = 'input/regions-selected.shp';

	% required matlab toolboxes
	meta.toolbox_C = {
		'bioinformatics_toolbox',  'Bioinformatics Toolbox'
		'curve_fitting_toolbox',   'Curve Fitting Toolbox'
		'image_toolbox',           'Image Processing Toolbox'
		'map_toolbox',             'Mapping Toolbox'
		'signal_toolbox',          'Signal Processing Toolbox'
		'statistics_toolbox',      'Statistics and Machine Learning Toolbox'
		... % 'symbolic_toolbox',        'Symbolic Math Toolbox'
	};
end

