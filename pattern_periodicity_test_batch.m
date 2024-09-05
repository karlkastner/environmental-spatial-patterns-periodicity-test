% Fri 30 Jun 10:49:17 CEST 2023
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
%% batch script for reproducing the analysis, figures and tables
%
function pattern_periodicity_test_batch()
	% create library and output folder
	mkdir('./lib/');
	mkdir('./mat/');
	mkdir('./img/');
	mkdir('./output/');
	mkdir('./lib/auxiliar/');
	addpath('./lib/auxiliar');

	meta = pattern_periodicity_test_metadata();

	% check if required matlab toolboxes are installed
	for idx=1:size(meta.toolbox_C,1)
		if (~license('test',meta.toolbox_C{idx,1}))
			printf('%s is missing, execution will likely fail at a later point.\',meta.toolbox_C{idx,2});
		end
	end

	% fetch the script for fetching library files
	url  = 'https://raw.githubusercontent.com/karlkastner/auxiliar/master/dependencies_fetch.m';
	dest = './lib/auxiliar/dependencies_fetch.m';
	urlwrite(url,dest);

	% fetch library files
	% dependencies_determine(meta.filename.dependencies,meta.filename.profile,{'pattern_periodicity_test_batch'}) 
	dependencies_fetch(meta.url.repository,meta.filename.dependencies);

	% add library subfolders to path
	addpath_recursive('./lib');

	% set to true to save fitures to files 
	pflag      = false;
	meta.pflag = pflag;

	% minimum working example, note that dependencies have to be fetched first
	pattern_periodicity_test_minimum_example();

	% Figure 1: manmade vs natural patterns
	pattern_manmade_vs_natural_plot(meta);
	
	% Figure 2, 3: periodicity test 
	example_periodicity_test_amendments(pflag);
	
	% Figure 4b: example pattern delineation
	example_pattern_delineation_plot(meta);

	% Figure 4c: example tested frequency range
	example_frequency_range_of_test(pflag);

	% section 2.1: example of Renshaw & Fords (1984) test
	example_periodicity_test_renshaw_ford();

	% fetching of satellite images of patterns in the global dataset
	vegetation_patterns_fetch();

	% analying of patterns in the global dataset
	vegetation_patterns_analyze_2d(meta);

	% Figure 5: plot results of the periodicity test of patterns in the dataset
	vegetation_patterns_plot_periodicity_test_results(meta);

	% Figure 6: example of periodic and stochastic patterns
	pattern_synthetic_periodic_plot(meta);

	% Figure
	vegetation_patterns_plot_joint_density_and_acf();

	% Figure
	vegetation_patterns_plot_goodness_of_fit();

	% Figure SI 1:
	experiment_patch_variation();

	% Figure SI 2:
	experiment_acf_bias();

	% Figure SI 3: 
	example_detection_of_periodic_components_in_stochastic_patterns(pflag);

	% Figure SI 4:
	experiment_periodicity_test_2d_bias();

	% for web-page generation run vegetation-patterns.rmd in R

end % function pattern_periodicity_test_batch

