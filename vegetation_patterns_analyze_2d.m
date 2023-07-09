% 2022-09-26 14:32:10.449630621 +0200
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
%
%% analyze the satellite images fetched from the google satellite server
%% Note: run  devilspie to avoid showing figures in foreground while processing
function sp_a = vegetation_patterns_analyze_2d(meta)
	if (nargin()<1)
		meta = pattern_periodicity_test_metadata();
	end

	type_C = {'anisotropic','isotropic'};
		
	for tdx=1:length(type_C)
	type = type_C{tdx};

	iname = sprintf([meta.filename.patterns_sampling_interval,'.shp'],type);
	obase = sprintf(meta.filename.patterns_analyzed,type);

	shp      = Shp.read(iname);
	centroid = Shp.centroid(shp);
	% for a quick check, choose only a subset of patterns, such as 1:10:l
	fdx      = 1:length(shp);
	spa	 = Spatial_Pattern_Array(centroid(fdx,:),[shp(fdx).dx_sample]);
	spa.type = type;

	% assign paterns to world regions
	spa.assign_regions(meta.filename.region_shp);

	% analyze patterns, including the periodicity test
	spa.analyze();

	% export pattern analysis results to shape file
	spa.export_shp([obase,'.shp']);


	% save result to mat file
	save([obase,'.mat'],'-v7.3','spa');

	% save stat only
	for idx=1:length(spa.sp_a)
		spa.sp_a(idx).clear_1d_properties();
	end

	% save result to mat file
	save([obase,'-stat.mat'],'-v7.3','spa');

	if (nargout() > 0)
		sp_a(idx) = spa
	end
	end % for tdx (pattern type)
	
end % function

