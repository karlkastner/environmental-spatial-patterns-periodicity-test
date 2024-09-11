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

	type_C = meta.type_C;
		
	for tdx=1:length(type_C)
		type = type_C{tdx};
	
		ifilename = sprintf([meta.filename.patterns_sampling_interval,'.shp'],meta.date_str,type);
		ofilename = sprintf(meta.filename.patterns_analyzed,meta.date_str,type);
		ofilename_stat_mat = sprintf([meta.filename.patterns_analyzed,'-stat'],meta.date_str,type);
	
		% TODO why is the spa here recreated from the shp instead of being
		%      loaded directly, is this to exclude invalid or skipped patterns?
		shp      = Shp.read(ifilename);
		centroid = Shp.centroid(shp);
		spa	 = Spatial_Pattern_Array(centroid,[],[shp.dx_sample]);
		spa.type = type;
	
		% assign paterns to world regions
		spa.assign_regions(meta.filename.region_shp);
	
		% analyze patterns, including the periodicity test
		spa.analyze();
		
		% export pattern analysis results to shape file
		spa.export_shp(ofilename);
	
		% save result to mat file
		save(ofilename,'-v7.3','spa');
	
		% save stat only
		spa.clear_1d_properties();
	
		% save stat to mat file
		save(ofilename_stat_mat,'-v7.3','spa');

		% link output shapefiles
		% prj projection file is not recreated and hence not copied 
		% this cannot be symlinks as otherwise git does not
		% store the file contents
		system(['cd output && ' ...
			'ln -f  ../mat/',meta.date_str,'/vegetation-patterns-*analyzed*.shp' ...
			'       ../mat/',meta.date_str,'/vegetation-patterns-*analyzed*.shx' ...
			'       ../mat/',meta.date_str,'/vegetation-patterns-*analyzed*.dbf' ...
			'       .' ...
		       ]);

		if (nargout() > 0)
			sp_a(idx) = spa
		end
	end % for tdx (pattern type)
	
end % function vegetation_patterns_analyze_2d

