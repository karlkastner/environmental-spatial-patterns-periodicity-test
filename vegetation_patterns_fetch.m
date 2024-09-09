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
%% determine the sampling interval for fetching images from the Google satellite
%% server and later processing
%%
function vegetation_patterns_fetch(meta)
	if (nargin()<1)
		meta = pattern_periodicity_test_metadata();
	end

	type_C = meta.type_C;

	% determine if higher frequency spectrum of sampled patterns is well resolved
	% for anisotropic and isotropic patterns
	for tdx=1:length(type_C)
		type = type_C{tdx};

		iname = sprintf(meta.filename.patterns_shp,type);
		obase = sprintf(meta.filename.patterns_sampling_interval,type);

		spa = Spatial_Pattern_Array();
		spa.type = type;
		% for a quick testing: when skip is set to an integer larger 1,
		% only every skip-pattern is processed
		spa.opt.skip = meta.skip;

		% axuiliary file with hash of spectral quantiles
		spa.fetch(iname,[obase,'.shp']);

		% save result to mat file
		save([obase,'.mat'],'-v7.3','spa');
	end % for type

end % function vegetation_patterns_fetch

