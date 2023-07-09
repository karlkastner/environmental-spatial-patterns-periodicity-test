% Mon 15 May 16:22:50 CEST 2023
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
%% plot an example pattern delineated by a polygon
%
function example_pattern_delineation_plot(meta)
	if (nargin()<1)
		meta = pattern_periodicity_test_metadata();
	end
	pflag = meta.pflag;

	% coordinate of pattern	centrold
	X0 = [07.79111,+047.92607];
	% load pattern polygons
	shp = Shp.read(sprintf(meta.filename.patterns_shp,'anisotropic'));
	centroid = Shp.centroid(shp);
	% determine pattern closest to centroid
	d = hypot(centroid(:,1)-X0(2),centroid(:,2)-X0(1));
	%d = hypot([shp.X]-X0(2),[shp.Y]-X0(1));
	[~,mdx] = min(d);
	
	% reproject coordinate system	
	proj    = projcrs(900913);
	[xp,yp] = projfwd(proj,shp(mdx).Y,shp(mdx).X);
	s = 1e-3;
	[xc,yc] = projfwd(proj,X0(1),X0(2));

	% read image file
	g = GeoImg();
	g.read('input/pattern_2d_+07.79111_+047.92607.png');
	b = g.img;	

	b   = imnormalize(b,0.05,0.999);
	x = g.x;
	y = g.y;

	figure(1);
	clf();
	imagesc(s*(x-xc),s*(y-yc),b);
	hold on
	plot(s*(xp-xc),s*(yp-yc),'r','linewidth',3)
	axis equal
	axis xy
	xlim(s*[-1,1]*2001);
	ylim(s*[-1,1]*2001);
	xlabel('Easting / km');
	ylabel('Northing / km');
	set(gca,'xtick',(-2:2));
	set(gca,'ytick',(-2:2));
	colormap gray;
	clim(quantile(b(:),[0.05,0.95]));
	
	if (pflag)
		ps = 2.5;
		pdfprint(1,'img/pattern-observed-polygon.pdf',ps);
	end
	
end % function
