% 2023-06-02 14:43:13.373394646 +0200
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
%% plot natural and manmade patterns
%
function sp = pattern_manmade_vs_natural_plot(meta)
	if (nargin()<1)
		meta = pattern_periodicity_test_metadata();
	end
	pflag  = meta.pflag;
	fflag  = pflag;
	folder = 'input';
	p      = 0.375;
	cmap   = flipud((1-p)*gray()+p*colormap_vegetation(256));
	file_C = {
		  'nature-anisotropic-mexico_-103.147245_27.656891.png'
		, 'nature-isotropic-mexico_-107.158326_31.337108.png'
		, 'plantation-striped-spain_-2.685893_37.565588.png'
		, 'plantation-hexagonal-spain_-5.356567_37.401112.png'
	};
	l = [4,4,3,3];

	% analyze patterns
	tab = table();
	sp = repmat(Spatial_Pattern(),length(file_C),1);
	for idx=1:length(file_C)
		% note, without this line, all elements in sp are kept equal
		sp(idx) = Spatial_Pattern();
		disp(['Processing pattern ', num2str(idx)]);
		sp(idx).imread([folder,filesep,file_C{idx}]);
		sp(idx).analyze_grid();

		q1 = 0.025;
	        q2 = 0.975;
		% read satellite image of pattern
		b = imread([folder,filesep,file_C{idx}]);
		b = imnormalize(b,q1,q2);
		if (idx ~= 3)
	        	sp(idx).b_square = b;
		else
	        	sp(idx).b_square = imrotate(b,90);
		end

		% pattern statistics
		tab.filename{idx} = file_C{idx};
		tab.p_periodic(idx) = sp(idx).stat.p_periodic;
		tab.isisotropic(idx) = sp(idx).stat.isisotropic;
		if (sp(idx).stat.isisotropic)
			tab.wavelength(idx) =  1./sp(idx).stat.fc.radial.con;
			tab.L_rel(idx)      =  sp(idx).L(1).*sp(idx).stat.fc.radial.con;
		else
			tab.wavelength(idx) =  1./sp(idx).stat.fc.x.con;
			tab.L_rel(idx)      = sp(idx).L(1).*sp(idx).stat.fc.x.con;
		end
		tab.intS_hp_sig(idx) = sp(idx).stat.stati.intS_hp_sig;
	end % for idx

	% plot patterns
	for idx=1:length(sp)
		% plot pattern
		splitfigure([2,2],[1,idx],fflag);
		cla();
		%sp(idx).b_square = imread([folder,filesep,file_C{idx}(1:end-4),'-contrast-stretch.png']);
		sp(idx).plot('b');
		if (idx~=1)
			xlim([0,10]);
			ylim([0,10]);
		end

		% plot periodogram
		splitfigure([2,2],[2,idx],fflag);
		cla();
		% the dots are otherwise too tiny on the plot
		if (3 == idx || 4 == idx)
				% filtering with nf=3 reduced the peak height to 1/4,
				% this is compensated here for by  multiplying with 4
				sp(idx).S.con = 4*ifftshift(trifilt2(fftshift(sp(idx).S.con),3));
				sp(idx).S.rot.con = 4*ifftshift(trifilt2(fftshift(sp(idx).S.rot.con),3));
		end
		if (idx == 4)
			c = sp(idx).plot('S.con');
		else
			c = sp(idx).plot('S.rot.con');
		end
		cl = clim();
		clim(0.5*cl);
		axis([-1,1,-1,1]*l(idx));
		%c = colorbar();
		title(c,'$\hat S/\lambda_c^2$','interpreter','latex');
		set(c,'location','east');
		colormap(cmap)

		% plot correlogram
		splitfigure([2,2],[3,idx],fflag);
		cla();
		%sp(idx).b_square = imread([folder,filesep,file_C{idx}(1:end-4),'-contrast-stretch.png']);
		if (idx == 4)
			c=sp(idx).plot('R.con');
		else
			c=sp(idx).plot('R.rot.con');
		end
		if (1 == idx)
			%caxis([-0.2413,0.66]);
			caxis(0.3*[-1,2])
		end
		if (idx==2)
			caxis(0.1*[-1,2])
			%caxis([-0.0627,0.5]);
		end
		if (idx==1 || idx == 2)
			xlim(1.5*[-1,1]);
			ylim(1.5*[-1,1]);
		else
			xlim(2.5*[-1,1]);
			ylim(2.5*[-1,1]);
		end
		colormap(cmap);
		%c = colorbar();
		title(c,'$\hat R$','interpreter','latex');
		set(c,'location','east');

		if (pflag)
			ps = 3.5;
			f = file_C{idx};
			pdfprint(10+idx,['img',filesep,f(1:end-4),'-pattern.pdf'],ps);
			pdfprint(20+idx,['img',filesep,f(1:end-4),'-periodogram-2d.pdf'],ps);
			pdfprint(30+idx,['img',filesep,f(1:end-4),'-correlogram-2d.pdf'],ps);
		end % if pflag

		end % for idx
		disp(tab);
end % function

