% 2023-06-02 14:43:13.373394646 +0200
%
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
function sp = pattern_manmade_vs_natural_plot(pflag)
	if (nargin()<1)
		pflag = 0;
	end
	fflag = pflag;
	folder = 'input';
	file_C = {
		 'nature-anisotropic-mexico_-103.147245_27.656891.png'
		,'nature-isotropic-mexico_-107.158326_31.337108.png'
		,'plantation-striped-spain_-2.685893_37.565588.png'
		,'plantation-hexagonal-spain_-5.356567_37.401112.png'
	};
	l = [3,3,3,3];
	
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
			tab.wavelength(idx) =  1./sp(idx).stat.fc.radial.bar;
			tab.L_rel(idx)      =  sp(idx).L(1).*sp(idx).stat.fc.radial.bar;
		else
			tab.wavelength(idx) =  1./sp(idx).stat.fc.x.bar;
			tab.L_rel(idx)      = sp(idx).L(1).*sp(idx).stat.fc.x.bar;
		end
		tab.intS_hp_sig(idx) = sp(idx).stat.stati.intS_hp_sig;
	end % for idx
	
	% plot patterns
	for idx=1:length(sp)

		% plot pattern		
		splitfigure([2,2],[1,idx],fflag);
		cla()
		%sp(idx).b_square = imread([folder,filesep,file_C{idx}(1:end-4),'-contrast-stretch.png']);
		sp(idx).plot('b');
		if (idx~=1)
			xlim([0,10]);
			ylim([0,10]);
		end
		
		% plot periodogram
		splitfigure([2,2],[2,idx],fflag);
		cla
		sp(idx).plot('S.rot.hp');
		c = clim();
		clim(0.5*c);
		axis([-1,1,-1,1]*l(idx))

		if (pflag)
			ps = 3.5;
			f = file_C{idx};
			pdfprint(10+idx,['img',filesep,f(1:end-4),'-pattern.pdf'],ps);
			pdfprint(20+idx,['img',filesep,f(1:end-4),'-periodogram-2d.pdf'],ps);
		end % if pflag
	
		end % for idx
		disp(tab)
end % function

