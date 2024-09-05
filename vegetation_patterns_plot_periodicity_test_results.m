% Mon  8 May 13:25:55 CEST 2023
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
%% function vegetation_patterns_plot_periodicity_test_results(meta)
%%
%% plot and tabulate results of the periodicity test of the large global
%% dataset
%
function vegetation_patterns_plot_periodicity_test_results(meta)
	if (nargin()<1)
		meta = pattern_periodicity_test_metadata();
	end
	pflag  = meta.pflag;
	fflag  = pflag;
	% printscale
	ps     = 2;
	type_C = {'anisotropic','isotropic'};

	% significance level for accepting patterns to have periodic frequency components
	significance_level_a1 = meta.significance_level_a1;
	% significance_level_a1 = 0.01;
	% minimum number of patterns in a region for being displayed
	n_min = 50;
	% confindence levels for error bars in P	
	a_P   = 0.05;
	
	% plot extend
	ymax  = 0.2;

	nt = length(type_C);
	for tdx=1:nt
		% load analysis result
		type = type_C{tdx};
		clear spa
		iname = sprintf(meta.filename.patterns_analyzed_mat,type);
		%iname = sprintf(['mat',filesep,basename(meta.filename.patterns_analyzed),'-stat.mat'],type);
		load(iname,'spa');
		% prefetch
		stat.area_msk       = spa.area_msk;
		stat.p_periodic     = spa.p_periodic;
		stat.region_id      = spa.region_id;
		stat.intS_hp_sig    = spa.intS_hp_sig;
		stat.quality_check  = spa.quality_check;

		% number of regions
		nr = length(spa.region_C);
		if (1==tdx)
			area_m2   = zeros(nt,nr+1);
			p_med     = zeros(nt,nr+1);
			L_eff     = zeros(nt,nr+1);
			w_p_per   = zeros(nt,nr+1);
			N_per     = zeros(nt,nr+1);
			N         = zeros(nt,nr+1);
			N_qc      = zeros(nt,nr+1);
			ci        = zeros(nt,nr+1);
			intS_hp_sig = zeros(nt,nr+1);
		end

		% numbner of patterns
		n = spa.n;
		%length(stat.p_periodic);
	
		switch (type)
		case {'isotropic'}
			wavelength = spa.wavelength_r;
			L_eff_     = spa.L_eff_r;
		case {'anisotropic'}
			wavelength = spa.wavelength_x;
			L_eff_     = spa.L_eff_x;
		end % switch type

		% determine fraction of periodic patterns per region
		% process other last
		for rdx=[2:nr+1,1]
			switch (rdx)
			case {nr+1} % all regions, exclude patterns in other regions
				fdx = stat.region_id ~= 1; 
				%fdx = true(n,1);
			otherwise
				% select patterns in current region
				fdx = (stat.region_id == rdx);
			end
			% number of patterns in the region before quality check 
			N(tdx,rdx) = sum(fdx(1:n));
			
			fdx =   cvec(fdx) & stat.quality_check;
%		              & (rdx ~= 1); % exclude other regions

			% number of patterns after quality check
			N_qc(tdx,rdx) = sum(fdx);

			if (N_qc(tdx,rdx) < n_min)
				N_qc(tdx,rdx) = 0;
				% assign to other
	%			stat.region_id(fdx) == 0;
				continue;
			end
			area_m2(tdx,rdx) = nanmedian(stat.area_msk(fdx)); %arrayfun(@(x) x.stat.area_msk, stat_a(fdx)));
			p_med(tdx,rdx)    = median(stat.p_periodic(fdx));
			N_per(tdx,rdx)    = sum(stat.p_periodic(fdx)<=significance_level_a1);
			w_p_per(tdx,rdx)  = wmean(stat.area_msk(fdx),stat.p_periodic(fdx)<=significance_level_a1);
	 		fdx_ = cvec(stat.intS_hp_sig)>0;
			intS_hp_sig(tdx,rdx) = mean(stat.intS_hp_sig(fdx & fdx_)); 
			L_eff(tdx,rdx) = median(L_eff_(fdx)./wavelength(fdx));
		end % for rdx (each region) 
	end % for tdx (pattern type)
	
	% fractions of patterns with significant frequency components
	P_periodic = N_per./N_qc;
	if (0)
		% approximate confidence intervals
		s_pp  = sqrt(P_periodic.*(1-P_periodic))./sqrt(N_qc);
		s     = norminv(1-a_P/2);
		ci    = s*s_pp;
		ci(:,:,2) = ci;
	else
		% note, makedist does not set the covariance right, so pseudo resampling and fitdist
		pd = cell(size(P_periodic));
		for idx=1:size(P_periodic,1)
		for jdx=1:size(P_periodic,2)
			%p = p_per(idx,jdx);
			n_per = N_per(idx,jdx);
			n = N_qc(idx,jdx);
			x = zeros(n,1);
			x(1:n_per)=1;
			if (n_per > 0)
			pd{idx,jdx} = fitdist(x,'binomial');
			ci_ = paramci(pd{idx,jdx},'alpha',a_P); 
			ci(idx,jdx,:) = ci_(:,2)-P_periodic(idx,jdx);
			else
			ci(idx,jdx,1:2) = NaN;
			end
		end
		end
	end
	
	for tdx=1:length(type_C)
		sdx = find(N_qc(tdx,:)>0);
	
		splitfigure([2,2],[1,1+2*(tdx-1)],fflag);
		cla();
		plot(1:length(sdx),p_med(tdx,sdx),'*');
		%lab = {sp.region_C{:},'all'};
		lab = {spa.region_C{:},'.                             all'};
		set(gca,'xtick',1:length(sdx),'xticklabel',lab(sdx),'xticklabelrot',45)
		xlim([0.5,length(sdx)+0.75]);
		hline(significance_level_a1);
		ylabel('median($p$)','interpreter','latex');
	
		splitfigure([2,2],[1,2+2*(tdx-1)],fflag);
		cla();
		errorbar(1:length(sdx),P_periodic(tdx,sdx),ci(tdx,sdx,1),ci(tdx,sdx,2),'*','linewidth',1);
		set(gca,'xtick',1:length(sdx),'xticklabel',lab(sdx),'xticklabelrot',45)
		xlim([0.5,length(sdx)+0.75]);
		ylim([0,ymax]);
		hline(significance_level_a1);
		ylabel(sprintf('E$[p\\le %0.2g]$',significance_level_a1),'interpreter','latex');
		for idx=1:length(sdx)
			text(idx, P_periodic(tdx,sdx(idx)),['  ',num2str(N_qc(tdx,sdx(idx)))]); 
		end
		if (~pflag)
			title(type_C{tdx})
		end
	
		if (pflag)
			pdfprint(12+2*(tdx-1),['img/pattern-',type_C{tdx},'-p-periodic-threshold-',num2str(significance_level_a1),'.pdf'],ps);
		end
	end

	if (0)
	% number of pixels
	npxl = [];	
	npxl = 0.5*n;
	an = significance_level_a1./nt
	critical = 1/2*chi2inv(1 - an,2)./(nt);
	end	

	disp(N)
	disp(N_qc)
	disp(P_periodic)
	fprintf('                                                  aniso iso\n');
	fprintf('Number of patterns                                %4d %4d\n',N(:,end));
	fprintf('Number of patterns retained after quality check   %4d %4d\n',N_qc(:,end));
	fprintf('Fraction of patterns retained after quality check %3.1f %3.1f %%\n',N_qc(:,end)./N(:,end)*100);
	fprintf('Median p-value                                    %1.2f %1.2f\n',p_med(:,end));
	fprintf('Number of patterns with p < a1                    %4d %4d\n',N_per(:,end));
	fprintf('Expected number of patterns with p < a1           %3.1f %3.1f\n',N_qc(:,end)*significance_level_a1);
	fprintf('Fraction of patterns with p < a1                  %1.2f %1.2f %%\n',N_per(:,end) ./ N_qc(:,end)*100)
	fprintf('Fraction of patterns with p < a1, area weighted   %1.2f %1.2f %%\n',w_p_per(:,end)*100)
	fprintf('Expected fraction of patterns with p < a1         %3.2f %3.2f %%\n',significance_level_a1*100,significance_level_a1*100);
	fprintf('Median area per pattern                           %3.2g %3.2g km^2\n', area_m2(:,end)/1e6);
	%fprintf('Median number of grid cells per pattern       %g %g\n', round(area_m2(:,end)/dx^2));
	fprintf('Median effective pattern length L_eff                    %3.2g %3.2g\n',L_eff(:,end));
	% 	fprintf('Median required fraction of spectral energy for pattern to be accepted as periodic: %f %%\n',100*critical);
	fprintf('Mean fraction of spectral energy contained in siginificant components %3.2g%% %3.2g%%\n', 100*intS_hp_sig(:,end));
	
end % function

