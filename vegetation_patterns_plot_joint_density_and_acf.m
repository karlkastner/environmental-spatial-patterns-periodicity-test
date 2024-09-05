% Sun 26 Nov 21:38:23 CET 2023
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

meta = pattern_periodicity_test_metadata();
if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;
ps = 4;

field = meta.field;

q = 0.05;
L = 5;
xi = linspace(0,L,20*L)';
df = 1/40;
fi = linspace(0,10,10/df)';
Rcq = [];
Scq = [];
jc = fzero(@(x) besselj(1,2*pi*x),1);

type_C = {'anisotropic','isotropic'};
f_C = {'gamma','logn','normpdf_wrapped','periodic'};
f_legend_C = {'Gamma','Log-Normal','Normal','Periodic'};
%f_legend_C = {'gamma','lognormal','wrapped normal','periodic'};

if (~exist('loaded','var') || ~loaded)
loaded = true;

Rc_a = [];
Sc_a = [];
lc_a = [];
r2_C = {};
Ri_C = {};
Si_C = {}; 

for jdx=1:length(type_C)

filename = sprintf(meta.filename.patterns_analyzed_mat,type_C{jdx});
%filename = ['mat/vegetation-patterns-',type_C{jdx},'-analyzed.mat'];
load(filename);

n  = length(spa.sp_a);

Ri_C{jdx} = NaN(length(xi),n);
Si_C{jdx} = NaN(length(fi),n);

% max density
Sc_a = padd2d(Sc_a,[n,2],NaN);
% first max of acf
Rc_a = padd2d(Rc_a,[n,2],NaN);
% wavelength
lc_a = padd2d(lc_a,[n,2],NaN);
r2_C{jdx} = NaN(n,length(f_C));

lc__a = [];
% weigths
w = zeros(length(spa.sp_a),1);
for idx=1:round(n)
	%dx = spa.dx_sample(idx);
	switch (type_C{jdx})
	case {'anisotropic'}
	if (0 == spa.sp_a(idx).stat.isisotropic)
		x   = spa.sp_a(idx).x;
		x   = x-x(1);
		R   = spa.sp_a(idx).R.rot.x.(field);
		fc  = spa.sp_a(idx).stat.fc.x.(field);
		lc  = 1./fc;
		lc_a(idx,jdx) = lc;
		fdx = find(x > 0.5*lc & x < 1.5*lc);	
		fx  = spa.sp_a(idx).f.x;

		if (sum(fdx)>0)
		[Rc_a(idx,jdx),mdx] = max(R(fdx));
		lc_ = x(fdx(mdx));
		w(idx) = spa.sp_a(idx).stat.L_eff.x;
		Ri = interp1(x/lc_,R,xi,'linear');

		
		fdx = fx>=0;
		% C = cumsum([0;cvec(spa.sp_a(idx).S.rot.xp.(field)(fdx))]);
		C = cumsum([0;cvec(spa.sp_a(idx).S.rot.x.(field)(fdx))]);
		Ci = interp1(inner2outer(fx(fdx))/fc,C,inner2outer(fi),'linear');
		Ci = Ci/Ci(round(end/2));
		Si = diff(Ci);
		Si = Si/(sum(Si)*(fi(2)-fi(1)));
		%Si = interp1(fx,spa.sp_a(idx).S.rot.x.con/lc,fi,'linear');

		%R_(:,jdx) = R_(:,jdx) + Ri;
		Ri_C{jdx}(:,idx) = Ri;
		Si_C{jdx}(:,idx) = Si;
		%Sc_a(idx,jdx) = spa.sp_a(idx).stat.Sc.x.(field)*fc;
		Sc_a(idx,jdx) = spa.sp_a(idx).stat.Sc.xp.(field)*fc;
%		Sc(idx,1) = spa.sp_a(idx).stat.Sc.xp.con*fc;
		try
		for kdx=1:length(f_C)
			r2_C{jdx}(idx,kdx)=spa.sp_a(idx).stat.fit.x.(f_C{kdx}).stat.goodness.r2;
		end % for kdx
		catch
		end
		end 
	end
	case {'isotropic'}
	if (1 == spa.sp_a(idx).stat.isisotropic)
		fc  = spa.sp_a(idx).stat.fc.radial.(field);
		lc  = 1./fc;
		lc_a(idx,2) = lc;
		r   = spa.sp_a(idx).r;
		R    = spa.sp_a(idx).R.radial.(field);
		Jhat = sqrt(besselj(0,2*pi*r/lc).^2+besselj(1,2*pi*r/lc).^2);
		R    = R./Jhat;
		fdx = find(r > (jc-0.5)*lc & r < (jc+0.5)*lc);	
		if (lc<(spa.sp_a(idx).r(end)))
			[~,mdx] = max(R(fdx)); %.*sqrt(r(fdx))); % x/lc
			Rc_a(idx,jdx) = R(fdx(mdx));
			lc_ = r(fdx(mdx));
			lc__a(idx) = lc_;
			Ri = interp1(r/(lc_/jc),R,xi,'linear');
			%Ri = interp1(r/(lc),R,xi,'linear');
			S = spa.sp_a(idx).S.radial.(field);
			C = cumsum(spa.sp_a(idx).S.radial.(field));
			C = C/C(end);
			fr = spa.sp_a(idx).f.r;
			Ci = interp1(fr/fc,C,fi,'linear');
			Si = cdiff(Ci);
			Si = Si/(sum(Si)*(fi(2)-fi(1)));
			Ri_C{jdx}(:,idx) = Ri;
			Si_C{jdx}(:,idx) = Si;
			Sc_a(idx,jdx) = spa.sp_a(idx).stat.Sc.radial.(field)*fc;
			try
			for kdx=1:length(f_C)
				r2_C{jdx}(idx,kdx)=spa.sp_a(idx).stat.fit.radial.(f_C{kdx}).stat.goodness.r2;
			end % for kdx
			catch
			end
		end
	end
	end
	%end
end % for idx
end % for jdx
end

for jdx=1:length(type_C)
	Rcq(:,jdx) = quantile(Rc_a(:,jdx),[q,0.5,1-q]);
	Scq(:,jdx) = quantile(Sc_a(:,jdx),[q,0.5,1-q]);
	w=w/sum(w);
	fdx = (w>0);
	R_ = [nanmean(Ri_C{jdx},2)];
	% minimize the hellinger distance
	S_ = [nanmean(sqrt(Si_C{jdx}),2).^2];
	S_ = S_./sum(S_*(fi(2)-fi(1)));

	splitfigure([2,3],[1,(jdx-1)*3+1],fflag);
	cla;
	if (0) %2 == jdx)
	plot(xi,R_); %.*(pi*sqrt(xi))]);
	hold on;
	%Rcq_ = (pi)*Rcq;
	errorbar(1,Rcq_(2,jdx),Rcq_(2,jdx)-Rcq_(1,jdx),Rcq_(3,jdx)-Rcq_(2,jdx),'k*');
	else
	plot(xi,R_,'linewidth',1);
	hold on;
	if (1==jdx)
		errorbar(1,Rcq(2,jdx),Rcq(2,jdx)-Rcq(1,jdx),Rcq(3,jdx)-Rcq(2,jdx),'k*');
	else
		errorbar(jc,Rcq(2,jdx),Rcq(2,jdx)-Rcq(1,jdx),Rcq(3,jdx)-Rcq(2,jdx),'k*');
	end
	end

	if (1==jdx)
	xlabel('Distance $x/\lambda_c$','interpreter','latex');
	ylabel('Autocorrelation $R_x$','interpreter','latex');
	else
	xlabel('Distante $r/\lambda_c$','interpreter','latex');
	%ylabel('Autocorrelation $R_r/\hat J(0,2 \pi x/\lambda_c)$','interpreter','latex');
	ylabel('Autocorrelation $R_r/\hat J_0(k_c\,r)$','interpreter','latex');
	end
	ylim([-0.4,1]);
	xlim([0,3]);
	axis square

	splitfigure([2,3],[1,(jdx-1)*3+2],fflag);
	%if (1 == jdx) cla; end
	cla();
	plot(fi,S_,'linewidth',1);
	hold on
	errorbar(1,Scq(2,jdx),Scq(2,jdx)-Scq(1,jdx),Scq(3,jdx)-Scq(2,jdx),'k*');
	ylim([0,1.2]);
	xlim([0,3.5]);
	axis square
	if (1==jdx)
	xlabel('Wavenumber $k_x/k_c$','interpreter','latex');
	ylabel('Density $S_x^+/\lambda_c$','interpreter','latex');
	else
	xlabel('Wavenumber $k_r/k_c$','interpreter','latex');
	ylabel('Density $S_r/\lambda_c$','interpreter','latex');
	end

	r2_q = quantile(r2_C{jdx},[0.05,0.5,0.95])
	splitfigure([2,3],[1,(jdx-1)*3+3],fflag);
	cla();
	errorbar(1:length(f_C),r2_q(2,:),r2_q(2,:)-r2_q(1,:),r2_q(3,:)-r2_q(2,:),'*','linewidth',1);
	xlim([0.5,length(f_C)+0.5]);
	set(gca,'xtick',1:length(f_C),'xticklabel',f_legend_C)
	hline(0,'color',[0.5,0.5,0.5]);
	ylabel('Goodness of fit $R^2$','interpreter','latex');
	axis square
	ylim([-1,1]);


	figure(2);
	for idx=1:6;
	fdx = (spa.region_id==idx);
	%S_ = nanmean(Si_C{jdx}(:,fdx),2);
	S_ = [nanmean(sqrt(Si_C{jdx}(:,fdx)),2).^2];
	subplot(2,6,idx+(jdx-1)*6);
	cla();
	plot(fi,S_);
	title(sum(fdx));
	ylim([0,1]);
	Scq_ = cvec(quantile(Sc_a(fdx,jdx),[q,0.5,1-q]));
	hold on
	errorbar(1,Scq_(2,1),Scq_(2,1)-Scq_(1,1),Scq_(3,1)-Scq_(2,1),'k*');
	end
	 figure(3);
	 for idx=1:6;
	 fdx = (spa.region_id==idx);
	 R_ = nanmean(Ri_C{jdx}(:,fdx),2);
	 subplot(2,6,idx+(jdx-1)*6);
	 cla
	 plot(xi,R_);
	 title(sum(fdx));
	 Rcq_ = cvec(quantile(Rc_a(fdx),[0.05,0.5,0.95]));
	 hold on
%	 errorbar(1,Rcq_);
	if (jdx==1)
	 errorbar(1,Rcq_(2,1),Rcq_(2,1)-Rcq_(1,1),Rcq_(3,1)-Rcq_(2,1),'k*');
	else
	 errorbar(jc,Rcq_(2,1),Rcq_(2,1)-Rcq_(1,1),Rcq_(3,1)-Rcq_(2,1),'k*');
	end
	 end
	 ylim([-0.5,1]);

	% fit to average density
	if (0)
	splitfigure([2,3],[1,(jdx-1)*3+3],fflag);
	hold on;
	pdf_C = {@gampdf,@lognpdf,@normpdf_wrapped};
	for kdx=1:length(pdf_C);
		l = length(fi);
		switch (type_C{jdx})
		case {'anisotropic'}
			w = ones(l,1);
		case {'isotropic'}
			w = (0:l-1)';
		end
		% TODO weight with r
		[a,b,out] = fit_spectral_density(fi,S_,w,pdf_C{kdx},[1,1]);
		r2_(kdx,jdx) = out.goodness.r2;
	end
	plot(1:3,r2_(:,2),'ro');
	end
end
%
figure(5);
subplot(2,2,1)
	loglog(lc_a,Sc_a./lc_a,'.');

%subplot(2,2,3);
%errorbar(1:2,Rcq(2,:),Rcq(2,:)-Rcq(1,:),Rcq(3,:)-Rcq(2,:),'*');
%ylim([0,1]);
%for idx=1:2
%	%text(jdx,Rcq(2,jdx),['  ',num2str(nc(jdx))]);
%end

%	for idx=1:2
%	splitfigure([2,2],[6,idx],fflag);
%	q = quantile(r2{idx},[0.05,0.5,0.95]);
%	barplot(q);
%	ylabel('R^2');
%	set(gca,'xtick',1:3,'xticklabel',f_C);
%	end
%	splitfigure([2,2],[6,1],fflag);
%	q = quantile(r2{1},[0.05,0.5,0.95]);
%	barplot(q);
%	set(gca,'xtick',1:3,'xticklabel',f_C);
	

if (pflag)
	pdfprint(11,'img/anisotropic-patterns-autocorrelation.pdf',ps);
	pdfprint(12,'img/anisotropic-patterns-density.pdf',ps);
	pdfprint(13,'img/r2-density-fit-anisotropic.pdf',ps);	
	pdfprint(14,'img/isotropic-patterns-autocorrelation.pdf',ps);
	pdfprint(15,'img/isotropic-patterns-density.pdf',ps);
	pdfprint(16,'img/r2-density-fit-isotropic.pdf',ps);	
end

