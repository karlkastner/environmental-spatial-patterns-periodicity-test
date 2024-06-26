% Mon 11 Dec 10:21:54 CET 2023
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
if (~exist('pflag','var'))
	pflag = 0;
end
fflag = pflag;

type_C = {'anisotropic','isotropic'};

if (~exist('r2','var'))
r2 = {};

for tdx=1:2

load(['mat/vegetation-patterns-',type_C{tdx},'-analyzed.mat']);

n=length(spa.sp_a);
f_C={'gamma','logn','periodic'};
for idx=1:n
for jdx=1:length(f);
	 try
		 r2{tdx}(idx,jdx)=spa.sp_a(idx).stat.fit.radial.(f_C{jdx}).stat.goodness.r2;
	 catch e
		e
	 	r2{tdx}(idx,jdx)=NaN;
	end
end % for jdx
end % for idx
end % for tdx
end % if ~exist

for tdx=1:2
q = quantile(r2{tdx},[0.05,0.5,0.95])
splitfigure([2,2],[1,tdx],fflag);
cla();
qt=q;
hold on;
errorbar(1:3,qt(2,:),qt(2,:)-qt(1,:),qt(3,:)-qt(2,:),'*','linewidth',1);
xlim([0.5,3.5]);
set(gca,'xtick',1:length(f_C),'xticklabel',f_C)
hline(0,'color',[0.5,0.5,0.5]);
ylabel('R^2');
ylim([-0.6,1]);
axis square

end % for tdx


if (pflag)
	ps = 3.5;
	pdfprint(11,'img/r2-density-fit-anisotropic.pdf',ps);	
	pdfprint(12,'img/r2-density-fit-isotropic.pdf',ps);	
end

