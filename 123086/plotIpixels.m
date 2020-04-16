%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an implementation of CONFIGR,
% as described in CAS/CNS Technical Report TR-2007-016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programmed by Chaitanya Sai (August 2007)
%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The notation follows the article's notation,
% as does the headers for each step of the
% algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function ret_val = plotIpixels(Ipixels,Iinterpol,filled_in_figr_rects_store,fig_num_in)

ImActSize=size(Ipixels)/5;


 
 


rects_num_unique=union(filled_in_figr_rects_store,filled_in_figr_rects_store,'rows');
rects_num_unique=rects_num_unique(find(rects_num_unique(:,1)==fig_num_in),:);

Ihigh=zeros(size(Ipixels));

for i=1:size(rects_num_unique,1);
    tVect=rects_num_unique(i,:);
    Ihigh(tVect(2):tVect(3),tVect(4):tVect(5))=1;
end


Ireverted=imresize((Ipixels==1),ImActSize);
Iinterpol_reverted=imresize(Iinterpol,ImActSize);
Ifigure_reverted=imresize((Ipixels==.5),ImActSize);
Iground_reverted=imresize((Ipixels==.25),ImActSize);
Ihigh_reverted=imresize(Ihigh,ImActSize);

folderplusfilenameBW='C:\MATLAB6p5p1\work\CONFIGR_RandomDots_Audit_Aug16\Figure'
close all
SimpleRoadImageFunky(Ireverted+Ifigure_reverted*.5+Ihigh_reverted*.042);
axis equal
axis tight
set(get(gcf,'Children'),'xtick',[])
set(get(gcf,'Children'),'ytick',[])


if size(rects_num_unique,1)~=0
for j=1:size(rects_num_unique,1)
    text(round((rects_num_unique(j,1)+rects_num_unique(j,2))/2*5),round((rects_num_unique(j,3)+rects_num_unique(j,4))/2*5),['Figure filling In number: ',num2str(j)],'fontsize',10)
end

xlabel(['Total number of figure filling in iteration: ',num2str(size(rects_num_unique,1))],'fontsize',18);
else
xlabel(['No figure filling in'],'fontsize',18);
end

title(['Iteration Number: ', num2str(fig_num_in-1)],'fontsize',18)

exportfig(gcf,[folderplusfilenameBW,num2str(fig_num_in),'.jpg'],'format','jpeg','resolution',1000,'height',3,'color','cmyk');
close all
ret_val=''