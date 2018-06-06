fdr = 0.05; % False discovery rate for multiple comparison

dataf = 'classes';  % Which data to use: 'classes' or 'subclasses'

yrs = 1980:2011;
nyr = numel(yrs);

folder = ['../data/',dataf];
obsfitness = dlmread([folder,'/rowSums_1980-2011_4.tsv'],'\t',0,1);
codes99 = textread([folder,'/ipc_class_description_4.tsv'],'%s');
i99 = cellfun('length',regexp(codes99,'[A-H]+99')) == 1; % Remove 99
codes = codes99(~i99);
nnd = numel(codes); % Number of labels
obsfitness(:,i99) = []; % Remove X99 codes
obsfitness = obsfitness';

maxeig = zeros(1,nyr);
ncore = zeros(1,nyr);
naftercore = zeros(1,nyr);
acssize = zeros(1,nyr);
nodecore = zeros(nnd,nyr);
fitness = zeros(nnd,nyr);
pfe = zeros(nnd,nyr);
sedge = zeros(nnd,nnd,nyr);


for i = 1:nyr
  [ maxeig(i), pfe(:,i), ncore(i), naftercore(i), acssize(i), ...
    nodecore(:,i), fitness(:,i), A, sedge(:,:,i) ] = ...
      getstatyr5( yrs(i), fdr, folder, codes );
end


% Export for drawing the network
csvwrite([dataf,'_nodecore.csv'], nodecore);
csvwrite([dataf,'_sedge.csv'], reshape(sedge,[nnd*nnd*nyr,1]));

ntypeprop = zeros(3,nyr); % Core, Periphery, Rest
for i = 1:nyr
  ntypeprop(1,i) = mean(nodecore(:,i) == 2);
  ntypeprop(2,i) = mean(nodecore(:,i) == 1);
  ntypeprop(3,i) = mean(nodecore(:,i) == 0);
end

avgofit = zeros(3,nyr); % Core, ACS, Rest
for i = 1:nyr
    avgofit(1,i) = mean(obsfitness(nodecore(:,i) == 2,i));
    avgofit(2,i) = mean(obsfitness(nodecore(:,i) >= 1,i));
    avgofit(3,i) = mean(obsfitness(nodecore(:,i) == 0,i));
end

totofit = zeros(3,nyr); % Core, ACS, Rest
for i = 1:nyr
    totofit(1,i) = sum(obsfitness(nodecore(:,i) == 2,i));
    totofit(2,i) = sum(obsfitness(nodecore(:,i) >= 1,i));
    totofit(3,i) = sum(obsfitness(nodecore(:,i) == 0,i));
end

shareofit = zeros(3,nyr); % Core, ACS, Rest
for i = 1:nyr
    shareofit(1,i) = sum(obsfitness(nodecore(:,i) == 2,i));
    shareofit(2,i) = sum(obsfitness(nodecore(:,i) == 1,i));
    shareofit(3,i) = sum(obsfitness(nodecore(:,i) == 0,i));
end
shareofit = shareofit./repmat(sum(shareofit),[3,1]);

%modfitness = fitness./(max(max(fitness))).*max(max(obsfitness));
modfitness = pfe;
avgmfit = zeros(3,nyr); % Core, ACS, Rest
for i = 1:nyr
    avgmfit(1,i) = mean(modfitness(nodecore(:,i) == 2,i));
    avgmfit(2,i) = mean(modfitness(nodecore(:,i) >= 1,i));
    avgmfit(3,i) = mean(modfitness(nodecore(:,i) == 0,i));
end


secs = zeros(nnd,1); % Sections
nsec = zeros(8,1); % Number of classes in each section
for i = 1:8
    ii = cellfun('length',regexp(codes,[char(64+i) '[0-9][0-9]'])) == 1;
    nsec(i) = sum(ii);
    secs(ii) = i;
end
nseccs = cumsum(nsec);
propnode = zeros(8,nyr,3);
for k = 1:nyr
    for j = 1:8
        propnode(j,k,1) = mean(nodecore(secs == j,k) == 2);
        propnode(j,k,2) = mean(nodecore(secs == j,k) == 1);
        propnode(j,k,3) = mean(nodecore(secs == j,k) == 0);
    end
end
totsec = zeros(8,nyr,3);
for k = 1:nyr
    for j = 1:8
        totsec(j,k,1) = sum(secs(nodecore(:,k) == 2) == j);
        totsec(j,k,2) = sum(secs(nodecore(:,k) == 1) == j);
        totsec(j,k,3) = sum(secs(nodecore(:,k) == 0) == j);
    end
end

linksSec = zeros(2,nyr); % Within, Between
tmp = arrayfun(@ones, nsec, 'UniformOutput', false);
inSec = blkdiag(tmp{:});
for k = 1:nyr
  inACS = double(nodecore(:,k) >= 1); % In ACS
  sedgeIN = sedge(:,:,k).*inSec; % Links in the same section
  nSec = sum(sum(sedgeIN)); % Number of links in same section
  sedgeIN = sedgeIN .* (inACS*inACS'); % Only links part of ACS
  nACS = sum(sum(sedgeIN)); % Number of links in same section and ACS
  linksSec(1,k) = nACS ./ nSec;
  sedgeOT = sedge(:,:,k).*(1-inSec); % Links between sections
  nSec = sum(sum(sedgeOT)); % Number of links between sections
  sedgeOT = sedgeOT .* (inACS*inACS'); % Only links part of ACS
  nACS = sum(sum(sedgeOT)); % Number of links between sections and ACS
  linksSec(2,k) = nACS ./ nSec;
end

csvwrite([dataf,'_OBSCnt.csv'], totsec(:,:,1) + totsec(:,:,2));
csvwrite([dataf,'_nsec.csv'], nsec);

%% Plots
fg = figure('position',[10 10 550 500]);
subplot(2,2,1)
plot(yrs, maxeig)
%xlabel('Year')
title('Max eigenvalue')
axis tight

subplot(2,2,2)
plot(yrs, ncore)
%xlabel('Year')
title('Size of core')
axis tight

subplot(2,2,3)
plot(yrs, naftercore)
%xlabel('Year')
title('Size of periphery')
axis tight

subplot(2,2,4)
plot(yrs, acssize)
%xlabel('Year')
title('Size of ACS')
axis tight

saveas(fg,[dataf,'_acs.png'],'png')
close(fg)


%% Core indicator sorted
fg = figure('position',[10 10 900 500]);

myp = newplot;
myp.Layer = 'top';
myp.Box = 'on';
myp.YDir = 'reverse';
myp.TickLength = [0 0];

%ii=1:121;
[~,ii] = sortrows(nodecore);
imagesc(yrs,1:nnd,nodecore(ii,:));
colormap([1 1 1; 0 0 1; 1 0 0])
axis tight
title('Node type sorted: core (red); periphery (blue); weak (white)')

saveas(fg,[dataf,'_coreinds.png'],'png')
%close(fg)



%% Core indicator unsorted
%fg = figure('position',[10 10 900 500]);

[~,ii] = sortrows([secs nodecore]);
imagesc('XData',yrs,'YData',1:nnd,'CData',nodecore(ii,:));
colormap([1 1 1; 0 0 1; 1 0 0])
axis tight
title('Node type unsorted: core (red); periphery (blue); weak (white)')
hold on
for i = 1:7
    plot([yrs(1)-.5,yrs(end)+.5],[nseccs(i),nseccs(i)]+.5,'k','LineWidth',2)
end
% yticks(nseccs - .5*nsec)
% yticklabels({'A','B','C','D','E','F','G','H'})
set(gca,'YTick',nseccs - .5*nsec)
set(gca,'YTickLabel',{'A','B','C','D','E','F','G','H'})
hold off

saveas(fg,[dataf,'_coreindu.png'],'png')
close(fg)

%% Average Observed Fitness
fg = figure('position',[10 10 600 500]);

myp = plot(yrs,avgofit');
set(myp,{'color'},{[1 0 0]; [0 0 1]; [1 .8 0]})
legend({'Core', 'ACS', 'Rest'}, 'Location', 'northeast')
axis tight
title('Average observed fitness for each node type')

saveas(fg,[dataf,'_avgofit.png'],'png')
close(fg)

%% Total Observed Fitness
fg = figure('position',[10 10 600 500]);

myp = plot(yrs,totofit');
set(myp,{'color'},{[1 0 0]; [0 0 1]; [1 .8 0]})
legend({'Core', 'ACS', 'Rest'}, 'Location', 'east')
axis tight
title('Total observed fitness for each node type')

saveas(fg,[dataf,'_totofit.png'],'png')
close(fg)

%% Shared Observed Fitness
fg = figure('position',[10 10 600 500]);

myp = area(yrs, shareofit');
set(myp,{'FaceColor'},{[1 0 0]; [0 0 1]; [1 .8 0]})
legend({'Core', 'Periphery', 'Rest'}, 'Location', 'southeast')
% hold on
% plot(yrs,avgmfit','--')
% hold off
% legend({'Core O', 'ACS O', 'Rest 0', 'Core M', 'ACS M', 'Rest M'}, 'Location', 'northwest')
axis tight
ylim([0,1])
title('Proportion of observed fitness for each node type')

saveas(fg,[dataf,'_shareofit.png'],'png')
close(fg)

%% Number of nodes
fg = figure('position',[10 10 600 500]);

myp = area(yrs, ntypeprop');
set(myp,{'FaceColor'},{[1 0 0]; [0 0 1]; [1 .8 0]})
legend({'Core', 'Periphery', 'Rest'}, 'Location', 'southeast')
axis tight
ylim([0,1])
title('Proportion of nodes of each type')

saveas(fg,[dataf,'_ntypeprop.png'],'png')
close(fg)


%% Prop'n of types
fg = figure('position',[10 10 1800 500]);

subplot(1,3,1)
myp = plot(yrs,propnode(:,:,1)','LineWidth',2);
legend({'A','B','C','D','E','F','G','H'}, 'Location', 'southeast')
title('Proportion of nodes in core for each year')
set(myp,{'color'},num2cell(jet(8),2))

subplot(1,3,2)
myp = plot(yrs,propnode(:,:,1)' + propnode(:,:,2)','LineWidth',2);
legend({'A','B','C','D','E','F','G','H'}, 'Location', 'southeast')
title('Proportion of nodes in ACS for each year')
set(myp,{'color'},num2cell(jet(8),2))

subplot(1,3,3)
myp = plot(yrs,propnode(:,:,3)','LineWidth',2);
legend({'A','B','C','D','E','F','G','H'}, 'Location', 'southeast')
title('Proportion of nodes outside the ACS for each year')
set(myp,{'color'},num2cell(jet(8),2))

saveas(fg,[dataf,'_propnode.png'],'png')
close(fg)

%% Prop'n of types separately
fg = figure('position',[10 10 750 900]);
seclb = {'A','B','C','D','E','F','G','H'};
seccols = num2cell(jet(8),2);

for i = 1:8
subplot(4,2,i)
myp = plot(yrs,propnode(i,:,1) + propnode(i,:,2),'LineWidth',2);
title(seclb(i))
set(myp,{'color'},seccols(i))
xlim([1980,2011])
ylim([0,1])
end

saveas(fg,[dataf,'_propnodesep.png'],'png')
close(fg)

%% Prop'n of sections
fg = figure('position',[10 10 1800 500]);

subplot(1,3,1)
mosaic(propnode(:,:,1)./repmat(sum(propnode(:,:,1)),[8,1]), {}, ...
       num2cell(jet(8),2))
title('Proportion of each section in core')
set(gca, 'XTick', ((0:5:30)+.5)./32)
set(gca,'XTickLabel',num2cell(1980:5:2010))
box

subplot(1,3,2)
mosaic(sum(propnode(:,:,1:2),3)./repmat(sum(sum(propnode(:,:,1:2),3)),[8,1]),{}, ...
       num2cell(jet(8),2))
title('Proportion of each section in ACS')
set(gca, 'XTick', ((0:5:30)+.5)./32)
set(gca,'XTickLabel',num2cell(1980:5:2010))
box

subplot(1,3,3)
mosaic(propnode(:,:,3)./repmat(sum(propnode(:,:,3)),[8,1]), {}, ...
       num2cell(jet(8),2))
title('Proportion of each section not in ACS')
set(gca, 'XTick', ((0:5:30)+.5)./32)
set(gca,'XTickLabel',num2cell(1980:5:2010))
box

saveas(fg,[dataf,'_propsecsc.png'],'png')
close(fg)



fg = figure('position',[10 10 1800 500]);

subplot(1,3,1)
mosaic(propnode(:,:,1)./(nsec*sum(propnode(:,:,1))), {}, ...
       num2cell(jet(8),2))
title('Proportion of each section in core (Unscaled)')
set(gca, 'XTick', ((0:5:30)+.5)./32)
set(gca,'XTickLabel',num2cell(1980:5:2010))
box

subplot(1,3,2)
mosaic(sum(propnode(:,:,1:2),3)./(nsec*sum(sum(propnode(:,:,1:2),3))), {}, ...
       num2cell(jet(8),2))
title('Proportion of each section in ACS (Unscaled)')
set(gca, 'XTick', ((0:5:30)+.5)./32)
set(gca,'XTickLabel',num2cell(1980:5:2010))
box

subplot(1,3,3)
mosaic(propnode(:,:,3)./(nsec*sum(propnode(:,:,3))), {}, ...
       num2cell(jet(8),2))
title('Proportion of each section not in ACS (Unscaled)')
set(gca, 'XTick', ((0:5:30)+.5)./32)
set(gca,'XTickLabel',num2cell(1980:5:2010))
box

saveas(fg,[dataf,'_propsecun.png'],'png')
close(fg)


fg = figure('position',[10 10 600 500]);

myp = plot(yrs,linksSec');
legend({'Within sections', 'Between sections'}, 'Location', 'northwest')
axis tight
ylim([0,1])
title('Proportion of links part of the ACS')

saveas(fg,[dataf,'_linksSec.png'],'png')
close(fg)
