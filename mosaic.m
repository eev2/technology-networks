% Copyright (c) 2014, Samuel Cheng
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the  nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


function []=mosaic(data,labels,colors)

if min(data(:)) < 0
    error('data has to be non-negative');
end

xs = sum(data);
if (any(xs == 0))
  error('Reguire positive totals in each column');
end
data = data./repmat(xs, [size(data,1),1]);

xs=sum(data);
xs=xs/sum(xs);
ys=data*diag(1./sum(data));

xsp = min(min(xs(xs>0)));
ysp = min(min(ys(ys>0)));

gap=min([xsp,ysp])/4;
gap=min([gap,0.01]);

xs=[0 cumsum(xs)];
ys=cumsum(ys);
ys=[zeros(1,size(ys,2)); ys];

if ~exist('colors','var')
    colors={'r','b','g','k','c','m','y'};
end

if ~exist('labels','var')
   labels={};
end

for id=1:length(xs)-1
   plot_rectangles(xs(id),xs(id+1),ys(:,id),gap,labels,colors);
end

end

function plot_rectangles(x1,x2,ys,gap,labels,colors)

if ~exist('colors','var')
    colors={'r','b','g','k','c','m','y'};
end

if ~exist('gap','var')
    gap=0.01;
end

g2=gap/2;
xtmp=[x1+g2, x2-g2, x2-g2, x1+g2, x1+g2];
hold on;
for id=1:length(ys)-1
   if (ys(id+1) > ys(id))
     ytmp=[ys(id)+g2, ys(id)+g2, ys(id+1)-g2, ys(id+1)-g2, ys(id)+g2];
     plot(xtmp,ytmp,'Color',colors{mod(id-1,length(colors))+1});

     if exist('labels','var')
       if length(labels)>=id
         h=text((x1+x2)/2,(ys(id)+ys(id+1))/2,labels{id});
         set(h,'horizontalalignment','center');
         set(h,'verticalalignment','middle');
      end
     end
   end
end
hold off;
end
