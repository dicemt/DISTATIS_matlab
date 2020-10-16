function varargout = distatis(data,nPCs,group3,groupcolor3,groupsym3,groupsiz3,grouplabel3,compcoding,compcolor,compsym,compsiz,complabel)
% Perform DISTATIS on a 3D distance (dissimilarity) matrix
% FORMAT [eigval,eigvector,fscore,eigval3,eigvector3,fscore3] = distatis(data)
% INPUTS:
% data           - Distance (dissimilarity) matrix [n x n x m]
%
% OPTIONAL INPUTS (for plots):
% nPCs           - The number of PCs to plot (default: 2)
% group3         - Group coding for dim 3 (eg, [1 1 2 2] % 2 groups)
% groupcolor3    - Group color for dim 3 (eg, [1 0 0; 0 0 1] % red, blue)
% groupsym3      - Group symbol for dim 3 (eg, 'x.' % x symbol, filled circle)
% groupsiz3      - Group size for dim 3 (eg, [4 6] % 2 groups)
% grouplabel3    - Group label for dim 3 (eg, {'pixels','measures','ratings','pairwise'})
% compcoding     - Compromise coding for dim 1/2 (eg, [1 1 1 2 2 2] % 2 groups)
% compcolor      - Compromise color for dim 1/2 (eg, [0 1 0; 1 0 0] % green, red)
% compsym        - Compromise symbol for dim 1/2 (eg, 'x.' % x symbol, filled circle)
% compsiz        - Compromise size for dim 1/2 (eg, [4 8] % 2 groups)
% complabel      - Compromise label for dim 1/2 (eg, {'M1','M2','M3','F1','F2','F3'})
% 
% OUTPUTS:
% eigval         - Eigenvalues for compromise
% eigvector      - Eigenvector for compromise
% fscore         - Factor score for compromise
% eigval3        - Eigenvalues for dim3 data
% eigvector3     - Eigenvector for dim3 data
% fscore3        - Factor score for dim3 data
% 
% REFERENCE:
% Abdi H, Toole AJO, Valentin D, Edelman B (2005)
% DISTATIS: The analysis of multiple distance matrices.
% IEEE CVPR, 42-47.
%
% Abdi H, Valentin D, Chollet S, Chrea C (2007)
% Analyzing assessors and products in sorting tasks: DISTATIS, theory and applications.
% Food Qual Prefer 18:627-640.
%
% REQUIREMENTS:
% Statistics and Machine Learning Toolbox
% plot_arrow.m
% https://www.mathworks.com/matlabcentral/fileexchange/3345-plot-arrowhead
%__________________________________________________________________________
% Copyright (C) 2020 Daisuke MATSUYOSHI
% $Id: distatis.m 0006 2020-10-16Z $

%% Data check
[nx,ny,nz] = size(data);
if nx ~= ny
    error('Invalid dimension (input matrix)!');
end

%% Settings
if nargin < 2 || isempty(nPCs)
    nPCs = 2;
end

if nargin < 3 || isempty(group3)
    group3 = 1:nz;
end

if nargin < 4 || isempty(groupcolor3)
    groupcolor3 = [];
end

if nargin < 5 || isempty(groupsym3)
    groupsym3 = [];
end

if nargin < 6
    groupsiz3 = [];
end

if nargin < 7 || isempty(grouplabel3)
    grouplabel3 = cellstr(num2str(transpose(1:nz)));
end

if nargin < 8 || isempty(compcoding)
    compcoding = 1:nx;
end

if nargin < 9 || isempty(compcolor)
    compcolor = [];
end

if nargin < 10 || isempty(compsym)
    compsym = [];
end

if nargin < 11
    compsiz = [];
end

if nargin < 12 || isempty(complabel)
    complabel = cellstr(num2str(transpose(1:nx)));
end



%% Initialize
matI = eye(nx);
m = ones(nx,1)/nx;
M = matI - ones(nx,1) * m';

%% Calc
% Cross product
Snorm = zeros(size(data));
for i=1:nz
    S = (-1/2) * M * data(:,:,i) * M';
    eigSt = eig(S);
    eigS = sort(eigSt,'descend');
    Snorm(:,:,i) = S/eigS(1);
end
X = reshape(Snorm,[],size(Snorm,3)); % complete data, normalized cross-product vectors

% Compromise matrix
A = X' * X; % scalar product matrix
N = diag(diag(A));
n = zeros(size(X,2),1);
for i=1:size(X,2)
    n(i) = norm(X(:,i),2); % L2
end
C = A ./ (n(:)*n(:).'); % Cosine matrix % N^(-1/2)*A*N^(-1/2)
Rv = A ./ sqrt(diag(A)*diag(A).'); % Congruence coefficients

% PCA
[V,D] = eig(C);
[d,idx] = sort(diag(D),'descend');
Ds = D(idx,idx);
Vs = V(:,idx);
G = Vs * diag(d.^(1/2));
varargout{4} = d;
varargout{5} = Vs;
varargout{6} = G;

% Optimal weights
P = Vs;
alpha = P ./ pinv((ones(length(P),1)/length(P))'*P)';
Splus = reshape(X * alpha(:,1),nx,nx);
[V,D] = eig(Splus);
[dSplus,idx] = sort(diag(D),'descend');
eigSplusD = D(idx,idx);
eigSplusV = V(:,idx);
Fscore = eigSplusV * diag(dSplus.^(1/2));
varargout{1} = dSplus;
varargout{2} = eigSplusV;
varargout{3} = Fscore;

%% Plots
% PCA Plot1 (PC1 & PC2) | Dim 3
PCx = 1; PCy = 2;
figure
hold on
gscatter(G(:,1),G(:,2),group3,groupcolor3,groupsym3,groupsiz3,'off') % weighted with eigenvalue
text(G(:,1), G(:,2), grouplabel3);
plot([0,0],ylim,':k','HandleVisibility','off')
plot([0,0],ylim,':k','HandleVisibility','off') % draw twice to plot at full length
plot(xlim,[0,0],':k','HandleVisibility','off')
xlabel(sprintf('PC%d (%.1f%%)',PCx,d(1)/sum(d)*100))
ylabel(sprintf('PC%d (%.1f%%)',PCy,d(2)/sum(d)*100))
hold off

% PCA Plot2 (PC1 & PC2) | Compromise Dim 1/2
figure
hold on
gscatter(eigSplusV(:,PCx),eigSplusV(:,PCy),compcoding,compcolor,compsym,compsiz,'off')
text(eigSplusV(:,PCx), eigSplusV(:,PCy), complabel);
plot([0,0],ylim,':k','HandleVisibility','off')
plot([0,0],ylim,':k','HandleVisibility','off') % draw twice to plot at full length
plot(xlim,[0,0],':k','HandleVisibility','off')
xlabel(sprintf('PC%d (%.1f%%)',PCx,dSplus(PCx)/sum(dSplus)*100))
ylabel(sprintf('PC%d (%.1f%%)',PCy,dSplus(PCy)/sum(dSplus)*100))
hold off

% PCA Plot3 (PC1 & PC2) | Compromise Dim 1/2 & 3
Fs = zeros(size(data));
for i=1:size(Snorm,3)
    Fs(:,:,i) = Snorm(:,:,i) * eigSplusV * diag(dSplus.^(-1/2));
end

if isempty(groupcolor3)
    cols = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 0 0 0];
    if nz < 15
        cols = [cols; lines];
    else
        cols = colorcube;
    end
else
    Gidx = unique(group3);
    nGroup = numel(unique(group3));
    nPerGroup = zeros(nGroup,1);
    for i=1:nGroup
        nPerGroup(i) = sum(group3==Gidx(i));
    end
    cols = zeros(nz,3);
    for i=1:nGroup
        cols((i-1)*nPerGroup(i)+1:i*nPerGroup(i),:) = repmat(groupcolor3(i,:),nPerGroup(i),1);
    end
end

figure
hold on
gscatter(eigSplusV(:,PCx),eigSplusV(:,PCy),compcoding,compcolor,compsym,compsiz,'off')
text(eigSplusV(:,PCx), eigSplusV(:,PCy), complabel);
plot(xlim,[0,0],':k','HandleVisibility','off')
for j=1:size(Fs,3)
    for i=1:length(eigSplusV)
        x0 = eigSplusV(i,PCx);
        y0 = eigSplusV(i,PCy);
        x1 = Fs(i,PCx,j);
        y1 = Fs(i,PCy,j);
        plot_arrow(x0, y0, x1, y1,'color',cols(j,:),'facecolor',cols(j,:));
    end
end
plot([0,0],ylim,':k','HandleVisibility','off')
plot([0,0],ylim,':k','HandleVisibility','off') % draw twice to plot at full length
xlabel(sprintf('PC%d (%.1f%%)',PCx,dSplus(PCx)/sum(dSplus)*100))
ylabel(sprintf('PC%d (%.1f%%)',PCy,dSplus(PCy)/sum(dSplus)*100))
hold off

if nPCs == 3
    % PCA Plot | Compromise Dim 1/2
    PCx = 1; PCy = 3;
    figure
    hold on
    gscatter(eigSplusV(:,PCx),eigSplusV(:,PCy),compcoding,compcolor,compsym,compsiz,'off')
    text(eigSplusV(:,PCx), eigSplusV(:,PCy), complabel);
    plot([0,0],ylim,':k','HandleVisibility','off')
    plot([0,0],ylim,':k','HandleVisibility','off') % draw twice to plot at full length
    plot(xlim,[0,0],':k','HandleVisibility','off')
    xlabel(sprintf('PC%d (%.1f%%)',PCx,dSplus(PCx)/sum(dSplus)*100))
    ylabel(sprintf('PC%d (%.1f%%)',PCy,dSplus(PCy)/sum(dSplus)*100))
    hold off

    % PCA Plot | Compromise Dim 1/2 & 3
    figure
    hold on
    gscatter(eigSplusV(:,PCx),eigSplusV(:,PCy),compcoding,compcolor,compsym,compsiz,'off')
    text(eigSplusV(:,PCx), eigSplusV(:,PCy), complabel);
    plot(xlim,[0,0],':k','HandleVisibility','off')
    for j=1:size(Fs,3)
        for i=1:length(eigSplusV)
            x0 = eigSplusV(i,PCx);
            y0 = eigSplusV(i,PCy);
            x1 = Fs(i,PCx,j);
            y1 = Fs(i,PCy,j);
            plot_arrow(x0, y0, x1, y1,'color',cols(j,:),'facecolor',cols(j,:));
        end
    end
    plot([0,0],ylim,':k','HandleVisibility','off')
    plot([0,0],ylim,':k','HandleVisibility','off') % draw twice to plot at full length
    xlabel(sprintf('PC%d (%.1f%%)',PCx,dSplus(PCx)/sum(dSplus)*100))
    ylabel(sprintf('PC%d (%.1f%%)',PCy,dSplus(PCy)/sum(dSplus)*100))
    hold off

    % PCA Plot | Compromise Dim 1/2
    PCx = 2; PCy = 3;
    figure
    hold on
    gscatter(eigSplusV(:,PCx),eigSplusV(:,PCy),compcoding,compcolor,compsym,compsiz,'off')
    text(eigSplusV(:,PCx), eigSplusV(:,PCy), complabel);
    plot([0,0],ylim,':k','HandleVisibility','off')
    plot([0,0],ylim,':k','HandleVisibility','off') % draw twice to plot at full length
    plot(xlim,[0,0],':k','HandleVisibility','off')
    xlabel(sprintf('PC%d (%.1f%%)',PCx,dSplus(PCx)/sum(dSplus)*100))
    ylabel(sprintf('PC%d (%.1f%%)',PCy,dSplus(PCy)/sum(dSplus)*100))
    hold off

    % PCA Plot | Compromise Dim 1/2 & 3
    figure
    hold on
    gscatter(eigSplusV(:,PCx),eigSplusV(:,PCy),compcoding,compcolor,compsym,compsiz,'off')
    text(eigSplusV(:,PCx), eigSplusV(:,PCy), complabel);
    plot(xlim,[0,0],':k','HandleVisibility','off')
    for j=1:size(Fs,3)
        for i=1:length(eigSplusV)
            x0 = eigSplusV(i,PCx);
            y0 = eigSplusV(i,PCy);
            x1 = Fs(i,PCx,j);
            y1 = Fs(i,PCy,j);
            plot_arrow(x0, y0, x1, y1,'color',cols(j,:),'facecolor',cols(j,:));
        end
    end
    plot([0,0],ylim,':k','HandleVisibility','off')
    plot([0,0],ylim,':k','HandleVisibility','off') % draw twice to plot at full length
    xlabel(sprintf('PC%d (%.1f%%)',PCx,dSplus(PCx)/sum(dSplus)*100))
    ylabel(sprintf('PC%d (%.1f%%)',PCy,dSplus(PCy)/sum(dSplus)*100))
    hold off
elseif nPCs > 3
    nLoop = floor((nPCs-2)/2);
    for i = 1:nLoop
        % PCA Plot | Compromise Dim 1/2
        PCx = i*2+1; PCy = i*2+2;
        figure
        hold on
        gscatter(eigSplusV(:,PCx),eigSplusV(:,PCy),compcoding,compcolor,compsym,compsiz,'off')
        text(eigSplusV(:,PCx), eigSplusV(:,PCy), complabel);
        plot([0,0],ylim,':k','HandleVisibility','off')
        plot([0,0],ylim,':k','HandleVisibility','off') % draw twice to plot at full length
        plot(xlim,[0,0],':k','HandleVisibility','off')
        xlabel(sprintf('PC%d (%.1f%%)',PCx,dSplus(PCx)/sum(dSplus)*100))
        ylabel(sprintf('PC%d (%.1f%%)',PCy,dSplus(PCy)/sum(dSplus)*100))
        hold off

        % PCA Plot | Compromise Dim 1/2 & 3
        figure
        hold on
        gscatter(eigSplusV(:,PCx),eigSplusV(:,PCy),compcoding,compcolor,compsym,compsiz,'off')
        text(eigSplusV(:,PCx), eigSplusV(:,PCy), complabel);
        plot(xlim,[0,0],':k','HandleVisibility','off')
        for j=1:size(Fs,3)
            for i=1:length(eigSplusV)
                x0 = eigSplusV(i,PCx);
                y0 = eigSplusV(i,PCy);
                x1 = Fs(i,PCx,j);
                y1 = Fs(i,PCy,j);
                plot_arrow(x0, y0, x1, y1,'color',cols(j,:),'facecolor',cols(j,:));
            end
        end
        plot([0,0],ylim,':k','HandleVisibility','off')
        plot([0,0],ylim,':k','HandleVisibility','off') % draw twice to plot at full length
        xlabel(sprintf('PC%d (%.1f%%)',PCx,dSplus(PCx)/sum(dSplus)*100))
        ylabel(sprintf('PC%d (%.1f%%)',PCy,dSplus(PCy)/sum(dSplus)*100))
        hold off
    end
end
