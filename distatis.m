function varargout = distatis(data,group3,groupcolor3,groupsym3)
% Perform DISTATIS on 3D data
% FORMAT [eigval,eigvector,fscore,eigval3,eigvector3,fscore3] = distatis(data)
% INPUTS:
% data           - 3D matrix [n x n x m]
%
% OPTIONAL INPUTS:
% group3         - Group coding for dim 3 (eg, [1 1 2 2] % 2 groups, dim3 is 4)
% groupcolor3    - Group color for dim 3 (eg, 'br' % blue, red)
% groupsym3      - Group symbol for dim 3 (eg, 'x.' % x point, filled circle)
% 
% OUTPUTS:
% eigval         - Eigenvalues for compromise
% eigvector      - Eigenvector for compromise
% fscore         - Factor score for compromise
% eigval3        - Eigenvalues for dim3 data
% eigvector3     - Eigenvector for dim3 data
% fscore3        - Factor score for dim3 data
% 
% REFERENCE
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
% $Id: distatis.m 0001 2020-06-18Z $

if nargin < 2
    group3 = 1:size(data,3);
end

if nargin < 3
    groupcolor3 = [];
end

if nargin < 4
    groupsym3 = [];
end

nrow = size(data,1);
matI = eye(nrow);
m = ones(nrow,1)/nrow;
M = matI - ones(nrow,1) * m';

%% Calc
% Cross product
Snorm = zeros(size(data));
for i=1:size(data,3)
    S = (-1/2) * M * data(:,:,i) * M';
    eigSt = eig(S);
    eigS = sort(eigSt,'descend');
    Snorm(:,:,i) = S/eigS(1);
end
X = reshape(Snorm,[],size(Snorm,3)); % complete data, normalized cross-product vectors

% Compromise matrix
A = X' * X; % scalar product matrix
N = diag(diag(A));
for i=1:size(X,2)
    n(i) = norm(X(:,i),2); % L2
end
C = A ./ (n(:)*n(:).'); % Cosine matrix
Rv = A ./ sqrt(diag(A)*diag(A).'); % congruence coefficients

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
Splus = reshape(X * alpha(:,1),nrow,nrow);
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
figure
hold on
gscatter(G(:,1),G(:,2),group3,groupcolor3,groupsym3,[],'off') % give each component the length corresponding to its eigenvalue
text(G(:,1), G(:,2), cellstr(num2str([1:length(G)]')));
plot([0,0],ylim,':k','HandleVisibility','off')
plot([0,0],ylim,':k','HandleVisibility','off') % draw twice to plot at a full length
plot(xlim,[0,0],':k','HandleVisibility','off')
xlabel(sprintf('PC1 (%.1f%%)',d(1)/sum(d)*100))
ylabel(sprintf('PC2 (%.1f%%)',d(2)/sum(d)*100))
hold off

% PCA Plot2 (PC1 & PC2) | Compromise Dim 1/2
figure
hold on
gscatter(eigSplusV(:,1),eigSplusV(:,2),1:length(eigSplusV),[],[],[],'off')
text(eigSplusV(:,1), eigSplusV(:,2), cellstr(num2str([1:length(eigSplusV)]')));
plot([0,0],ylim,':k','HandleVisibility','off')
plot([0,0],ylim,':k','HandleVisibility','off') % draw twice to plot at a full length
plot(xlim,[0,0],':k','HandleVisibility','off')
xlabel(sprintf('PC1 (%.1f%%)',dSplus(1)/sum(dSplus)*100))
ylabel(sprintf('PC2 (%.1f%%)',dSplus(2)/sum(dSplus)*100))
hold off

% PCA Plot3 (PC1 & PC2) | Compromise Dim 1/2 & 3
for i=1:size(Snorm,3)
    Fs(:,:,i) = Snorm(:,:,i) * eigSplusV * diag(dSplus.^(-1/2));
end
if size(data,3) > 7
    cols = colorcube;
else
    cols = lines;
end
figure
hold on
gscatter(eigSplusV(:,1),eigSplusV(:,2),1:length(eigSplusV),[],[],[],'off')
text(eigSplusV(:,1), eigSplusV(:,2), cellstr(num2str([1:length(eigSplusV)]')));
plot([0,0],ylim,':k','HandleVisibility','off')
plot([0,0],ylim,':k','HandleVisibility','off') % draw twice to plot at a full length
plot(xlim,[0,0],':k','HandleVisibility','off')
for j=1:size(Fs,3)
    for i=1:length(eigSplusV)
        x0 = eigSplusV(i,1);
        y0 = eigSplusV(i,2);
        x1 = Fs(i,1,j);
        y1 = Fs(i,2,j);
        plot_arrow(x0, y0, x1, y1,'color',cols(j,:),'facecolor',cols(j,:));
    end
end
xlabel(sprintf('PC1 (%.1f%%)',dSplus(1)/sum(dSplus)*100))
ylabel(sprintf('PC2 (%.1f%%)',dSplus(2)/sum(dSplus)*100))
hold off
