function output = makelr_0(nrows,ncols,k_use,epsilon_use,flag_nrm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [output] = makelr_0(nrows,ncols,k_use,epsilon_use,flag_nrm);
%
% This function generates a matrix A of size nrows-x-ncols which has ;
% numerical rank k_use with error epsilon_use. A is randomly generated to resemble a ;
% multivariate gaussian distribution. Specifically, the first k_use singular ;
% values of A are 1, whereas singular values k_use+1 and beyond equal epsilon_use, ;
% with the orientation of the principal components chosen uniformly at random. ;
% The inputs are nrows and ncols (integers) dictating the size of A, ;
%     the numerical rank k_use (an integer, default 1), ;
%     the error epsilon_use (a double betwen 0 and 1, default 0), ;
%     a normalization_flag flag_nrm (either 0 or 1, default 0). ;
% If flag_nrm==1, then then the values of A are normalized to match ;
%     variance with those in a uniformly gaussian random matrix. ;
% The output is A;
%
% test by running with no arguments:
% i.e., >> makelr_0();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% makes an epsilon_use-error k_use-rank block;
% tries to balance the variance ;
% test with: ;
%{

  N=512;nbins=128;dim=3;epsilon_use=0.1;
  nrows = N; ncols = 0.25*N;
  sigma = sqrt(dim*(1-epsilon_use.^2)/(nrows*ncols) + epsilon_use.^2/max(nrows,ncols));
  hbins = linspace(-6*sigma,+6*sigma,nbins);
  A = makelr_0(nrows,ncols,dim,epsilon_use);
  h = hist(A(:),hbins); h = h/sum(h)/mean(diff(hbins));
  he = 1/sqrt(2*pi*sigma.^2)*exp(-hbins.^2/2/sigma.^2);
  figure;cla;plot(hbins,h,'ro-',hbins,he,'k-'); xlim([min(hbins),max(hbins)]);
  disp(sprintf('std: %f vs %f',sigma,std(A(:))));
  N = 128;
  nrows = N; ncols = 0.5*N;
  A1 = makelr_0(nrows,ncols,1,0.1,1);
  A2 = makelr_0(nrows,ncols,3,0.1,1);
  A3 = makelr_0(nrows,ncols,1,1.0,1);
  A4 = makelr_0(nrows,ncols,3,1.0,1);
  A = [A1 , A2 , A3 , A4];
  figure;
  subplot(2,1,1);cla;imagesc(A);
  subplot(2,2,3);cla;plot(mean(A));title('mean');
  subplot(2,2,4);cla;plot(var(A));title('var');

 %}

if nargin<2;
disp(sprintf(' testing makelr_0: '));
disp(sprintf(' generating test data A of size 512x512'));
N=512;nbins=128;dim=3;epsilon_use=0.1;
nrows = N; ncols = 0.25*N;
sigma = sqrt(dim*(1-epsilon_use.^2)/(nrows*ncols) + epsilon_use.^2/max(nrows,ncols));
hbins = linspace(-6*sigma,+6*sigma,nbins);
A = makelr_0(nrows,ncols,dim,epsilon_use);
h = hist(A(:),hbins); h = h/sum(h)/mean(diff(hbins));
he = 1/sqrt(2*pi*sigma.^2)*exp(-hbins.^2/2/sigma.^2);
figure;cla;
subplot(2,2,1);
disp(sprintf(' plotting distribution of values of A.'));
plot(hbins,h,'ro-',hbins,he,'k-'); xlim([min(hbins),max(hbins)]);
title('distribution of A-values');
disp(sprintf(' sigma: %f vs distribution std %f',sigma,std(A(:))));
disp(sprintf(' generating test data A of size 256x256'));
N = 128;
nrows = N; ncols = 0.5*N;
A1 = makelr_0(nrows,ncols,1,0.1,1);
A2 = makelr_0(nrows,ncols,3,0.1,1);
A3 = makelr_0(nrows,ncols,1,1.0,1);
A4 = makelr_0(nrows,ncols,3,1.0,1);
A = [A1 , A2 , A3 , A4];
disp(sprintf(' plotting mean and variance for each column of A.'));
subplot(2,2,2);cla;imagesc(A);title('heatmap of A');
subplot(2,2,3);cla;plot(mean(A));title('mean');xlim([1,size(A,2)]);
subplot(2,2,4);cla;plot(var(A));title('var');xlim([1,size(A,2)]);
disp(sprintf(' printing figure (see test.jpg)'));
print('-djpeg','./test.jpg');
return;
end;%if nargin<2;

na=0;
if (nargin<1+na); nrows=[]; end; na=na+1;
if (nargin<1+na); ncols=[]; end; na=na+1;
if (nargin<1+na); k_use=[]; end; na=na+1;
if (nargin<1+na); epsilon_use=[]; end; na=na+1;
if (nargin<1+na); flag_nrm=[]; end; na=na+1;
if isempty(k_use); k_use=1; end;
if isempty(epsilon_use); epsilon_use=0; end;
if isempty(flag_nrm); flag_nrm=0; end;

output = randn(nrows,ncols);
[U,S,V] = svds(output,max(1024,4*k_use)); 
%S=diag(S); S(k_use+1:end) = epsilon_use*S(k_use+1:end); S=diag(S); 
S=diag(S); S(1:k_use)=1; S(k_use+1:end) = epsilon_use; S=diag(S); 
sigma=1;
if flag_nrm==1;
sigma = sqrt(k_use*(1-epsilon_use.^2)/(nrows*ncols) + epsilon_use.^2/max(nrows,ncols));
end;%if flag_nrm;
output = U*S*transpose(V)/sigma;

