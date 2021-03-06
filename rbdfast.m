function [SIc,SI] = rbdfast(X,Y,Index,M)

%RBD Matlab Code for Random Balance Designs
%For the estimation of first order indices
%
%   SIc = rbdfast(X,Y) estimation of first order indices according to the X
%   values
%
%   SIc = rbdfast([],Y,Index) estimation of first order indices according to
%   the permutation used to create the sampling design
%
%   SIc = rbdfast(X,Y,[],M) estimation of first order indices with a defined
%   number of harmonics
%
%   X = N-by-k matrix of model inputs
%   Y = N-by-l matrix of model output
%   SIc = k-by-l matrix of estimated first order sensitivity indices unbiased
%   SI = k-by-l matrix of estimated first order sensitivity indices biased
%
%   N = sample size = total number of model evaluations
%   k = number of model input
%   l = number of model output
%
%
%Author: S. Tarantola (JRC)
%Joint Research Centre All rights Reserved
%
%Update: M. Rabouille
%Add: Reordering Y according a random design X (EASI algorithm) from E Plischke.
%Add: Unbiased estimator from J-Y Tissot & C Prieur.
%Note: The estimate is less dependant on the M value which can be raised up to 10.
%
%References: 
%S. Tarantola, D. Gatelli and T. Mara (2006)
%Random Balance Designs for the Estimation of First Order 
%Global Sensitivity Indices, Reliability Engineering and System Safety, 91:6, 717-727
%
%Elmar Plischke (2010)
%An effective algorithm for computing global sensitivity indices (EASI)
%Reliability Engineering & System Safety, 95:4, 354-360. <10.1016/j.ress.2009.11.005>
%
%Jean-Yves Tissot, Clémentine Prieur (2012)
%Bias correction for the estimation of sensitivity indices based on random balance designs.
%Reliability Engineering and System Safety, Elsevier,  107, 205-213. <10.1016/j.ress.2012.06.010> <hal-00507526v2>


% Test function
if nargin==0
    rbdfast_test()
    return
end

% Number of harmonics considered for the Fast Fourier Transform
if nargin < 4
    M=10;
elseif ~isreal(M) || length(M)~=1 || M <= 0 || floor(M) ~= M
	error('M must be a positive integer.')
end

if isempty(X)
    if nargin < 3 || isempty(Index)
        error('A index of permutation must be defined.')
    end
    [N, k] = size(Index); 
    useindex = true;
else
    Index=[];
    [N, k] = size(X);
    useindex = false;
end

if size(Y,1)~=N
    error('Arguments dimensions are not consistent')
end
if N<2*(M+100)
    warning('RBD:lowSampleSize','Insufficient simulations for proper analysis');
end

SI =zeros(k,size(Y,2));
SIc=zeros(k,size(Y,2));
lamda=(2*M)/N;

for a=1:k	
    if useindex
        % ---- reordering of y wrt ith index
        [~,ind]=sort(Index(:,a),'descend');
        Yorg=Y(ind,:);
    else
        % ---- reordering of y wrt ith variable
        [~,ind]=sort(X(:,a));
        ind=ind([1:2:N N-mod(N,2):-2:1]);
        Yorg=Y(ind,:);
    end

    %-----calculation spe1 at integer frequency
    spectrum=(abs(fft(Yorg))).^2/N/(N-1);
    % Normalization by N-1 to match the definition of the unbiased variance
    % We thus have the same definition as var(Y)
    % var(Yorg) ~ sum(spectrum(2:end,:))
    % var(Yorg|Xi) ~ 2*sum(spectrum(2:M+1,:)) = sum(spectrum(2:M+1,:))+sum(spectrum(end:-1:(end+1-M),:))

    V=sum(spectrum(2:end,:));
    SI(a,:)=2*sum(spectrum(2:M+1,:))./V;
    
    SIc(a,:)=SI(a,:)-lamda/(1-lamda)*(1-SI(a,:));

end	
end	

function rbdfast_test()
% ISHIGAMI fonction
% Crestaux et al. (2007) and Marrel et al. (2009) use: a = 7 and b = 0.1.
% Sobol' & Levitan (1999) use: a = 7 and b = 0.05.
a = 7;	b = 0.05;

fonc = @(X) sin(X(:,1)) + a*(sin(X(:,2))).^2 + b*(X(:,3)).^4.*sin(X(:,1));
ninput = 3; % def: X=[x1,x2,x3] -> xi=U(-pi,pi)

E = a/2;
Vx1 = 1/2*(1+b*pi^4/5)^2;
Vx2 = a^2/8;
Vx13 = b^2*pi^8*8/225;
V = Vx1 + Vx2 + Vx13;
exact = [   Vx1/V,    0,      0;
            0,        Vx2/V,  0;
            Vx13/V,   0,      0];
exact= diag(exact);


rng shuffle


SIc =zeros(ninput,500);
SI =zeros(ninput,500);
warning('off','RBD:lowSampleSize')
for N=50:500
    X = -pi + 2*pi.*rand(N,ninput);
    [SIc(:,N),SI(:,N)] = rbdfast(X, fonc(X));
end
warning('on','RBD:lowSampleSize')

figure
plot(1:N,SIc,'b',1:N,SI,'r')
hold on
plot([1 N],[exact exact],'k')
hold off
title('Effect of biais')
ylabel('SI')
xlabel('Simulation Number')


SIc2 =zeros(ninput,500);
warning('off','RBD:lowSampleSize')
for N=50:500
    s0=(-pi:2*pi/N:pi)'; % N+1 values from -pi to pi
    Index=zeros(N,ninput);
    for z=1:ninput
        Index(:,z)=randperm(N); % 3 independent random indexes from 1 to N
    end
    s=s0(Index); % Assigning values to the index -> "random" values between [-pi, pi[
    X=.5+asin(sin(s))/pi; % Uniform sampling in [0, 1]
    X=-pi+2*pi.*X; % Rescaling the uniform sampling between [-pi, pi]
    SIc2(:,N) = rbdfast([], fonc(X), Index);
end
warning('on','RBD:lowSampleSize')

figure
plot(1:N,SIc,'b',1:N,SIc2,'r')
hold on
plot([1 N],[exact exact],'k')
hold off
title('Effect of sample organisation')
ylabel('SI')
xlabel('Simulation Number')


SIc =zeros(ninput,30);
SI =zeros(ninput,30);
X = -pi + 2*pi.*rand(500,ninput);
for M=1:30
    [SIc(:,M),SI(:,M)] = rbdfast(X, fonc(X),[],M);
end
figure
plot(1:M,SIc,'b',1:M,SI,'r')
hold on
plot([1 M],[exact exact],'k')
hold off
title('Effect of the M value')
ylabel('SI')
xlabel('M value')


end

