%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Janurary 3, 15:18:05 2018
%
% @Matlab Code: Mai Weixiong
%                Tang Liang
% @Email: maiweixiong@gmail.com
%         liang.tang@hotmail.com
%
% @AFD Author: Qian Tao
% @Home page: http://www.fst.umac.mo/en/staff/fsttq.html#afd
%
% Open two parameters: distn,contn for get_roots_new function
% 2018-06-13 09:53
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,rem1,S1,a,amp,pha,blaschke_z] = ...
    Unwinding_Blaschke(f,n,t,dist,cont,newtonOrder,tol,distn,contn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input : 'f' is real or analytic signal;
%         'n' is the maximal steps of the decomposition,the default n is 50;
%         't' is the interval [0,2*pi] which has the same size as f;
%         'dist' is the step distance and the default is 0.02;
%         'cont' is the contribution rate and the default is 0.95;
%         'newtonOrder' Newton-Cotes method,default value is 6
%
% Author: Mai Weixiong (maiweixiong@gmail.com)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
    if nargin == 1
        n = 50;
        t = linspace(0,2*pi,length(f));
        dist = 0.02;
        cont = 0.95;
        newtonOrder = 6;
        tol = 0.008;
        distn = 0.05;
        contn = 0.97;        
    end
    
    if nargin == 2
        t = linspace(0,2*pi,length(f));
        dist = 0.02;
        cont = 0.95;
        newtonOrder = 6;
        tol = 0.008;
        distn = 0.05;
        contn = 0.97;           
    end
    
    if nargin == 3
        dist = 0.02;
        cont = 0.95;
        newtonOrder = 6; 
        tol = 0.008;
        distn = 0.05;
        contn = 0.97;           
    end
    
    if nargin == 4
        cont = 0.95;
        newtonOrder = 6; 
        tol = 0.008;
    end    
    
    if nargin == 5
        newtonOrder = 6;  
        tol = 0.008;
        distn = 0.05;
        contn = 0.97;           
    end    
    
    if nargin == 6
        tol = 0.008;
        distn = 0.05;
        contn = 0.97;           
    end     
    
    if nargin == 7
        distn = 0.05;
        contn = 0.97;           
    end      
    
    if nargin == 8
        contn = 0.97;           
    end        
    
    if nargin > 9
        error('Unwinding: too many arguments...');
    end   

    [~,m] = size(f);
    G = zeros(n,m);
    a = zeros(n,1);
    SO = G;
    G(1,:) = f;
    x = exp(1i*t);
    S1 = zeros(n,1);
    C = Unit_Disk(dist,cont);
    [N,~] = size(C);
    f2 = Intg(f,f,newtonOrder);
    Weight = weight(length(f),newtonOrder);
    Base = zeros(length(C),m);
    B = zeros(n,length(t));
    pha = zeros(n,length(t));
    amp = zeros(n,length(t));
    blaschke_z = ones(n,100);

    for k = 1:N
        Base(k,:) = e_a(C(k),x);
    end
    r = get_roots_new(f,distn,contn,newtonOrder,tol);
    B(1,:) = blaschke1(r,f);
    if(isempty(r))
        pha(1,:) = 0;
        amp(1,:) = 0;
    else
        [pha(1,:),amp(1,:)] = disp_uw(r,B(1,:));
        blaschke_z(1,1:length(r)) = r;
    end

    SO(1,:) = f./B(1,:);
    in_prod = 1;
    a(1) = 0;
    S1(1) = Intg(SO(1,:),1,newtonOrder);
    in_prod = in_prod.*B(1,:);
    F(1,:) = S1(1).*Bn(a(1),exp(t.*1i)).*in_prod;
    fn = F(1,:);
    rem1 = abs(Intg(fn-f,fn-f,newtonOrder)/f2);
    expTemp = exp(t.*1i);
    for j=2:n
        G(j,:) = ((SO(j-1,:)-S1(j-1).*e_a(a(j-1),expTemp)).* ...
                        (1-conj(a(j-1)).*expTemp))./(expTemp-a(j-1));
        r = get_roots_new(G(j,:),distn,contn,newtonOrder,tol);
        B(j,:) = blaschke1(r,G(j,:));
        if(isempty(r))
            pha(j,:) = 0;
            amp(j,:) = 0;
        else
            [pha(j,:),amp(j,:)] = disp_uw(r,B(j,:));
            blaschke_z(j,1:length(r)) = r;
        end
        
        SO(j,:) = G(j,:)./B(j,:);
        S = conj(Base*(SO(j,:)'.*Weight));
        [~,I] = max(abs(S));
        S1(j) = S(I);
        a(j) = C(I);
        in_prod = in_prod.*B(j,:);
        F(j,:) = S1(j).*Bn(a(1:j),exp(t.*1i)).*in_prod;
        fn = fn + F(j,:);
        rem1 = abs(Intg(fn-f,fn-f,newtonOrder)/f2);
    end
end



function [ret] = Unit_Disk(dist,cont)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete the unit disk;
% input:
%           dist is the step distance and the default is 0.02;
%           cont is the contribution rate and the default is 0.95;
% Output:  
%           is an array;
% Author:   Mai Weixiong (maiweixiong@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin == 0
        dist = 0.02;
        cont = 0.95;
    end

    if nargin == 1
        cont = 0.95;
    end
    
    if nargin > 2
        error('Unit_Disk: too many arguments...');
    end  

    t = -1:dist:1;
    [~,n] = size(t);
    real = repmat(t',n,1);
    image = repmat(t',1,n);
    image = reshape(image',n^2,1);
    ret1 = real + 1i.*image;
    ret2 = ret1.*(abs(ret1) < cont);
    ret2 = ret2(abs(ret2) > 0);
    ret2 = sort(ret2); % Sorting does not have much effect on the result
    ret = [ret2;0];
end

function ret = Intg(f,g,newtonOrder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete the unit disk;
% input:
%           'f' is an array;
%           'g' is an array.
%           'newtonOrder'  Newton-Cotes method,default value is 6
%
% Author:   Mai Weixiong (maiweixiong@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin > 3
        error('intg: too many arguments...');
    end
    
    if nargin == 2
        newtonOrder = 6;
    end
    
    Weigth = weight(length(f),newtonOrder);
    ret = f*(g'.*Weigth);
end

function y = weight(n,newtonOrder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the integral by Newton-Cotes method
% input:
%           'n' is the length
%           'newtonOrder'  Newton-Cotes method,default value is 6
%
% Author:   Mai Weixiong (maiweixiong@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = zeros(n,1);
    if nargin == 1
        % Here we use the 6 Order
        newtonOrder = 6;              
    end
    
    % Newton-Cotes method
    Newton = zeros(newtonOrder+1,newtonOrder); 
    Newton(1:2,1) = [1/2;1/2];
    Newton(1:3,2) = [1/6;4/6;1/6];
    Newton(1:4,3) = [1/8;3/8;3/8;1/8];
    Newton(1:5,4) = [7/90;16/45;2/15;16/45;7/90];
    Newton(1:6,5) = [19/288;25/96;25/144;25/144;25/96;19/288];
    Newton(1:7,6) = [41/840;9/35;9/280;34/105;9/280;9/35;41/840];

    k = floor((n-1)/newtonOrder);
    if k > 0
        iter = 1:1:k;
        nonNewton = nonzeros(Newton(:,newtonOrder));
        exNewton = repmat(nonNewton(1:end-1),k,1);
        exNewton = [exNewton;nonNewton(end)];
        exIter = iter*newtonOrder + 1;
        exNewton(exIter) = exNewton(exIter) + nonNewton(end);
        exNewton(end) = exNewton(end) - nonNewton(end);
        y(1:k*newtonOrder+1) = y(1:k*newtonOrder+1) + exNewton;
    end
    
    y = y*newtonOrder/(n-1);
    nleft = n-k*newtonOrder-1;

    % for the left NLEFT points, use the corresponding Newton-Cotes method.
    if nleft > 0
        y((n-nleft):n) = y((n-nleft):n) + ...
            nonzeros(Newton(:,nleft))*nleft/(n-1);
    end
end

function [ ret ] = e_a(a,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the integral by Newton-Cotes method
% input:
%           'a' is a point in the unit disk;
%           'z' is a complex variable.
%
% Author:   Mai Weixiong (maiweixiong@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ret = ((1-abs(a)^2)^(0.5))./(1-conj(a).*z);
end

function [pha,amplitude] = disp_uw(a,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input : 'a' is a point in the unit disk;
%         'f' is real or analytic signal;
%
% Author:   Mai Weixiong (maiweixiong@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = linspace(0,2*pi,length(f));
    phase_d = 0;
    amplitude = abs(f);
    for j = 1:length(a)
        phase_d = (1-abs(a(j))^2)./(1-2.*abs(a(j)).*cos(t-angle(a(j)))+ ...
            abs(a(j))^2)+phase_d;
    end
    Mphase_d = max(max(phase_d));
    phase_d = phase_d/Mphase_d;
    pha = phase_d*Mphase_d;
end

function ret = Bn(a,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To calculate the Blaschke product with a;
% input : 'a' is a sequence of points in the unit disk of complex plane;
%         'z' is a complex variable.
%
% Author:   Mai Weixiong (maiweixiong@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n = length(a);
    if(n == 1)
        ret = ((sqrt(1-(abs(a(1))^2))))./(1-conj(a(1)).*z);
    else
        ret = ((1-abs(a(n))^2)^(0.5)./(1-conj(a(n)).*z));
        for j = 1:n-1
            ret = ret.*((z-a(j))./(1-conj(a(j)).*z));
        end
    end
end

function ret = get_roots_new(f,dist,cont,newtonOrder,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input : 'f' is real or analytic signal;
%         'dist' is the step distance and the default is 0.02;
%         'cont' is the contribution rate and the default is 0.95;
%         'newtonOrder' Newton-Cotes method,default value is 6
%
% Author:   Mai Weixiong (maiweixiong@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin > 5
        error('intg: too many arguments...');
    end
    
    if nargin == 4
        tol = 0.008;
    end
    
    if nargin == 3
        newtonOrder = 6;
        tol = 0.008;
    end
    
    if nargin == 2
        cont = 0.95;
        newtonOrder = 6;
        tol = 0.008;
    end
    
    if nargin == 1
        dist = 0.02;
        cont = 0.95;
        newtonOrder = 6;
        tol = 0.008;
    end    
    
    t = linspace(0,2*pi,length(f));
    x = exp(1i*t);
    C = Unit_Disk(dist,cont);
    N = length(C);
    Base = zeros(N,length(t));
    Weight = weight(length(f),newtonOrder);
    h = t(2)-t(1);
    f_der = zeros(size(f));
    
    for j = 3:length(f)-2
        f_der(j) = (f(j-2)-8*f(j-1)+8*f(j+1)-f(j+2))/(12*h);
    end
    
    f_der = -1i*f_der;
    num_zero = abs((sum(f_der(3:end-2)./f(3:end-2))/length(f(3:end-2))));
    num_zero = round(num_zero);
    r = [];
    
    for k = 1:N
        Base(k,:) = 1./(1-conj(C(k)).*x);
    end

    for j = 1:num_zero
        S = conj(Base*(f'.*Weight));
        [M,I] = min(abs(S));
        r(j) = C(I);
        f = f.*(1-conj(r(j))*x)./(x-r(j));
        if(M > tol)
            ret = r(1:j-1);
            break;
        end
        ret = r(1:j);
    end
    if(isempty(r))
        ret = [];
    end
end

function B=blaschke1(a,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input : 'f' is real or analytic signal;
%
% Author:   Mai Weixiong (maiweixiong@gmail.com)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = linspace(0,2*pi,length(f));
    x = exp(1i*t);
    B = ones(size(t));
    if(isempty(a))
        B = ones(size(f));
    else
        a1 = a(abs(a)~=0);
        for k = 1:length(a1)
            B = B.*(x-a1(k))./(1-conj(a1(k))*x).*(-conj(a1(k))/abs(a1(k)));
        end
        n = length(a)-length(a1);
        B = B.*x.^n;
    end
end

