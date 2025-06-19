classdef MaOP2 < PROBLEM
    % <multi/many> <real> <large/none> <expensive/none>
    % Benchmark MOP proposed by Deb, Thiele, Laumanns, and Zitzler
    
    %------------------------------- Reference --------------------------------
    % K. Deb, L. Thiele, M. Laumanns, and E. Zitzler, Scalable test problems
    % for evolutionary multiobjective optimization, Evolutionary multiobjective
    % Optimization. Theoretical Advances and Applications, 2005, 105-145.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    
    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M+4; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            nvar=obj.D;
            nobj=obj.M;
            x=PopDec;
            PopObj=zeros(size(x,1),nobj);
            
            for c=1:size(x,1)
                f = zeros(1, nobj);
                g = 0;
                tmp1 = prod(sin(0.5*pi*x(1:(nobj-1))));
                for n=nobj:1:nvar
                    if mod(n,5)==0
                        g = g + (x(n) - tmp1).^2;
                    else
                        g = g + (x(n) - 0.5).^2;
                    end;
                end
                g = g*200;    % The scale 200 is added to g.
                
                tmp2 = 1;
                for m = nobj:(-1):1
                    p = 2^(mod(m,2)+1);
                    if m==nobj
                        f(m) = (1 + g)*sin(0.5*x(1)*pi)^p;
                    elseif m<nobj&&m>=2
                        tmp2 = tmp2*cos(0.5*pi*x(nobj-m));
                        f(m) = (1 + g)*(tmp2*sin(0.5*pi*x(nobj-m+1)))^p;        % wrong pow function
                    else
                        f(m) = (1 + g)*(tmp2*cos(0.5*pi*x(nobj-1)))^p;          % wrong pow function
                    end
                end
                PopObj(c,:)=f;
            end
            %             nvar=obj.D;
            %             nobj=obj.M;
            %             x=PopDec;
            %             g=zeros(obj.N,1);
            %             f   = zeros(obj.N, nobj);
            %             %             tmp1 = ones(size(x,1),1);
            %
            %             for c=1:size(x,1)
            %                 tmp1(c) = prod(sin(0.5*pi*x(c,1:(nobj-1))));
            %             end
            %
            %             for c=1:size(x,1)
            %                 for n=nobj:1:nvar
            %                     if mod(n,5)==0
            %                         g(c) = g(c) + (x(c,n) - tmp1(c)).^2;
            %                     else
            %                         g(c) = g(c) + (x(c,n) - 0.5).^2;
            %                     end
            %                 end
            %             end
            %             g = g*200;    % The scale 200 is added to g.
            %
            %             %            tmp2 = 1;
            %             tmp2 = ones(size(x,1),1);
            %             for c=1:size(x,1)
            %                 for m = nobj:(-1):1
            %                     p = 2^(mod(m,2)+1);
            %                     if m==nobj
            %                         f(c,m) = (1 + g(c))*sin(0.5*x(c,1)*pi)^p;
            %                     elseif m<nobj&&m>=2
            %                         tmp2(c) = tmp2(c)*cos(0.5*pi*x(c,nobj-m));
            %                         f(c,m) = (1 + g(c))*(tmp2(c)*sin(0.5*pi*x(c,nobj-m+1)))^p;        % wrong pow function
            %                     else
            %                         f(c,m) = (1 + g(c))*(tmp2(c)*cos(0.5*pi*x(c,nobj-1)))^p;          % wrong pow function
            %                     end
            %                 end
            %             end
            %
            %             PopObj=f;
        end
        
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            name =([class(obj) '_F' num2str(obj.M) '.txt']);
            R=load(name);
            %             R = UniformPoint(N,obj.M)/2;
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            
            %            x =load(['MaOP' num2str(i) '_F' num2str(j) '.txt']);
            %               load(fullfile(fileparts(CallStack(1).file),'MaOP1_F3.txt'),'Dataset');
            %             if obj.M == 2
            %                 R = obj.GetOptimum(100);
            %             elseif obj.M == 3
            %                 a = linspace(0,1,10)';
            %                 R = {a*a'/2,a*(1-a')/2,(1-a)*ones(size(a'))/2};
            %             else
            R = [];
            %             end
        end
    end
end