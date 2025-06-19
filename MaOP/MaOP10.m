classdef MaOP10 < PROBLEM
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
                tmp = prod(sin(0.5*pi*x(1:(nobj-1))));
                for n=nobj:1:nvar
                    if mod(n,5)==0
                        g = g + (x(n) - tmp).^2;
                    else
                        g = g + (x(n) - 0.5).^2;
                    end;
                end
                
                g = g*100;   % The scale 100 is multiplied by g. 2018.5.8
                
                alpha = zeros(1, nobj);
                tau = sqrt(2)/2;
                
                alpha(1) = -(2*x(1)-1).^3 + 1;
                T = floor((nobj-1)/2);
                
                for i=1:T
                    z = 2*(2*x(i+1) - floor(2*x(i+1))) - 1;
                    if x(i+1)<0.5
                        p = 0.5 + x(1);
                    else
                        p = 1.5 - x(1);
                    end
                    alpha(2*i)    = x(1) + 2*x(i+1)*tau + tau*abs(z).^p;
                    alpha(2*i+1)  = x(1) - (2*x(i+1)-2)*tau + tau*abs(z).^p;
                end
                
                if mod(nobj,2)==0
                    alpha(nobj) = 1 - alpha(1);
                end
                
                f = (1 + g)*alpha;
                PopObj(c,:)=f;
            end
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