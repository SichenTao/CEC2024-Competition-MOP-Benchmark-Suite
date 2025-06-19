classdef MaOP6 < PROBLEM
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
                g = zeros(1, nobj);
                for m=1:nobj
                    if m<=3
                        g(m) =  max(0, 1.4*sin(4*x(1)*pi)) + sum(abs(x(3:nvar) - x(1)*x(2)).^2);
                    else
                        g(m) =  exp((x(m) - x(1)*x(2)).^2) - 1;
                    end;
                    g(m) = g(m)*10;
                end
                alpha1 = x(1)*x(2);
                alpha2 = (1-x(2))*x(1);
                alpha3 = (1-x(1));
                f(1) = (1 + g(1))*alpha1;
                f(2) = 2*(1 + g(2))*alpha2;
                f(3) = 6*(1 + g(3))*alpha3;
                for m = 4:nobj
                    f(m) = (1 + g(m))*(m*alpha1/nobj + (1-m/nobj)*alpha2 + sin(0.5*m*pi/nobj)*alpha3);  % n in g(n) -->g(m)
                end
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