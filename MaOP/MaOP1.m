classdef MaOP1 < PROBLEM
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
                f   = zeros(1, nobj);
                g   = sum((x(nobj:nvar) - 0.5).^2 + (1 - cos(20*pi*(x(nobj:nvar)-0.5))))/nvar;  %    '+' before cos is changed to '-' 2018.5.8
                tmp = 1;
                for m = nobj:(-1):1
                    id = nobj - m + 1;
                    if m>1
                        f(m) = (1 + g)*(1 - tmp*(1 - x(id)));
                        tmp  = tmp*x(id);
                    else
                        f(m) = (1 + g)*(1 - tmp);
                    end
                    f(m) = (0.1 + 10*m)*f(m);
                end
                PopObj(c,:)=f;
            end
            %             nvar=obj.D;
            %             nobj=obj.M;
            %             x=PopDec;
            %             g=zeros(size(x,1),1);
            %             f   = zeros(size(x,1), nobj);
            %             for c=1:size(x,1)
            %                 g(c)   = sum((x(c,nobj:nvar) - 0.5).^2 + (1 - cos(20*pi*(x(c,nobj:nvar)-0.5))))/nvar;  %    '+' before cos is changed to '-' 2018.5.8
            %             end
            %             tmp = ones(size(x,1),1);
            % %             if obj.N~=size(x,1)
            % %                 y=0
            % %             end
            %             for c=1:size(x,1)
            %                 for m = nobj:(-1):1
            %                     id = nobj - m + 1;
            %                     if m>1
            %                         %                     yy=tmp.*(1 - x(:,id));
            %                         f(c,m) = (1 + g(c)).*(1 - tmp(c).*(1 - x(c,id)));
            %                         tmp(c)  = tmp(c).*x(c,id);
            %                     else
            %                         f(c,m) = (1 + g(c)).*(1 - tmp(c));
            %                     end
            %                     f(c,m) = (0.1 + 10*m).*f(c,m);
            %                 end
            %             end
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