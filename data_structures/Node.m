classdef Node < handle
    properties 
        q;                   % coefficients of cuts
        s;                   % corresponding sos components
        lb;                  % lower bound of p*
        ub;                  % upper bound of p*
        x;                   % a point in the region
        region;              % corresponding region of the node
        vertices;            % corresponding vertices of the node
        delta;               % maximum distance between the vertices of the node
        cutting_edge;
    end
    
    methods
        
        % create a node
        function x = Node(q, sq, lb, ub, x0)
            if nargin > 5
                error('Too many input arguments.')
            end
            if nargin >= 5
                x.x = x0;
            else
                x.x = [];
            end
            if nargin >= 4
                x.ub = ub;
            else
                x.ub = inf;
            end
            if nargin >= 3
                x.lb = lb;
            else
                x.lb = -inf;
            end
            if nargin >= 2
                x.s = sq;
            else
                x.s = [];
            end
            if nargin >= 1
                x.q = q;
            else
                x.q = [];
            end
            x.region = [];
            x.vertices = [];
            x.delta = -inf;
            x.cutting_edge = [];
        end

        % compare the lower bounds of two nodes
        function flag = lt(x, y)
            flag = (x.lb < y.lb);
        end
    
    end
end