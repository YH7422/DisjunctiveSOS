classdef MinHeap_BnB < handle

    properties
        k;                  % current number of elements
        n;                  % heap capacity
        x;                  % heap array
    end

    %
    % Public methods
    %
    methods (Access = public)
        %
        % Constructor
        %
        function this = MinHeap_BnB(n, x0)  
            this.n = n;
            this.x = repmat(Node(), n, 1);
            if ((nargin == 2) && ~isempty(x0))
                % Insert given elements
                k0 = numel(x0);
                if (k0 > n)
                    % Heap overflow
                    this.OverflowError();
                else
                    this.x(1:k0) = x0(:);
                    this.SetLength(k0);
                end
            else
                % Empty heap
                this.SetLength(0);
            end
        end
        
        %
        % Check for empty heap
        %
        function bool = IsEmpty(this)      
            if (this.k == 0)
                bool = true;
            else
                bool = false;
            end
        end
        
        %
        % Check for full heap
        %
        function bool = IsFull(this)            
            if (this.k == this.n)
                bool = true;
            else
                bool = false;
            end
        end
        
        %
        % Clear the heap
        %
        function Clear(this)     
            this.SetLength(0);
        end
        
        %
        % Insert key
        %
        function InsertKey(this,key)
            %------------------------ InsertKey ---------------------------
            % Syntax:       H.InsertKey(key);
            %               
            % Inputs:       key is a node
            %               
            % Description:  Inserts key into H
            %--------------------------------------------------------------
            
            this.SetLength(this.k + 1);
            this.x(this.k) = Node([], [], inf);
            this.DecreaseKey(this.k,key);
        end
        
        %
        % Return minimum element
        %
        function min = ReturnMin(this)
            %------------------------ ReturnMin ---------------------------
            % Syntax:       min = H.ReturnMin();
            %               
            % Outputs:      min is the minimum key in H
            %               
            % Description:  Returns the minimum key in H
            %--------------------------------------------------------------
            
            if (this.IsEmpty() == true)
                min = [];
            else
                min = this.x(1);
            end
        end
        
        %
        % Extract minimum element
        %
        function min = ExtractMin(this)
            %------------------------ ExtractMin --------------------------
            % Syntax:       min = H.ExtractMin();
            %               
            % Outputs:      min is the minimum key in H
            %               
            % Description:  Returns the minimum key in H and extracts it
            %               from the heap
            %--------------------------------------------------------------
            
            this.SetLength(this.k - 1);
            min = this.x(1);
            this.x(1) = this.x(this.k + 1);
            this.MinHeapify(1);
        end
    end
    
    %
    % Private methods
    %
    methods (Access = private)
        %
        % Decrease key at index i
        %
        function DecreaseKey(this,i,key)
            if (i > this.k)
                % Index overflow error
                MinHeap.IndexOverflowError();
            elseif (this.x(i) < key)
                % Decrease key error
                MinHeap.DecreaseKeyError();
            end
            this.x(i) = key;
            while ((i > 1) && (this.x(i) < this.x(this.parent(i))))
                this.Swap(i,this.parent(i));
                i = this.parent(i);
            end
        end
        
        %
        % Build the min heap
        %
        function BuildMinHeap(this)
            for i = floor(this.k / 2):-1:1
                this.MinHeapify(i);
            end
        end
        
        %
        % Maintain the min heap property at a given node
        %
        function MinHeapify(this,i,size)
            % Parse inputs
            if (nargin < 3)
                size = this.k;
            end
            
            ll = this.left(i);
            rr = this.right(i);
            if ((ll <= size) && (this.x(ll) < this.x(i)))
                smallest = ll;
            else
                smallest = i;
            end
            if ((rr <= size) && (this.x(rr) < this.x(smallest)))
                smallest = rr;
            end
            if (smallest ~= i)
                this.Swap(i,smallest);
                this.MinHeapify(smallest,size);
            end
        end
    end

    %
    % Protected methods
    %
    methods (Access = protected)
        %
        % Swap elements
        %
        function Swap(this,i,j)
            val = this.x(i);
            this.x(i) = this.x(j);
            this.x(j) = val;
        end
        
        %
        % Set length
        %
        function SetLength(this,k)
            if (k < 0)
                this.UnderflowError();
            elseif (k > this.n)
                this.OverflowError();
            end
            this.k = k;
        end
    end
    
    %
    % Private static methods
    %
    methods (Access = private, Static = true)
        %
        % Parent node
        %
        function p = parent(i)
            p = floor(i / 2);
        end
        
        %
        % Left child node
        %
        function l = left(i)
            l = 2 * i;
        end
        
        % Right child node
        function r = right(i)
            r = 2 * i + 1;
        end
        
        %
        % Overflow error
        %
        function OverflowError()
            error('Heap overflow');
        end
        
        %
        % Underflow error
        %
        function UnderflowError()
            error('Heap underflow');
        end
        
        %
        % Decrease key error
        %
        function DecreaseKeyError()
            error('You can only decrease keys in MinHeap');
        end
        
        %
        % Index overflow error
        %
        function IndexOverflowError()
            error('MinHeap index overflow');
        end
    end
end
