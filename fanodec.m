function y = fanodec(code, trellis)

    % compute convolutional code parameters
    k = log2(trellis.numInputSymbols);
    n = log2(trellis.numOutputSymbols);

    % compute number of nodes
    numel = length(code) / n;
    
    % resahpe received codeword vector
    code = reshape(code, [n numel]);
    
    % declare and initialize path through code tree
    path = repmat(struct('metric',0,'state',0,'theta',1),numel,1);
    node = 1;
    
    while true
        
        % get output bits for current node in code tree
        x = code(node, :);
        
        % check whether node information or tail
        if node <= ninfo
            
            % create empty column vector for all paths emenating from
            % current node
            metric = zeros(trellis.numInputSymbols, 1);

            % look forward to best node
            for i = 0:trellis.numInputSymbols-1

                % get coded bits from curent node given input bits
                y = trellis.outputs(path(node), i);

                % compute branch metric
                metric(i) = (-2.*de2bi(y, n)+1) * x;

            end

            % find n:th best branch emanating from current node
            [B,I] = sort(metric);
            branch_metric = B(path(node).theta);
            branch = I(path(node).theta);
            
            
        else
            
            % get coded bits from current node given tail bits
            y = trellis.outputs(path(node).state, 1);
            
            % compute branch metric
            branch_metric = (-2.*de2bi(y,n)+1) * x;
            branch = 1;
            
        end
        
        % check threshold condition
        if path(node).metric + branch_metric >= threshold
            
            % check if first time we have visited this node
            if path(node).metric < threshold + delta
                
                % tighten threshold
                while path(node).metric + branch_metric >= threshold
                    threshold = threshold + delta;
                end
                
            end
            
            % store path and branch metric
            path(node+1).metric = path(node).metric + branch_metric;
            path(node+1).state = tellis.nextStates(path(node).state, ...
                branch);
            path(node+1).theta = 1;
            
            % move forward
            node = node + 1;
            
            continue;
            
        end
            
        while true

            % check backward path metric
            if node == 1 || path(node-1).metric < threshold
            
                % cannot back up, relax threshold and look forward
                threshold = threshold - delta;
                path(node).theta = 1;
                break;
                
            else

                % backward path metric is greater than threshold
                % move backward
                node = node - 1;
                
                % check whether previous node is part of information bits 
                % or whether it is part of tail bits, cannot follow
                % 'next best branch' out of tail node
                if node <= ninformation

                    % check if current node was reached by following
                    % worst branch emanating from previous node
                    if ~(path(node).theta == trellis.numInputSymbols)
                        path(node).theta = path(node).theta + 1;
                        break;
                    end
                    
                end

            end

        end
        
    end

end