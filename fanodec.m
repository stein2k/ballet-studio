function Y = fanodec(code, ConstraintLength, CodeGenerator)

    % declare, initialize fano sequential decoder
    threshold = 0;
    delta = 4;
    
    % calculate constant for Fano metric
    R = mean(abs(code));

    % compute convolutional code parameters
    nstates = 2^(sum(ConstraintLength)-size(ConstraintLength,2));
    numInputSymbols = 2^size(ConstraintLength,2);
    numOutputSymbols = 2^size(CodeGenerator,2);
    k = log2(numInputSymbols);
    n = log2(numOutputSymbols);

    % compute number of nodes
    numel = length(code) / n;
    ninfo = numel-ConstraintLength+1;
    
    % resahpe received codeword vector
    code = reshape(code, [n numel]);
    
    % declare and initialize path through code tree
    path = repmat(struct('metric',0,'state',0,'branch',0,'theta',1), ...
        numel+1,1);
    node = 1;
    
    while node <= numel
        
        % get output bits for current node in code tree
        x = code(:, node);
        
        % inplace poly2trellis
        g = oct2dec(CodeGenerator);
        outputs = zeros(numInputSymbols,1);
        nextStates = zeros(numInputSymbols,1);
        for i = 1:numInputSymbols
            register = path(node).state + bitshift(i-1, ...
                ConstraintLength(1)-1);
            outputs(i) = bi2de(reshape(mod(sum(de2bi( ...
                bitand(register,g)),2),2), [1 size(CodeGenerator,2)]), ...
                'left-msb');
            nextStates(i) = bitshift(register,-1);
        end
        
        % check whether node information or tail
        if node <= ninfo
            
            % create empty column vector for all paths emenating from
            % current node
            metric = zeros(numInputSymbols, 1);

            % look forward to best node
            for i = 1:numInputSymbols

                % get coded bits from curent node given input bits
                y = -2*de2bi(outputs(i), n, 'left-msb')+1;

                % compute branch metric
                metric(i) = y*x - R;

            end

            % find n:th best branch emanating from current node
            [B,I] = sort(metric, 'descend');
            branch_metric = B(path(node).theta);
            branch = I(path(node).theta);
            
            % update branch emanating from current node
            path(node).branch = branch;
            
        else
            
            % get coded bits from current node given tail bits
            y = -2*de2bi(outputs(1), n, 'left-msb') + 1;
            
            % compute branch metric
            branch_metric = y*x - R;
            branch = 1;
            
            % update branch emanating from current node
            path(node).branch = branch;
            
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
            path(node+1).state = nextStates(branch);
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
                if node <= ninfo

                    % check if current node was reached by following
                    % worst branch emanating from previous node
                    if ~(path(node).theta == numInputSymbols)
                        path(node).theta = path(node).theta + 1;
                        break;
                    end
                    
                end

            end

        end
        
    end
    
    % declare, initialize decoded bits array
    Y = zeros(k, ninfo);
    
    % walk decoded fano path, get decoded bit matrix
    for i = 1:ninfo
        y = de2bi(path(i).branch-1, k, 'left-msb');
        Y(:,i) = y(:);
    end
    
    % reshape decoded bit matrix to column vector
    Y = reshape(Y,k*ninfo,1); 

end