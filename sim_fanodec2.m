function y = sim_fanodec()
    
    % define system properties
    InformationBits = 128;
    
    codes = repmat(struct('CRCPolynomial', [], 'ConstraintLength', 0, ...
        'CodeGenerator', []), 2, 1);
    
    codes(1).CRCPolynomial = [1 1 1 0 1 0 1 0 1];
    codes(1).ConstraintLength = 3;
    codes(1).CodeGenerator = [7 5];
    
    codes(2).CRCPolynomial = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1];
    codes(2).ConstraintLength = 9;
    codes(2).CodeGenerator = [561 753];
    
    EbNodB = 0:10;
    BitErrorRate = zeros(length(EbNodB), 2);
    
    for i = 1:length(codes)
        
        Code = codes(i);
        Rate = 0.5 * (InformationBits/(InformationBits+ ...
            length(Code.CRCPolynomial)));
        
        for j = 4:length(EbNodB)
            
            % define channel properties
            EbNo = 10^(EbNodB(j)/10.0);
            EcNodB = EbNodB(j) + 10*log(Rate)/log(10);
            EcNo = 10^(EcNodB/10);
            No = 1./EcNo;

            hCRCGen = comm.CRCGenerator('Polynomial', fliplr(find( ...
                fliplr(Code.CRCPolynomial)==1)-1));
            hConEnc = comm.ConvolutionalEncoder('TrellisStructure', ...
                poly2trellis(Code.ConstraintLength, Code.CodeGenerator), ...
                'TerminationMethod', 'Terminated');
            hMod = comm.BPSKModulator;
            hChan = comm.AWGNChannel('NoiseMethod', ...
                'Variance', 'Variance', No);
            hDemod = comm.BPSKDemodulator('DecisionMethod', ...
                'Approximate log-likelihood ratio', 'Variance', No);
            hError = comm.ErrorRate('ComputationDelay',3,'ReceiveDelay', 34);

            num_errors = 0;
            num_bits = 0;
            
            while num_bits < InformationBits*1000 || num_errors < 100
                data = randi([0 1],InformationBits,1);
                encodedData = step(hConEnc, step(hCRCGen,data));
                modSignal = step(hMod, encodedData);
                receivedSignal = step(hChan, modSignal);
                demodSignal = step(hDemod, receivedSignal);
                receivedBits = fanodec2(demodSignal, Code.CRCPolynomial, ...
                    Code.ConstraintLength, Code.CodeGenerator);
                num_errors = num_errors + sum(data~=receivedBits(1:InformationBits));
                num_bits = num_bits + InformationBits;
            end
            
            BitErrorRate(j, i) = num_errors/num_bits;
            
            display([EbNodB(j) num_errors/num_bits])
            
        end
        
    end
    
    save('fano_results.mat', 'BitErrorRate', 'EbNodB')

end