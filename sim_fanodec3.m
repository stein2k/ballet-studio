function y = sim_fanodec3()
    
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
    
    EbNodB = 0:9;
    BitErrorRate = zeros(length(EbNodB), 2);
    BitErrorRate2 = zeros(length(EbNodB), 2);
    BitErrorRateAPP = zeros(length(EbNodB), 2);
    
    for i = 1:length(codes)
        
        Code = codes(i);
        Rate = 0.5 * (InformationBits/(InformationBits+ ...
            length(Code.CRCPolynomial)));
        
        for j = 1:length(EbNodB)
            
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
            hDecAPP = comm.APPDecoder('TrellisStructure', ...
                poly2trellis(Code.ConstraintLength, Code.CodeGenerator), ...
                'TerminationMethod', 'Terminated');
            hError = comm.ErrorRate('ComputationDelay',3,'ReceiveDelay', 34);

            num_errors = 0;
            num_errors2 = 0;
            num_errors_app = 0;
            num_bits = 0;
            
            while num_bits < InformationBits*1000 || num_errors < 100
                data = zeros(InformationBits, 1); %    randi([0 1],InformationBits,1);
                encodedData = step(hConEnc, step(hCRCGen,data));
                modSignal = step(hMod, encodedData);
                receivedSignal = step(hChan, modSignal);
                demodSignal = step(hDemod, receivedSignal);
                receivedBits = fanodec2(demodSignal, Code.CRCPolynomial, ...
                    Code.ConstraintLength, Code.CodeGenerator);
                receivedBits2 = (step(hDecAPP, zeros(InformationBits+ ...
                    length(Code.CRCPolynomial)+Code.ConstraintLength-2, ...
                    1), -demodSignal) > 0);
                receivedBits3 = fanodec(demodSignal, ...
                    Code.ConstraintLength, Code.CodeGenerator);
                num_errors = num_errors + sum(data~=receivedBits(1:InformationBits));
                num_errors2 = num_errors2 + sum(data~=receivedBits3(1:InformationBits));
                num_errors_app = num_errors_app + sum(data~=receivedBits2(1:InformationBits));
                num_bits = num_bits + InformationBits;
            end
            
            BitErrorRate(j, i) = num_errors/num_bits;
            BitErrorRate2(j, i) = num_errors2/num_bits;
            BitErrorRateAPP(j, i) = num_errors_app/num_bits;
            
            display([EbNodB(j) BitErrorRate(j, i) BitErrorRate2(j, i) BitErrorRateAPP(j,i)])
            
        end
        
    end
    
    save('fano_results.mat', 'BitErrorRate', 'BitErrorRate2', 'BitErrorRateAPP', 'EbNodB')

end