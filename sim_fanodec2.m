function y = sim_fanodec()
    
    % define system properties
    InformationBits = 100;
    Rate = 0.5 * (InformationBits/(InformationBits+8));
    
    % define channel properties
    EbNodB = 6.0;
    EbNo = 10^(EbNodB/10.0);
    EcNodB = EbNodB + 10*log(Rate)/log(10);
    EcNo = 10^(EcNodB/10);
    EcnO = Rate*EbNo;
    No = 1./EcNo;

    hCRCGen = comm.CRCGenerator('Polynomial', [8 7 6 4 2 0]);
    hConEnc = comm.ConvolutionalEncoder('TrellisStructure', ...
        poly2trellis(3,[7 5]), 'TerminationMethod', 'Terminated');
    hMod = comm.BPSKModulator;
    hChan = comm.AWGNChannel('NoiseMethod', ...
        'Variance', 'Variance', No);
    hDemod = comm.BPSKDemodulator('DecisionMethod', ...
        'Approximate log-likelihood ratio', 'Variance', No);
    hError = comm.ErrorRate('ComputationDelay',3,'ReceiveDelay', 34);
    
    for counter = 1:20
        data = zeros(InformationBits, 1); %randi([0 1],InformationBits,1);
        encodedData = step(hConEnc, step(hCRCGen,data));
        modSignal = step(hMod, encodedData);
        receivedSignal = step(hChan, modSignal);
        demodSignal = step(hDemod, receivedSignal);
        receivedBits = fanodec2(demodSignal, [1 1 1 0 1 0 1 0 1], ...
            3, [7 5]);
        % receivedBits = step(hDec, demodSignal);
%         errors = step(hError, data, receivedBits);
    end
    

end