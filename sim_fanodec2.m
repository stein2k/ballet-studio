function y = sim_fanodec()
    
    EsNo = 20;
    noiseVar = 10.^(-EsNo./10);

    hCRCGen = comm.CRCGenerator('Polynomial', [8 7 6 4 2 0]);
    hConEnc = comm.ConvolutionalEncoder('TrellisStructure', ...
        poly2trellis(3,[7 5]), 'TerminationMethod', 'Terminated');
    hMod = comm.BPSKModulator;
    hChan = comm.AWGNChannel('NoiseMethod', ...
        'Signal to noise ratio (Es/No)',...
        'EsNo', EsNo);
    hDemod = comm.BPSKDemodulator('DecisionMethod', ...
        'Approximate log-likelihood ratio', 'Variance', noiseVar);
    hError = comm.ErrorRate('ComputationDelay',3,'ReceiveDelay', 34);
    
    for counter = 1:20
        data = randi([0 1],30,1);
        encodedData = step(hConEnc, step(hCRCGen,data));
        modSignal = step(hMod, encodedData);
        receivedSignal = step(hChan, modSignal);
        demodSignal = step(hDemod, receivedSignal);
        receivedBits = fanodec2(demodSignal, [1 1 1 0 1 0 1 0 1], ...
            3, [7 5]);
        % receivedBits = step(hDec, demodSignal);
        errors = step(hError, data, receivedBits);
    end
    

end