function y = sim_fanodec()

    trellis = poly2trellis(7, [171 133]);
    
    EsNo = 20;
    noiseVar = 10.^(-EsNo./10);

    hConEnc = comm.ConvolutionalEncoder('TerminationMethod','Terminated');
    hMod = comm.BPSKModulator;
    hChan = comm.AWGNChannel('NoiseMethod', ...
        'Signal to noise ratio (Es/No)',...
        'EsNo', EsNo);
    hDemod = comm.BPSKDemodulator('DecisionMethod', ...
        'Approximate log-likelihood ratio', 'Variance', noiseVar);
    hDec = comm.ViterbiDecoder('InputFormat','Hard');
    hError = comm.ErrorRate('ComputationDelay',3,'ReceiveDelay', 34);
    
    for counter = 1:20
        data = randi([0 1],30,1);
        encodedData = step(hConEnc, data);
        modSignal = step(hMod, encodedData);
        receivedSignal = step(hChan, modSignal);
        demodSignal = step(hDemod, receivedSignal);
        receivedBits = fanodec(demodSignal, 7, [171 133]);
        % receivedBits = step(hDec, demodSignal);
        errors = step(hError, data, receivedBits);
    end
    

end