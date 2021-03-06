function y = sim_fanodec()
    
    EsNo = 20;
    noiseVar = 10.^(-EsNo./10);

    % hCRCGen = comm.CRCGenerator('Polynomial', [8 7 6 4 2 0]);
    %hConEnc = comm.ConvolutionalEncoder('TrellisStructure', ...
    %    poly2trellis(3,[7 5]), 'TerminationMethod', 'Terminated');
    hConEnc = comm.ConvolutionalEncoder('TrellisStructure', ...
        poly2trellis(7,[171 133]), 'TerminationMethod', 'Terminated');
    hMod = comm.BPSKModulator;
    hChan = comm.AWGNChannel('NoiseMethod', ...
        'Signal to noise ratio (Es/No)',...
        'EsNo', EsNo);
    hDemod = comm.BPSKDemodulator('DecisionMethod', ...
        'Approximate log-likelihood ratio', 'Variance', noiseVar);
    hDec = comm.ViterbiDecoder('InputFormat','Hard');
    hError = comm.ErrorRate('ComputationDelay',3,'ReceiveDelay', 34);
    
%     g1 = mod(conv(de2bi(oct2dec(7),3,'left-msb'), [1 1 1 0 1 0 1 0 1]),...
%         2);
%     g2 = mod(conv(de2bi(oct2dec(5),3,'left-msb'), [1 1 1 0 1 0 1 0 1],...
%         2);
    
    for counter = 1:20
        data = randi([0 1],30,1);
%         encodedData = step(hConEnc, step(hCRCGen,data));
        encodedData = step(hConEnc, data);
        modSignal = step(hMod, encodedData);
        receivedSignal = step(hChan, modSignal);
        demodSignal = step(hDemod, receivedSignal);
        receivedBits = fanodec(demodSignal, 7, [171 133]);
        % receivedBits = step(hDec, demodSignal);
        errors = step(hError, data, receivedBits);
    end
    

end