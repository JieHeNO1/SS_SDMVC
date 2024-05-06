function [FreqsIdx,F_vec] = MixFreqs(Freqs,num_Harm,Fs,totaltime,maxFreq)
    Freqs = sort(Freqs,'descend');
    Freqs = Freqs(:);
    num_Freqs = numel(Freqs);
    maxFreq = min(maxFreq,Fs/2);
    num_Harmonics = min(floor(maxFreq./Freqs),50);
    num_Harmonics = min(num_Harmonics(:),num_Harm(:));
    F_vec = Freqs(1).*(1:num_Harmonics(1)).';
    F_vec = F_vec(:);
    if num_Freqs>1
        for i=2:num_Freqs
            if Freqs(i-1)>3*Freqs(i)
                F_vec_add = Freqs(i).*(-num_Harmonics(i):num_Harmonics(i)).';
            else 
                F_vec_add = (Freqs(i-1)-Freqs(i)).*(-num_Harmonics(i):num_Harmonics(i)+1).';
            end
            F_vec_add = F_vec_add(:);
            [F_mat,F_mat_add] = ndgrid(F_vec,F_vec_add);
            F_vec = unique(F_mat(:)+F_mat_add(:));
            F_vec(F_vec<=0) = [];
            F_vec(F_vec>=floor(Fs/2)) = [];
        end
    end
    FreqsIdx = totaltime.*F_vec;
    Idx_remove = abs(FreqsIdx-round(FreqsIdx))>eps;
    FreqsIdx(Idx_remove) = [];
    F_vec(Idx_remove) = [];
    FreqsIdx = FreqsIdx+1;
end