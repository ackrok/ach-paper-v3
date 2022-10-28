load('C:\Users\Anya\Desktop\FP_LOCAL\rev_A\baseline\cannulaToneFull')

%%
% for x = 1:length(base)
%     ii = find(strcmp({c3.rec},base(x).rec));
%     c3(ii).nbFPbase = base(x).nbFP; 
%     c3(ii).onRestbase = base(x).onRest;  
%     c3(ii).offRestbase = base(x).offRest;
% end

%% BASELINE + CANNULA MODE
lblInf = fliplr(unique({c3.rx}));
out = cell(length(lblInf),1);
winInf = [0 10];
for x = 1:length(c3)
    Fs = c3(x).Fs;
    demod = c3(x).nbFP{1}; % Extract demod (demodulated, non-baselined) photometry signal for cannula recording
    idxImm = extractEventST([1:length(demod)]', c3(x).onRest, c3(x).offRest,1); % Sample index during immobility
    idxImm(idxImm < winInf(1)*60*Fs | idxImm > winInf(2)*60*Fs) = []; % Remove sample indices that are outside of infusion window
    demodImm = demod(idxImm); % Restrict to immobility + infusion window
    if ~isempty(c3(x).nbFPbase)
        demodBase = c3(x).nbFPbase{1};
        demodBase = demodBase(extractEventST([1:length(demodBase)]', c3(x).onRestbase, c3(x).offRestbase, 1));  % Extract demod vector for baseline recording, restrict to immobility
    else; demodBase = nan;
    end
    tmp = [mode(demodImm), mode(demodBase)];
    out{strcmp(lblInf, c3(x).rx)} = [out{strcmp(lblInf, c3(x).rx)}; tmp];
end
frac = nan(max(cellfun(@length, out)),length(out));
for x = 1:length(out)
    tmp = out{x}(:,1)./out{x}(:,2);
    frac(1:length(tmp),x) = tmp;
end
fprintf('Done extracting percent baseline.\n');

%% FRACTION of BASELINE
figure; hold on
plot(frac', '.k', 'MarkerSize', 20);
errorbar(0.25+[1:length(out)],nanmean(frac,1),SEM(frac,1),'.b', 'MarkerSize', 20);
xlim([0.5 0.5+length(out)]); xticks([1:length(out)]); xticklabels(lblInf);
ylabel('% of baseline'); % ylim([0 1]); yticks([0:0.25:1]);
% p = []; for x = 1:4; [~,p(x)] = ttest(b(:,1),b(:,x+1)); end
% title(sprintf('IMM. a/n %1.4f, a/scop %1.4f',p(2),p(4)));

