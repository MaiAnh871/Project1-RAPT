doback=0;   % don't do backwards DP for now
fs = 8000;
% read in parameters
PP=voicebox;
f0min=50;            % Min F0 (Hz)                               []
f0max=500;            % Max F0 (Hz)                               []
tframe=0.01;          % frame size (s)                            []
tlpw=0.005;              % low pass filter window size (s)           []
tcorw=0.0075;            % correlation window size (s)               [0.0075]
candtr=0.3;          % minimum peak in NCCF                      []
lagwt=0.3;            % linear lag taper factor                   []
freqwt=0.02;          % cost factor for F0 change                 []
vtranc=0.005;          % fixed voice-state transition cost         []
vtrac=0.5;            % delta amplitude modulated transition cost []
vtrsc=0.5;            % delta spectrum modulated transition cost  [0.5]
vobias=0.0;          % bias to encourage voiced hypotheses       []
doublec=0.35;        % cost of exact doubling or halving         []
absnoise=0;      % absolute rms noise level                  []
relnoise=2.0;      % rms noise level relative to noise floor   []
signoise=0.001;      % ratio of peak signal rms to noise floor   []
ncands=20;          % max hypotheses at each frame              []
trms=0.03;              % window length for rms measurement         []
dtrms=0.02;            % window spacing for rms measurement        []
preemph=-7000;        % s-plane position of preemphasis zero      []
nfullag=7;        % number of full lags to try (must be odd)  []

% derived parameters (mostly dependent on sample rate fs)

krms=round(trms*fs);            % window length for rms measurement
kdrms=round(dtrms*fs);          % window spacing for rms measurement
rmswin=hanning(krms).^2;
kdsmp=round(0.25*fs/f0max);
hlpw=round(tlpw*fs/2);          % force window to be an odd length
blp=sinc((-hlpw:hlpw)/kdsmp).*hamming(2*hlpw+1).';
fsd=fs/kdsmp;
kframed=round(fsd*tframe);      % downsampled frame length
kframe=kframed*kdsmp;           % frame increment at full rate
rmsix=(1:krms)+floor((kdrms-kframe)/2); % rms index according to Talkin; better=(1:krms)+floor((kdrms-krms+1)/2)
minlag=ceil(fsd/f0max);
maxlag=round(fsd/f0min);        % use round() only because that is what Talkin does
kcorwd=round(fsd*tcorw);        % downsampled correlation window
kcorw=kcorwd*kdsmp;             % full rate correlation window
spoff=max(hlpw-floor(kdsmp/2),1+kdrms-rmsix(1)-kframe);  % offset for first speech frame at full rate
sfoff=spoff-hlpw+floor(kdsmp/2); % offset for downsampling filter
sfi=1:kcorwd;                   % initial decimated correlation window index array
sfhi=1:kcorw;                   % initial correlation window index array
sfj=1:kcorwd+maxlag;
sfmi=repmat((minlag:maxlag)',1,kcorwd)+repmat(sfi,maxlag-minlag+1,1);
lagoff=(minlag-1)*kdsmp;        % lag offset when converting to high sample rate
beta=lagwt*f0min/fs;            % bias towards low lags
log2=log(2);
lpcord=2+round(fs/1000);        % lpc order for itakura distance
hnfullag=floor(nfullag/2);
jumprat=exp((doublec+log2)/2);  % lag ratio at which octave jump cost is lowest
ssq=s.^2;
csssq=cumsum(ssq);
sqrt(min(csssq(kcorw+1:end)-csssq(1:end-kcorw))/kcorw);
afact=max([absnoise^2,max(ssq)*signoise^2,min(csssq(kcorw+1:end)-csssq(1:end-kcorw))*(relnoise/kcorw)^2])^2*kcorw^2;

% downsample signal to approx 2 kHz to speed up autocorrelation calculation
% kdsmp is the downsample factor

sf=filter(blp/sum(blp),1,s(sfoff+1:end));
sp=filter([1 exp(preemph/fs)],1,s); % preemphasised speech for LPC calculation
sf(1:length(blp)-1)=[];         % remove startup transient
sf=sf(1:kdsmp:end);             % downsample to =~2kHz
nsf=length(sf);                 % length of downsampled speech
ns=length(s);                   % length of full rate speech

% Calculate the frame limit to ensure we don't run off the end of the speech or decimated speech:
%   (a) For decimated autocorrelation when calculating sff():  (nframe-1)*kframed+kcorwd+maxlag <= nsf
%   (b) For full rate autocorrelation when calculating sfh():  max(fho)+kcorw+maxlag*kdsamp+hnfllag <= ns
%   (c) For rms ratio window when calculating rr            :  max(fho)+rmsix(end) <= ns
% where max(fho) = (nframe-1)*kframe + spoff

nframe=floor(1+min((nsf-kcorwd-maxlag)/kframed,(ns-spoff-max(kcorw-maxlag*kdsmp-hnfullag,rmsix(end)))/kframe));

% now search for autocorrelation peaks in the downsampled signal

cost=zeros(nframe,ncands);      % cumulative cost
prev=zeros(nframe,ncands);      % traceback pointer
mcands=zeros(nframe,1);         % number of actual candidates excluding voiceless
lagval=repmat(NaN,nframe,ncands-1);    % lag of each voiced candidate
tv=zeros(nframe,3);             % diagnostics: 1=voiceless cost, 2=min voiced cost, 3:cumulative voiceless-min voiced
if doback
    costms=cell(nframe,1);
end

% Main processing loop for each 10 ms frame
nframe_test = 20;
best = zeros (nframe_test, 1);   %create array save the best candidate of the
%nframe_test th frame and the frame converged.
f0 = zeros (nframe_test, 1);
f0(1) = 1;
frmcnt = zeros (nframe_test, 1); %create array which count the number of trace
%back time of each nframe_test;
cmpthF = zeros (nframe_test, 1); %create array save the frame which all
%cands of the nframe_test th frame converged.
%can consider nframe_test as buffer
for iframe = 1: nframe_test       % loop for each frame (~10 ms)

    % Find peaks in the normalized autocorrelation of subsampled (2Khz) speech
    % only keep peaks that are > 30% of highest peak

    fho=(iframe-1)*kframe+spoff;
    sff=sf((iframe-1)*kframed+sfj);
    sffdc=mean(sff(sfi));       % mean of initial correlation window length
    sff=sff-sffdc;              % subtract off the mean
    nccfd=normxcor(sff(1:kcorwd),sff(minlag+1:end));
    [ipkd,vpkd]=v_findpeaks(nccfd,'q');

    % Debugging: execute the line below to plot the autocorrelation peaks.
    v_findpeaks(nccfd,'q'); xlabel(sprintf('Lag = (x+%d)*%g ms',minlag-1,1000*kdsmp/fs)); ylabel('Normalized Cross Correlation'); title (sprintf('Frame %d/%d',iframe,nframe));

    vipkd=[vpkd ipkd];
    vipkd(vpkd<max(vpkd)*candtr,:)=[];          % eliminate peaks that are small
    if size(vipkd,1)
        if size(vipkd,1)>ncands-1
            vipkd=sortrows(vipkd);
            vipkd(1:size(vipkd,1)-ncands+1,:)=[];   % eliminate lowest to leave only ncands-1
        end
        lagcan=round(vipkd(:,2)*kdsmp+lagoff);        % convert the lag candidate values to the full sample rate
        nlcan=length(lagcan);
    else
        nlcan=0;
    end

    % If there are any candidate lag values (nlcan>0) then refine their accuracy at the full sample rate

    if nlcan
        laglist=reshape(repmat(lagcan(:)',nfullag,1)+repmat((-hnfullag:hnfullag)',1,nlcan),nfullag*nlcan,1);
        sfh=s(fho+(1:kcorw+max(lagcan)+hnfullag));
        sfhdc=mean(sfh(sfhi));
        sfh=sfh-sfhdc;
        e0=sum(sfh(sfhi).^2);                     % energy of initial correlation window (only needed to store in tv(:,6)
        lagl2=repmat(lagcan(:)',nfullag+kcorw-1,1)+repmat((1-hnfullag:hnfullag+kcorw)',1,nlcan);
        nccf=normxcor(sfh(1:kcorw),sfh(lagl2),afact);

        [maxcc,maxcci]=max(nccf,[],1);
        vipk=[maxcc(:) lagcan(:)+maxcci(:)-hnfullag-1];
        vipk=vipk(:,[1 2 2]);
        maxccj=maxcci(:)'+nfullag*(0:nlcan-1);    % vector index into nccf array
        msk=mod(maxcci,nfullag-1)~=1 & 2*nccf(maxccj)-nccf(mod(maxccj-2,nfullag*nlcan)+1)-nccf(mod(maxccj,nfullag*nlcan)+1)>0;  % don't do quadratic interpolation for the end ones
        if any(msk)
            maxccj=maxccj(msk);
            vipk(msk,3)=vipk(msk,3)+(nccf(maxccj+1)-nccf(maxccj-1))'./(2*(2*nccf(maxccj)-nccf(maxccj-1)-nccf(maxccj+1)))';
        end
        vipk(maxcc<max(maxcc)*candtr,:)=[];          % eliminate peaks that are small
        if size(vipk,1)>ncands-1
            vipk=sortrows(vipk);
            vipk(1:size(vipk,1)-ncands+1,:)=[];   % eliminate lowest to leave only ncands-1
        end

        % vipk(:,1) has NCCF value, vipk(:,2) has integer peak position, vipk(:,3) has refined peak position

        mc=size(vipk,1);
    else
        mc=0;
    end

    % We now have mc lag candidates at the full sample rate

    mc1=mc+1;               % total number of candidates including "unvoiced" possibility
    mcands(iframe)=mc;      % save number of lag candidates (needed for pitch consistency cost calculation)
    if mc
        lagval(iframe,1:mc)=vipk(:,3)';
        cost(iframe,1)=vobias+max(vipk(:,1));   % voiceless cost
        cost(iframe,2:mc1)=1-vipk(:,1)'.*(1-beta*vipk(:,3)');   % local voiced costs
        tv(iframe,2)=min(cost(iframe,2:mc1));
    else
        cost(iframe,1)=vobias;          % if no lag candidates (mc=0), then the voiceless case is the only possibility
    end
    tv(iframe,1)=cost(iframe,1);
    if iframe>1                         % if it is not the first frame, then calculate pitch consistency and v/uv transition costs
        mcp=mcands(iframe-1);
        costm=zeros(mcp+1,mc1);         % cost matrix: rows and cols correspond to candidates in previous and current frames (incl voiceless)

        % if both frames have at least one lag candidate, then calculate a pitch consistency cost

        if mc*mcp
            lrat=abs(log(repmat(lagval(iframe,1:mc),mcp,1)./repmat(lagval(iframe-1,1:mcp)',1,mc)));
            costm(2:end,2:end)=freqwt*min(lrat,doublec+abs(lrat-log2));  % allow pitch doubling/halving
        end

        % if either frame has a lag candidate, then calculate the cost of voiced/voiceless transition and vice versa

        if mc+mcp
            rr=sqrt((rmswin'*s(fho+rmsix).^2)/(rmswin'*s(fho+rmsix-kdrms).^2)); % amplitude "gradient"
            ss=0.2/(distitar(lpcauto(sp(fho+rmsix),lpcord),lpcauto(sp(fho+rmsix-kdrms),lpcord),'e')-0.8);   % Spectral stationarity: note: Talkin uses Hanning instead of Hamming windows for LPC
            costm(1,2:end)= vtranc+vtrsc*ss+vtrac/rr;   % voiceless -> voiced cost
            costm(2:end,1)= vtranc+vtrsc*ss+vtrac*rr;
            tv(iframe,4:5)=[costm(1,mc1) costm(mcp+1,1)];
        end
        costm=costm+repmat(cost(iframe-1,1:mcp+1)',1,mc1);  % add in cumulative costs
        [costi,previ]=min(costm,[],1);
        cost(iframe,1:mc1)=cost(iframe,1:mc1)+costi;
        prev(iframe,1:mc1)=previ;
    else                            % first ever frame
        costm=zeros(1,mc1); % create a cost matrix in case doing a backward recursion
    end
    if mc
        tv(iframe,3)=cost(iframe,1)-min(cost(iframe,2:mc1));
        tv(iframe,6)=5*log10(e0*e0/afact);
    end
    if doback
        costms{iframe}=costm; % need to add repmatted cost into this
    end
    
% %%do a traceback
       for iiframe = iframe:-1:2 %run respectively from frame nframe_test to the 2nd frame
        count = 0;
        for k = 2: 1: (mcands(iiframe) + 1)
            if prev (iiframe, 1) ~= prev (iiframe, k) || prev (iiframe, 1) == 1
                frmcnt (iframe) = frmcnt (iframe) + 1;
                break;
            else
                count = count + 1;
            end
        end
        % ok
        if count < (mcands(iiframe) - 1)
            continue;
        end
        % ok
        if count == (mcands(iiframe))
            cmpthF (iframe) = iiframe - 1;
            best (iframe) = prev (iiframe, 1);
%  %                     temp = fs/lagval(cmpthF(iframe),mbest(cmpthF(iframe)) - 1);
%  %                     f0(iframe) = temp;
            f0(iframe) = fs/lagval(cmpthF(iframe), best(iframe) - 1);
            fprintf ('Frame %d traceback %d frames!\n Its frequency is %.2f\n\n ', iframe, frmcnt(iframe) + 1, f0(iframe));
            break;
        end
    end   
end
