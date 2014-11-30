function d = edit_data(d,p,values)
% function d = edit_data(d,p,values)
%
% perform data editing (e.g. sidelobes, previous-ping interference, &c)
%
% version 0.7	last change 16.11.2012

% original A.Thurnherr/M.Visbeck

% changed sidelobe application bins               GK, Sep 2007    0.1-->0.2
% introduced masking of the last bin              GK, Mar 2008    0.2-->0.3
% moving maxbinrange to here                      GK, 24.07.2008  0.3-->0.4
% more than last bin masking                      GK, 13.05.2011  0.4-->0.5
% renamed cosd and sind to cos_d and sin_d        GK, 13.07.2012  0.5-->0.6
% use sfigure instead of figure                   GK, 16.11.2012  0.6-->0.7   

%
% general function info
%
disp(' ')
disp('EDIT_DATA:  apply various editing algorithms to velocities')


%
% store old data
%
d.ts_edited = d.ts;


%
% interpret maxbinrange for edit_mask_x_bins
%
if length(p.maxbinrange)>1 | p.maxbinrange(1)~=0
  max_bin_dn = d.down.Depth_Cells;
  if isfield(d,'up')
    max_bin_up = d.up.Depth_Cells;
    if length(p.maxbinrange)>1
      p.edit_mask_up_bins = [p.maxbinrange(2)+1:max_bin_up];
    else
      p.edit_mask_up_bins = [p.maxbinrange(1)+1:max_bin_up];
    end
  end
  p.edit_mask_dn_bins = [p.maxbinrange(1)+1:max_bin_dn];
  disp('>   Applying maxbinrange to edit_mask')
  disp(['>   No bins higher than ',int2str(p.maxbinrange),' will be used.'])
end


%----------------------------------------------------------------------
% Bin Masking
%----------------------------------------------------------------------
if ~isempty(p.edit_mask_dn_bins) | ~isempty(p.edit_mask_up_bins)

  nbad = 0;
  if length(d.zu) > 0
    for bi=1:length(p.edit_mask_up_bins)
      bn = length(d.zu)+1 - p.edit_mask_up_bins(bi);
      nbad = nbad + length(find(isfinite(d.weight(bn,:))));
      d.weight(bn,:) = NaN; 
      d.ts_edited(bn,:) = NaN;
    end
  end
  if length(d.zd) > 0
    for bi=1:length(p.edit_mask_dn_bins)
      bn = length(d.zu) + p.edit_mask_dn_bins(bi);
      nbad = nbad + length(find(isfinite(d.weight(bn,:))));
      d.weight(bn,:) = NaN; 
      d.ts_edited(bn,:) = NaN;
    end
  end

  disp(sprintf('    Bin masking               : set %d weights to NaN',nbad));

end % if bin masking enabled

%----------------------------------------------------------------------
% Side-Lobe Contamination
%----------------------------------------------------------------------

if p.edit_sidelobes

  nbad = 0;
  
  % first, the uplooker: d.z is -ve distance of ADCP from surface;
  % Cell_length is in cm, i.e. 0.015*Cell_length is 1.5 x bin size
  % in m --- the same value used by Firing's software
  
  if length(d.zu > 0)
  
%    for b=1:length(d.zu)+length(d.zd)
% this loops over all bins (up and down), changing that
% to only up, GK Sep 2007
    zlim = d.izm;
    for b=d.izu
      zlim(b,:) = (1 - cos_d(d.up.Beam_angle)) * d.z ...
		- 0.015*d.up.Cell_length;
    end
    ibad = find(d.izm > zlim);
    nbad = nbad + length(find(isfinite(d.weight(ibad))));
    d.weight(ibad) = NaN; 
    d.ts_edited(ibad) = NaN;
  
  end

  % now, the downlooker: p.zbottom is the +ve depth of the sea bed; therefore,
  % -d.z - p.zbottom is the -ve distance from the sea bed
 
% same as above, changing to only down, GK Sep 2007
  zlim = d.izm;
  for b = d.izd
    zlim(b,:) = -p.zbottom ...
	      + (1 - cos_d(d.down.Beam_angle)) * (d.z+p.zbottom) ...
	      + 0.015*d.down.Cell_length;
  end
  ibad = find(d.izm < zlim);
  nbad = nbad + length(find(isfinite(d.weight(ibad))));
  d.weight(ibad) = NaN; 
  d.ts_edited(ibad) = NaN;
  
  disp(sprintf('    Side-lobe contamination   : set %d weights to NaN',nbad));

end %if p.edit_sidelobes

%----------------------------------------------------------------------
% Time-Domain Spike Filter
%----------------------------------------------------------------------

if ~isnan(p.edit_spike_filter_max_curv)

  nbad = 0;
  for b=[1:size(d.ts,1)]
    dummy = d.ts(b,:)/rms(d.ts(b,:));
    ibad = find(diff(diff(dummy)) < -1*p.edit_spike_filter_max_curv) + 1;
    nbad = nbad + length(find(isfinite(d.weight(b,ibad))));
    d.weight(b,ibad) = NaN; 
    d.ts_edited(b,ibad) = NaN;
  end
  disp(sprintf('    Spike filter              : set %d weights to NaN',nbad));

end %if p.edit_spike_filter


%----------------------------------------------------------------------
% Previous-Ping Interference
%----------------------------------------------------------------------

if p.edit_PPI

  % NB: at present, PPI filtering is only implemented for the downlooker
  
  nbad = 0;
  
  % calc ping-intervals; dt(1) contains the difference (in seconds)
  % between the first two pings (t(2) - t(1)).
  
  dt = diff(d.time_jul)*86400;
  
  % use the mean sound speed below the approximate expected PPI depth;
  % this is anal but not very expensive; using 1500m/s would be nearly
  % as good.
  
  if isfield(d,'ctdprof_z') & isfield(d,'ctdprof_ss')
    guess_z = -p.zbottom + 1500 * nmean(dt) / 2;
    SS = nmean(d.ctdprof_ss(find(d.ctdprof_z > -guess_z)));
    if ~isfinite(SS), 
      SS = 1500; 
    end
  else
    SS = 1500;
  end
  
  % calculate the depth limits to remove the PPI for all (but the first)
  % ensembles; the beam-angle limits were found to be too conservative
  % and were replaced by a nominal layer_tickness.
  
  %PPI_min_beam_angle = 0.0 * d.down.Beam_angle;
  %PPI_max_beam_angle = 1.2 * d.down.Beam_angle;
  %PPI_max_z = -p.zbottom + SS * dt/2 * cos(pi/180 * PPI_min_beam_angle);
  %PPI_min_z = -p.zbottom + SS * dt/2 * cos(pi/180 * PPI_max_beam_angle);
  
  PPI_hab = SS * dt/2 * cos_d(d.down.Beam_angle);
  PPI_hab(find(PPI_hab > p.edit_PPI_max_hab)) = inf;
  
  PPI_min_z = -p.zbottom + PPI_hab - p.edit_PPI_layer_thickness / 2;
  PPI_max_z = PPI_min_z + p.edit_PPI_layer_thickness;
  
  % remove the contaminated data from the downlooker bins
  
  for b=d.izd
    ibad = find(d.izm(b,2:end) > PPI_min_z & d.izm(b,2:end) < PPI_max_z) + 1;
    nbad = nbad + length(find(isfinite(d.weight(b,ibad))));
    d.weight(b,ibad) = NaN; 
    d.ts_edited(b,ibad) = NaN;
  end
  disp(sprintf('    Previous-ping interference: set %d weights to NaN',nbad));

end %if p.edit_PPI

%----------------------------------------------------------------------
% Ensemble Skipping
%----------------------------------------------------------------------

if ~isempty(p.edit_skip_ensembles)

  nskipped = 0;
  iskip = [];
  for i=1:length(p.edit_skip_ensembles)
    if p.edit_skip_ensembles(i)
      iskip = [iskip i:length(p.edit_skip_ensembles):length(d.time_jul)];
    end
  end

  for b=1:length(d.zd)+length(d.zu)
    nskipped = nskipped + length(find(isfinite(d.weight(b,iskip))));
    d.weight(b,iskip) = NaN; 
    d.ts_edited(b,iskip) = NaN;
  end

  disp(sprintf('    Ensemble skipping         : set %d weights to NaN',nskipped));

end % if p.edit_skip_ensembles enabled


%----------------------------------------------------------------------
% Last Bin Masking
%----------------------------------------------------------------------

if p.edit_mask_last_bin(1)~=0

  if length(p.edit_mask_last_bin)>1
    mask_dn = p.edit_mask_last_bin(1);
    mask_up = p.edit_mask_last_bin(2);
  else
    mask_dn = p.edit_mask_last_bin;
    mask_up = p.edit_mask_last_bin;
  end
 
  
  nbad = 0;
  dummy = d.weight;
  for n=1:size(d.weight,2)
%    ind = find(~isnan(d.weight(:,n)));
    ind = find(d.weight(:,n)>0);
    if ~isempty(ind) 
%      if ind(1)~=1
%        d.weight(ind(1),n) = nan;
%      end
%      if ind(end)~=size(d.weight,1)
%        d.weight(ind(end),n) = nan;
%      end
      for m=1:mask_up
        if m<=length(ind)
          d.weight(ind(m),n) = nan;
        end
      end
      for m=1:mask_dn
        if m<=length(ind)
          d.weight(ind(end+1-m),n) = nan;
        end
      end
    end
  end
else
  if mask_up>0
    mask = d.weight(d.izu,:)*0+1;
    for n=1:size(mask,2)
      ind = find(~isnan(mask(:,n)));
      if ~isempty(ind)
        ind = ind - mask_up;
        if ind>0
          mask(ind+1:end,n) = nan;
        end
      end
    end
    d.weight(d.izu,:) = d.weight(d.izu,:).*mask;
  end
  if mask_dn>0
    mask = d.weight(d.izd,:)*0+1;
    for n=1:size(mask,2)
      ind = find(~isnan(mask(:,n)));
      if ~isempty(ind)
        ind = ind - mask_dn;
        if ind>0
          mask(ind+1:end,n) = nan;
        end
      end
    end
    d.weight(d.izd,:) = d.weight(d.izd,:).*mask;
  end
  nbad = length(find(isnan(dummy)-isnan(d.weight)));

  disp(sprintf('    Last Bin masking           : set %d weights to NaN',nbad));

end % if bin masking enabled


%----------------------------------------------------------------------
% Smoothed minimum correlation threshold
%----------------------------------------------------------------------

if any(p.minimum_correlation_threshold>0)
  if values.up==0
    cm = d.cm;
    bad = (cm<p.minimum_correlation_threshold(1));
    firstbad = size(bad,1)*ones(1,size(bad,2));
    for n=1:size(bad,2)
      dummy = find(bad(:,n));
      if ~isempty(dummy)
        firstbad(n) = dummy(1);
      end
    end
    firstbad = round(meanfilt(firstbad,50));
    newweight = d.weight;
    for n=1:length(firstbad)
      if firstbad(n)<=size(bad,1)
        newweight(firstbad(n):end,n) = nan;
      end
    end
    nbad = sum(isnan(newweight(:)))-sum(isnan(d.weight(:)));
    d.weight = newweight;
  else
    if length(p.minimum_correlation_threshold)==1
      p.minimum_correlation_threshold = p.minimum_correlation_threshold*[1,1];
    end
    % first do down-part
    cm = d.cm(d.izd,:);
    bad = (cm<p.minimum_correlation_threshold(1));
    firstbad = size(bad,1)*ones(1,size(bad,2));
    for n=1:size(bad,2)
      dummy = find(bad(:,n));
      if ~isempty(dummy)
        firstbad(n) = dummy(1);
      end
    end
    firstbad = round(meanfilt(firstbad,50));
    newweight = d.weight(d.izd,:);
    for n=1:length(firstbad)
      if firstbad(n)<=size(bad,1)
        newweight(firstbad(n):end,n) = nan;
      end
    end
    nbad1 = sum(isnan(newweight(:)))-sum(sum(isnan(d.weight(d.izd,:))));
    d.weight(d.izd,:) = newweight;
    % then do up part
    cm = d.cm(d.izu,:);
    bad = (cm<p.minimum_correlation_threshold(2));
    firstbad = size(bad,1)*ones(1,size(bad,2));
    for n=1:size(bad,2)
      dummy = find(bad(:,n));
      if ~isempty(dummy)
        firstbad(n) = dummy(1);
      end
    end
    firstbad = round(meanfilt(firstbad,50));
    newweight = d.weight(d.izu,:);
    for n=1:length(firstbad)
      if firstbad(n)<=size(bad,1)
        newweight(firstbad(n):end,n) = nan;
      end
    end
    nbad = nbad1 + sum(isnan(newweight(:)))-sum(sum(isnan(d.weight(d.izu,:))));
    d.weight(d.izu,:) = newweight;
  end
  disp(sprintf('    Minimum corr. threshold   : set %d weights to NaN',nbad));
end

%----------------------------------------------------------------------
% Plot Results of Editing
%----------------------------------------------------------------------

bin_no = [0];
if length(d.zu) > 0 
  bin_no = [-length(d.zu):1 bin_no]; 
end
if length(d.zd) > 0 
  bin_no = [bin_no 1:length(d.zd)]; 
end

sfigure(2);
clf;
orient landscape;
colormap([[1 1 1]; jet(128)]);

subplot(2,1,1);
imagesc([1:size(d.ts,2)],bin_no,...
	[d.ts(1:length(d.zu),:); ...
	 ones(1,size(d.ts,2))*NaN; ...
	 d.ts(size(d.ts,1)-length(d.zd)+1:end,:)...
        ]);
csc = caxis;
colorbar
xlabel('Ensemble #');
ylabel('Bin #');
title('Before Data Editing');

subplot(2,1,2);
imagesc([1:size(d.ts,2)],bin_no,...
	[d.ts_edited(1:length(d.zu),:); ...
	 ones(1,size(d.ts,2))*NaN; ...
	 d.ts_edited(size(d.ts,1)-length(d.zd)+1:end,:)...
        ]);
csc = caxis;
colorbar
xlabel('Ensemble #');
ylabel('Bin #');
title('After Data Editing');

streamer([p.name,'  Figure 14']);
hgsave('tmp/14')


ind = find(isnan(d.ts_edited));
d.cm_edited = d.cm;
d.cm_edited(ind) = nan;
sfigure(2);
clf;
orient landscape;
colormap([[1 1 1]; jet(128)]);

subplot(2,1,1);
imagesc([1:size(d.cm,2)],bin_no,...
	[d.cm(1:length(d.zu),:); ...
	 ones(1,size(d.cm,2))*NaN; ...
	 d.cm(size(d.cm,1)-length(d.zd)+1:end,:)...
        ]);
csc = caxis;
colorbar
xlabel('Ensemble #');
ylabel('Bin #');
title('Before Data Editing');

subplot(2,1,2);
imagesc([1:size(d.cm,2)],bin_no,...
	[d.cm_edited(1:length(d.zu),:); ...
	 ones(1,size(d.cm,2))*NaN; ...
	 d.cm_edited(size(d.cm,1)-length(d.zd)+1:end,:)...
        ]);
csc = caxis;
colorbar
xlabel('Ensemble #');
ylabel('Bin #');
title('After Data Editing');

streamer([p.name,'  Figure 15']);
hgsave('tmp/15')
