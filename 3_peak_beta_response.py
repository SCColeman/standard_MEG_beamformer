# -*- coding: utf-8 -*-
"""
Extract the location and timecourse of peak beta ERD related to button presses.

@author: Sebastian C. Coleman, ppysc6@nottingham.ac.uk
"""

import os.path as op
import mne
import numpy as np 
from matplotlib import pyplot as plt

### need the following line so 3d plotting works, for some reason
mne.viz.set_3d_options(depth_peeling=False, antialias=False)

#%% set up paths
root = r"R:\DRS-PSR\Seb\standard_MEG_beamformer\example_data"
data_path = op.join(root, "data")
deriv_path = op.join(root, "derivatives")

# data filenames
subject = "example"
data_fname = op.join(deriv_path, subject + "-raw.fif")
fwd_fname = op.join(deriv_path, subject + "-fwd.fif")

# freesurfer filenames
subjects_dir = op.join(root, "subjects_dir")
fs_subject = "fsaverage"

#%% load data

raw = mne.io.read_raw_fif(data_fname)
fwd = mne.read_forward_solution(fwd_fname)
events = mne.find_events(raw, stim_channel="UPPT001")
raw.pick('mag')

#%% epoch based on trigger

event_id = [1, 32]  # trigger of interest
tmin, tmax = -0.5, 1.5
epochs = mne.Epochs(
    raw,
    events,
    event_id,
    tmin,
    tmax,
    baseline=(-0.4, -0.1),
    preload=True,
    reject=dict(mag=4e-12),
    reject_by_annotation=True)

#%% compute covariance of all data, unfiltered

cov = mne.compute_covariance(epochs)
cov.plot(epochs.info)

#%% compute inverse (beamformer weights)

filters = mne.beamformer.make_lcmv(
    raw.info,
    fwd,
    cov,
    reg=0.05, # regularisation, 5% is fairly standard
    noise_cov=None,
    pick_ori="max-power",
    weight_norm="unit-noise-gain",
    rank=None)

#%% get pseudo T

# filter epochs
fband = [8, 13]   # alpha
epochs_filt = epochs.copy().filter(fband[0], fband[1])

# compute active and control covariance of filtered data        
act_min, act_max = -0.1, 0.1
con_min, con_max = 1, 1.2

# make pseudo T of only right button presses (index epochs with trigger val)
active_cov = mne.compute_covariance(epochs_filt["1"], tmin=act_min, 
                                    tmax=act_max, method="shrunk")
control_cov= mne.compute_covariance(epochs_filt["1"], tmin=con_min, 
                                    tmax=con_max, method="shrunk")

stc_active = mne.beamformer.apply_lcmv_cov(active_cov, filters)
stc_control = mne.beamformer.apply_lcmv_cov(control_cov, filters)
pseudoT = (stc_active - stc_control) / (stc_active + stc_control)

# plot
pseudoT.plot(src=fwd['src'], subject=fs_subject,
            subjects_dir=subjects_dir,
            surface="inflated",
            views=["lat", "med"],
            size=600,
            hemi="split",
            smoothing_steps=10,
            time_viewer=False,
            show_traces=False,
            colorbar=True)

#%% peak timecourse in label

# create LCMV generator (index epochs with trigger of interest)
stc_epochs = mne.beamformer.apply_lcmv_epochs(epochs["1"], filters,
                                              return_generator=True)
# get labels from parcellation
parc = "HCPMMP1_combined"
labels = mne.read_labels_from_annot("fsaverage", parc=parc, subjects_dir=subjects_dir)
labels = labels[2:]   # first few labels are not real

# morph labels to subject (obviously this does nothing while you're already
# using fsaverage)
labels = mne.morph_labels(labels, fs_subject, "fsaverage", subjects_dir)

# make atlas label from combination of other labels (i.e., make a bigger ROI)
label_list = [26, 32, 36]
hemi = "lh"
label_name = hemi + " somatosensory"   # CHANGE THIS TO MATCH LABEL_LIST!!!!!!

# combine vertices and pos from labels
vertices = []
pos = []
for l in label_list:
    vertices.append(labels[l].vertices)
    pos.append(labels[l].pos)
vertices = np.concatenate(vertices, axis=0)   
pos = np.concatenate(pos, axis=0)

# sort vertices and pos
vert_order = np.argsort(vertices)
vertices_ordered = vertices[vert_order]
pos_ordered = pos[vert_order,:]

new_label = mne.Label(vertices_ordered, pos_ordered, hemi=hemi, 
                      name=label_name, subject=fs_subject)

# extract peak from pseudoT in your new label
stc_inlabel = pseudoT.in_label(new_label)
label_peak = stc_inlabel.get_peak(mode="abs", vert_as_index=True)[0]

# extract timecourse of peak
n_epochs = len(epochs_filt["1"])
epoch_len = np.shape(epochs[0])[2]
epoch_peak_data = np.zeros((n_epochs,1,epoch_len))
for s,stc_epoch in enumerate(stc_epochs):
    stc_epoch_label = mne.extract_label_time_course(stc_epoch, new_label, 
                                                    fwd['src'], mode=None)
    epoch_peak_data[s,0,:] = stc_epoch_label[0][label_peak,:]
    
# make source epoch object
ch_names = ["peak"]
ch_types = ["misc"]
source_info = mne.create_info(ch_names=ch_names, sfreq=epochs.info["sfreq"],
                              ch_types=ch_types)
source_epochs = mne.EpochsArray(epoch_peak_data, source_info,
                                tmin=epochs.tmin)

# TFR
baseline = (1, 1.2)
freqs = np.arange(1,35)
n_cycles = freqs/2
power = mne.time_frequency.tfr_morlet(source_epochs, freqs=freqs, n_cycles=n_cycles,
                                           use_fft=True, picks="all"
                                           )
power[0].plot(picks="all", baseline=baseline)

# timecourse
source_epochs_filt = source_epochs.copy().filter(fband[0], fband[1], picks="all")
source_epochs_hilb = source_epochs_filt.copy().apply_hilbert(envelope=True, picks="all")
peak_timecourse = source_epochs_hilb.average(picks="all").apply_baseline(baseline)

plt.figure()
plt.plot(peak_timecourse.times, peak_timecourse.get_data()[0], color="black")
plt.ylabel("Oscillatory Power (A.U)")
plt.xlabel("Time (s)")
plt.title(new_label.name)

#%% morph pseudoT to fsaverage (this allows you to average over subjects)

# obviously this will do nothing if you're already using fsaverage, 
# but usually fs_subject would be a specific subject FS reconstruction
fname_fsaverage_src = op.join(subjects_dir, "fsaverage", "bem", 
                              "fsaverage-ico-5-src.fif")
src_to = mne.read_source_spaces(fname_fsaverage_src)

morph = mne.compute_source_morph(
    pseudoT,
    subject_from=fs_subject,
    subject_to="fsaverage",
    src_to=src_to,
    subjects_dir=subjects_dir,
)
pseudoT_morphed = morph.apply(pseudoT)

# plot
pseudoT_morphed.plot(src=src_to, subject="fsaverage",
            subjects_dir=subjects_dir,
            surface="inflated",
            views=["lat", "med"],
            size=600,
            hemi="split",
            smoothing_steps=10,
            time_viewer=False,
            show_traces=False,
            colorbar=True)