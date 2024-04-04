# -*- coding: utf-8 -*-
"""
Pre-process CTF MEG data, from raw (.ds), for a choice-reaction time task. 

@author: Sebastian C. Coleman, ppysc6@nottingham.ac.uk
"""

import os.path as op
import numpy as np
import mne

#%% set up paths

root = r"R:\DRS-PSR\Seb\standard_MEG_beamformer\example_data"
data_path = op.join(root, "data")
deriv_path = op.join(root, "derivatives")

subject = "example"
data_fname = op.join(data_path, subject + ".ds")

#%% load data

raw = mne.io.read_raw_ctf(data_fname, preload=True)
print(raw.info)

#%% third order gradiometer

raw.apply_gradient_compensation(grade=3)

#%% broadband filter

raw.filter(l_freq=1, h_freq=45)

#%% get movement parameters and annotate high movement

# get hpi info
chpi_locs = mne.chpi.extract_chpi_locs_ctf(raw, verbose=False)
head_pos = mne.chpi.compute_head_pos(raw.info, chpi_locs, verbose=False)
head_movement_fig = mne.viz.plot_head_positions(head_pos, mode="traces")

# calculate gradient of head movement
head_pos_grad = head_pos.copy()
head_pos_grad[:,1:] = np.gradient(head_pos[:,1:], axis=0)
head_movement_fig = mne.viz.plot_head_positions(head_pos_grad, mode="traces")

# mark bad head movements greater than 1mm / s
thresh = 0.001  # 1mm / s, adjust accordingly
bad_head_bool = head_pos_grad[:,1:4] > thresh
bad_head_bool = np.sum(bad_head_bool, axis=1) > 0

# create annotation
hpi_sfreq = head_pos[1,0] - head_pos[0,0]
bad_head_duration = 5*hpi_sfreq   
bad_head_onset = head_pos[bad_head_bool,0] - 2*hpi_sfreq
bad_head_description = "BAD_head"
bad_head_annot = mne.Annotations(bad_head_onset, bad_head_duration, 
                                 bad_head_description, 
                                 orig_time=raw.info['meas_date'])

#%% mark SQUID resets

thresh = 2e-12   # this usually works, CONSIDER ADJUSTING
squid_annot, bad_chan = mne.preprocessing.annotate_amplitude(
                                    raw, peak=dict(mag=thresh), picks='meg',
                                    bad_percent=5, min_duration=1/raw.info["sfreq"])
squid_annot.onset = squid_annot.onset - 2  # remove 2 s before 
squid_annot.duration = squid_annot.duration + 4  # stretch out annotation
raw.info["bads"].extend(bad_chan)

#%% annotate smaller muscle artefacts etc 

thresh = 5e-13   # this usually works, CONSIDER ADJUSTING
muscle_annot, bad_chan = mne.preprocessing.annotate_amplitude(
                                    raw, peak=dict(mag=thresh), picks='meg',
                                    bad_percent=5, min_duration=1/raw.info["sfreq"])
muscle_annot.onset = muscle_annot.onset - 0.2
muscle_annot.duration = muscle_annot.duration + 0.4

#%% add up annotations

raw.set_annotations(bad_head_annot + squid_annot + muscle_annot)
raw.plot()

#%% plot psd

raw.copy().plot_psd(fmax=45, picks='meg').show()

#%% remove any bad channels by LEFT CLICKING

raw.plot()

#%% ICA (LEFT CLICK ON CARDIAC AND BLINKING TIMECOURSES)

ica = mne.preprocessing.ICA(n_components=20)
ica.fit(raw, reject_by_annotation=True)
ica.plot_components()
ica.plot_sources(raw)

#%% remove bad components (THE ONES YOU CLICKED ON IN PREVIOUS PLOT)

ica.apply(raw)

#%% save out

preproc_fname = subject + "-raw.fif"
raw.save(op.join(deriv_path, preproc_fname), overwrite=True)
