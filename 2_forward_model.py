# -*- coding: utf-8 -*-
"""
Calculate forward solution for CTF MEG data, based on FreeSurfer reconstruction.

@author: Sebastian C. Coleman, ppysc6@nottingham.ac.uk
"""

import os.path as op
import mne

### need the following line so 3d plotting works, for some reason
mne.viz.set_3d_options(depth_peeling=False, antialias=False)

#%% set up paths
root = r"R:\DRS-PSR\Seb\standard_MEG_beamformer\example_data"
data_path = op.join(root, "data")
deriv_path = op.join(root, "derivatives")

subject = "example"
data_fname = op.join(deriv_path, subject + "-raw.fif")

#%% load data

raw = mne.io.read_raw_fif(data_fname)
raw.pick('meg')
print(raw.info)

#%% Get FS reconstruction for subject or use fsaverage for quick testing

# directory containing all freesurfer outputs
subjects_dir = op.join(root, "subjects_dir")

# Name of the subject directory in FS subjects_dir
fs_subject = "fsaverage"

plot_bem_kwargs = dict(
    subject=fs_subject,
    subjects_dir=subjects_dir,
    brain_surfaces="white",
    orientation="coronal",
    slices=[50, 100, 150, 200])

mne.viz.plot_bem(**plot_bem_kwargs)

#%% automated coregistration (although ideally use the GUI with fsaverage)

plot_kwargs = dict(
    subject=fs_subject,
    subjects_dir=subjects_dir,
    surfaces="head-dense",
    dig=True,
    meg="sensors",
    show_axes=True,
    coord_frame="meg",
)

coreg = mne.coreg.Coregistration(raw.info, fs_subject, 
                            subjects_dir=subjects_dir)
mne.viz.plot_alignment(raw.info, trans=coreg.trans, **plot_kwargs)
coreg.fit_fiducials()
coreg.set_grow_hair(0)
coreg.fit_icp(10)
coreg.omit_head_shape_points(5 / 1000)
coreg.fit_icp(10)
mne.viz.plot_alignment(raw.info, trans=coreg.trans, **plot_kwargs)

trans_fname = subject + "-trans.fif"
coreg.trans.save(op.join(deriv_path, trans_fname), overwrite=True)

#%% compute source space

# can change oct5 to other surface source space
src = mne.setup_source_space(
    fs_subject, spacing="oct6", add_dist=False, subjects_dir=subjects_dir)
src.plot(subjects_dir=subjects_dir)

src_fname = subject + "-src.fif"
src.save(op.join(deriv_path, src_fname), overwrite=True)

#%% single shell conduction model

conductivity = (0.3,)
model = mne.make_bem_model(
    subject=fs_subject, ico=4,
    conductivity=conductivity,
    subjects_dir=subjects_dir)
bem = mne.make_bem_solution(model)

bem_fname = subject + "-bem.fif"
mne.write_bem_solution(op.join(deriv_path, bem_fname), 
                       bem, overwrite=True)

#%% forward solution

fwd = mne.make_forward_solution(
    raw.info,
    trans=coreg.trans,
    src=src,
    bem=bem,
    meg=True,
    eeg=False
    )
print(fwd)

fwd_fname = subject + "-fwd.fif"
mne.write_forward_solution(op.join(deriv_path, fwd_fname), 
                       fwd, overwrite=True)
