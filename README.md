LEAD-DBS
========

LEAD-DBS is *NOT* intended for clinical use!

### About the adaption we made for our paper
Distance to white matter tracts is associated with deep brain stimulation motor outcome in Parkinson’s disease
Prent et al. 2019 https://doi.org/10.3171/2019.5.JNS1952
1. We added a button to import ExploreDTI data.
2. We added functionality to the visualisation tool within Lead DBS.
   - visualise DTI tracts loaded from ExploreDTI *.mat file
   - calculate distance from contact points to ExploreDTI tracts
   - color DTI tracts according to contact points
   - create histograms to show overview of minimal tract-electrode distance for each contact point.
(code is not yet cleaned up enough to put it on github and contains no manual yet; contact me on w.v.potters@amsterdamumc.nl if you want me to put up the current version and/or if you want to skype/mail about how to use it)

### About Lead-DBS

LEAD-DBS is a MATLAB-toolbox facilitating the: 

- reconstruction of deep-brain-stimulation (DBS) electrodes in the human brain on basis of postoperative MRI and/or CT imaging
- the visualization of localization results in 2D/3D
- a group-analysis of DBS-electrode placement results and their effects on clinical results
- simulation of DBS stimulations (calculation of volume of activated tissue – VAT)
- diffusion tensor imaging (DTI) based connectivity estimates and fiber-tracking from the VAT to other brain regions (connectomic surgery)

LEAD-DBS builds on SPM8/12, especially regarding warping and segmentation procedures.

### Installation

Usually, Lead-DBS can be downloaded from our website (www.lead-dbs.org) in fully functional form.
Alternatively, especially in case you wish to modify and contribute to Lead-DBS, you can

- Clone the Lead-DBS repository from here.
- Download the [necessary datafiles](http://www.lead-dbs.org/release/download.php?id=data) using this link and unzip the downloaded folder into the cloned git repository.
- We’d love to implement your improvements into Lead-DBS – please contact us for direct push access to Github or feel free to add pull-requests to the Lead-DBS repository.
