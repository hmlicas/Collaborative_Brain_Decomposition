# Collaborative Brain Decomposition and Brain Decoding

# Reference
Hongming Li, Theodore D. Satterthwaite, and Yong Fan. "Large-scale sparse functional networks from resting state fMRI." NeuroImage 156 (2017): 1-13.
Hongming Li and Yong Fan. "Interpretable, highly accurate brain decoding of subtly distinct brain states fromfunctional MRI using intrinsic functional networks and long short-term memory recurrent neural networks”, NeuroImage 202, 116059
Hongming Li and Yong Fan. “Brain Decoding from Functional MRI Using Long Short-Term Memory Recurrent NeuralNetworks”. In: Frangi A., Schnabel J., Davatzikos C., Alberola-López C., Fichtinger G. (eds) Medical Image Computing and Computer Assisted Intervention – MICCAI 2018. MICCAI 2018. Lecture Notes in Computer Science, vol 11072. Springer, Cham
Hongming Li, Xiaofeng Zhu, Yong Fan. “Identification of Multi-scale Hierarchical Brain Functional NetworksUsing Deep Matrix Factorization”. In: Frangi A., Schnabel J., Davatzikos C., Alberola-López C., Fichtinger G. (eds) Medical Image Computing and Computer Assisted Intervention – MICCAI 2018. MICCAI 2018. Lecture Notes in Computer Science, vol 11072. Springer, Cham

# Introduction
We develop a data-driven method for detecting subject-specific functional networks (FNs) while establishing group level correspondence. Our method simultaneously computes subject-specific FNs for a group of subjects regularized by group sparsity, to generate subject-specific FNs that are spatially sparse and share common spatial patterns across subjects. Our method is built upon nonnegative matrix decomposition techniques, enhanced by a data locality regularization term that makes the decomposition robust to imaging noise and improves spatial smoothness and functional coherences of the subject specific FNs.
  This method has been extended for identifying subject-specific FNs at multiple spatial scales with a hierarchical organization from resting state fMRI data. Our method is built upon a deep semi-nonnegative matrix factorization framework to jointly detect the FNs at multiple scales with a hierarchical organization, enhanced by group sparsity regularization that helps identify subject-specific FNs without loss of inter-subject comparability.
  The subject specific FNs have been adopted in a deep learning framework for brain decoding. Particularly, subject-specific FNs are computed from resting-statefMRI data and are used to characterize functional signals of task fMRI data with a compact representation for building brain decoding models, and long short-term memory (LSTM) recurrent neural networks (RNNs) are adopted to learn brain decoding mappings between functional profiles and brain states.

# Code
For usage, please refer to how_to_use.txt.
