We recognize that integrating MRSI into a clinical workflow is challenging and often requires specialized expertise that may not be widely available.

This pipeline is intended to **bridge the gap between research and clinical application**. We provide a validated workflow designed to help the broader clinical community adopt multi-voxel MRSI more consistently and efficiently.

## Recommended MRSI preprocessing steps

1) **Frequency and phase correction (alignment)**  
Because the scanner’s B0 magnetic field is subject to temporal drift—leading to frequency and phase shifts in typical MRS experiments—we recommend applying frequency and phase correction to minimize peak broadening, spectral lineshape distortion, and loss of SNR. This step should be performed for both the **water reference** and **metabolite** acquisitions (unless no misalignment is observed).  
This can be performed using `fsl_mrs_proc mrsi-align` (FSL-MRS toolbox).

Before applying this step, we recommend visually inspecting the data to confirm the presence of frequency and/or phase drift (e.g., shifting peak positions or inconsistent phase across dynamics).

2) **Eddy current correction (ECC)**  
Eddy currents arise from rapid gradient switching and can produce short-lived fluctuations in the B0 field. These fluctuations may introduce time-dependent spectral distortions and reduce spectral quality if left uncorrected.  
This step can be performed using `fsl_mrs_proc ecc` (FSL-MRS toolbox) and typically uses the **unsuppressed water reference** as the correction reference.

3) **Residual water removal**
We recommend applying residual water removal to the metabolite data to reduce any remaining water signal that could obscure nearby metabolite peaks or interfere with spectral fitting and quantification.
