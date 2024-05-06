Jie He, Haoran Zhang, Yimeng Li, Guanghui Li, Siao Lei, Zhumei Qian, Fei Xiong, Yuan Feng, Tao Zhu, Yu An, and Jie Tian. "Sequential Scan-Based Single-Dimension Multi-Voxel System Matrix Calibration for Open-Sided Magnetic Particle Imaging," in IEEE Transactions on Medical Imaging, 2024.

**************************************************************************************************************************************

Simulation
DataGeneration_... : Codes for generate simulated data
Reco_... : Case codes for reconstruct images with SS_SDMVC methods

**************************************************************************************************************************************

Experiment: Please download from https://drive.google.com/file/d/1mpblqE9ghoqaY4F_hny1AXXuMLgRWIIE/view?usp=drive_link
2D data
Sz_reco: System matrix for reconstruction, 200 frequencies were used for populate the system matrix. Before reconstruction, frequency selection was required.
Uz_reco1-Uz_reco4: Signal vectors for reconstruction.
dz_reco1-dz_reco4: Noise vectors for frequency selection with SNR>3.

3D data
Sz_reco_1-Sz_reco_16: 2D system matries for |z[i3]-z[j3]| ranging from 0 to 30 mm with 2 mm steps. 200 frequencies were used for population the system matrix. Before reconstruction frequency selection was required.
Uz_reco: Signal vector for reconstruction.
dz_reco: Noise vector for frequency selection with SNR>3.

**************************************************************************************************************************************

For questions/comments please send an email to: jieh@buaa.edu.cn
