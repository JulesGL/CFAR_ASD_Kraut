# 📡 Detection Class research project


Welcome to the README for my research project focusing on the article titled "The CFAR Adaptive Subspace Detector is a Scale-Invariant GLRT", written by S.Kraut and L.L. Scharf and published in IEEE Transactions on Signal Processing in 1999.

This project is conducted as part of my coursework in Detection within the Signals and Systems module at ISAE-Supaero, under the guidance of Professor Stéphanie Bidon.

## 📝 Project Overview 

In this research, the goal is to study a precise detector, the Adaptive Subspace Detector (ASD), particularly explored by L.L. Scharf and L.T. McWhorter in their conference paper "Adaptive matched subspace detectors and adaptive coherence estimators" from 1996 before being developed in the article at stake in this report. It was a particular breakthrough in detection cases with unknown covariance matrices expanding Kelly's detector to an unknown noise scaling factor. 

I wrote some examples of basic detectors (Neyman-Pearson, GLRT), implemented the method of the article, developed the computation involved in it and studied its performance.


## 🏗️ Project Structure 

The project is structured as follows:

- **Examples** This folder contains a few example files following cases of scalar and vectorial Neyman-Pearson and Generalized Likelihood Ratio (GLR) test.

- **CFAR_ASD_GLR_main_code** The main code of the simulation of ASD

- **Report** The PDF containing my report on this detection project containing the context studied, some details about the article (development of the test statistic and link to GLRT) and the results found after simulation.

## 🤝 How to Contribute 

Contributions and suggestions are welcome. If you have insights, ideas, or resources that could enhance the project, please feel free to reach out.

## 🙏 References 

Computations made using MatLab 2022b

S. Kraut and L. L. Scharf, "The CFAR adaptive subspace detector is a scale-invariant GLRT," in IEEE Transactions on Signal Processing, vol. 47, no. 9, pp. 2538-2541, Sept. 1999, doi: 10.1109/78.782198.

L. L. Scharf and L. T. McWhorter, "Adaptive matched subspace detectors and adaptive coherence estimators," Conference Record of The Thirtieth Asilomar Conference on Signals, Systems and Computers, Pacific Grove, CA, USA, 1996, pp. 1114-1117 vol.2, doi: 10.1109/ACSSC.1996.599116.

Kelly, E. (1986). An Adaptive Detection Algorithm. IEEE Transactions on Aerospace and Electronic Systems, AES-22(2), 115-127.

K. B. Petersen, & M. S. Pedersen. (2012). The Matrix Cookbook.

Wozencraft, J., & Jacobs, I. (1965). Principles of Communication Engineering by John M. Wozencraft and Irwin Jacobs. Wiley.

Tang, B., Liu, J., Huang, Z., Wang, G., & Fan, F. (2020). Adaptive Target Detection in Gaussian Clutter Edges. IEEE Transactions on Aerospace and Electronic Systems, 56(2), 1662-1673.

D. Xu, P. Addabbo, C. Hao, J. Liu, D. Orlando, & A. Farina (2021). Adaptive strategies for clutter edge detection in radar. Signal Processing, 186, 108127.

