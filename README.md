# KPROJ: A Band Unfolding Program
This program is currently interfaced to VASP, quantum espresso, Ab-init, ABACUS and Phonopy. It can do band unfoldings for both bulk and interface systems modelled in supercells. For interfaces, it uses a technique combining FFT and back FFT, which accelerates the calculation significantly. The program allows us to investigate properties of electronic states in any spatial region (defined by zlay1 and zlay2) including the vacuum region for slabs, which is helpful for understanding ARPES and STM/STS experiments. 

Contact: Mingxing Chen (mxchen@hunnu.edu.cn), School of Physics and Electronics, Hunan Normal University, Changsha, Hunan 410081, China.
         Jiaxin Chen (jxchen@hunnu.edu.cn), School of Physics and Electronics, Hunan Normal University, Changsha, Hunan 410081, China.

KPROJ was released under GPL V3. Use of KPROJ should reference:

Layer k-projection and unfolding electronic bands at interfaces, Mingxing Chen and M. Weinert, Phys. Rev. B 98, 245421 (2018).
KPROJ: A Program for Unfolding Electronic and Phononic Bands, Jiaxin Chen and Mingxing Chen, arXiv:2410.10910.

Installation:
cd ~/kproj/src
First, install the FFTW and add its include directory to the Makefile of KPROJ.
then:
1) make fox
2) make 
If the program has been succesfully compiled, an executable named kproj will be found. 

