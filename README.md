# HybridMD
HybridMD solves the Poisson equation for closely-packed dielectric sphere composites, and is applied in a standard MD for dynamics of dielectric nanoparticles under shift-truncated LJ and electrostatic forces (with polarization effect included).
The hybrid method solves the interface problem for Poisson by combining the analytical methods of image charges and image multipoles with the
spectrally accurate mesh-free method of moments. The resulting linear system is well conditioned and
requires many fewer unknowns on material interfaces as compared with standard methods.
We further apply fast multipole method which reduces the cost from O(N^2) to O(N), where N is the number elements.

   Authors:  
   - Zecheng Gan  (zecheng@nyu.edu) 
   - Shidong Jiang  (shidong.jiang@njit.edu) 
   - Erik Luijten  (luijten@northwestern.edu)
   - Zhenli Xu  (xuzl@sjtu.edu.cn)


## References
Please refer to the following references for more background:

  - Gan, Z., Jiang, S., Luijten, E., & Xu, Z. (2016). A hybrid method for systems of closely spaced dielectric spheres and ions. SIAM Journal on Scientific Computing, 38(3), B375-B395.
   
  - Gan, Z., Wang, Z., Jiang, S., Xu, Z., & Luijten, E. (2019). Efficient dynamic simulations of charged dielectric colloids through a novel hybrid method. The Journal of chemical physics, 151(2), 024112.


## Build Instructions
TBA

## Examples

TBA

## License
Released under the GNU General Public License v3.0.


## Disclaimer
This material is based upon work supported under NSFC Grants 91130012, 11571236 and 21773165; NSF Grants DMR-1121262, DMR-1310211, DMR-1610796 and DMS-1418918; and also the HPC Center of Shanghai Jiao Tong Univeristy. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the authors and do not necessarily reflect the views of the NSFC and NSF.
