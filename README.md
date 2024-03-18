# Direct Greens Function

### Program for calculations of the hole Green's function (single particle correlation function) in the t-Jz model on the 2D square lattice within self-avoiding walks approximation.

The main idea of this project is to restrict the Hilbert space of the problem by considering only the self-avoiding paths of the propagating hole. This reduces the self-consistent equations for the Green's function to recursive ones. The benefit of this approach is that much larger system sizes can be investigated (even several hundred sites compared to around 40 achievable through exact diagonalization). At the same time, approximations are minimal for the 2D t-Jz model, where results are not only qualitatively exact but also quantitatively precise. 

### Local correlation function and correlation function with rotational degrees of freedom

This project is split into two parts:

**S**elf-**a**voiding **W**alks (**saw**) -- for calculations of diagonal coefficients of the Green's function. Allows access to e.g. the local spectral function of a single hole in 2D t-Jz model. 

Results for 2D t-Jz model and further readings: [P. Wrzosek, K. Wohlfeld, Phys. Rev. B **103**, 035113 (2021)](https://doi.org/10.1103/PhysRevB.103.035113) ([arXiv:2012.00395](https://arxiv.org/abs/2012.00395)).

**Rot**ational **S**elf-**a**voiding **W**alks (**rotsaw**) -- access to off-diagonal coefficients through the extension of the **saw** case. Allows e.g. for calculations of spectral function with rotational degrees of freedom. 
