/* FastEddy®: SRC/HYDRO_CORE/CUDA/cuda_coriolisDevice_cu.h 
* ©2016 University Corporation for Atmospheric Research
* 
* This file is licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
* 
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef _NRAD_CUDADEV_CU_H
#define _NRAD_CUDADEV_CU_H

/*nRAD return codes */
#define CUDA_NRAD_SUCCESS    0


/*##############------------------- NRAD submodule variable declarations ---------------------#################*/
extern __constant__ int nRADSelector_d;   /* nRAD selector, (0 = none, 1 = nRAD activated)*/
extern __constant__ float D_turb_d; /*Turbine rotor diameter in nRAD model*/
extern __constant__ float D_hub_d; /*Turbine hub diameter in nRAD model*/
extern __constant__ float z_hh_d; /*Turbine hub height*/
extern __constant__ float x_turb_d; /*Turbine x-location in domain*/
extern __constant__ float y_turb_d; /*Turbine y-location in domain*/
extern __constant__ int nForcesnRAD_d; /*Number of components in force field from the nRAD*/

/*Forces and additional parameters for nRAD parameterization*/
extern float *forces_nRAD_d; /*Forces acting on the flow in nRAD model*/
extern float *sphere_nRAD_d; /*Flag showing possible turbine location in grid (turbine yaws, so this is a maybe)*/
extern float *dist_nRAD_d; /* Distance perpendicular to nRAD */


/*##############-------------- CORIOLIS_CUDADEV submodule function declarations ------------------############*/

/*----->>>>> int cuda_coriolisDeviceSetup();       ---------------------------------------------------------
* Used to cudaMalloc and cudaMemcpy parameters and coordinate arrays, and for the CORIOLIS_CUDA submodule.
*/
extern "C" int cuda_nRADDeviceSetup(); // cuda_coriolisDeviceSetup();

/*----->>>>> extern "C" int cuda_coriolisDeviceCleanup();  -----------------------------------------------------------
* Used to free all malloced memory by the CORIOLIS submodule.
*/
extern "C" int cuda_nRADDeviceCleanup(); //cuda_coriolisDeviceCleanup();

/*----->>>>> __device__ void  cudaDevice_calcCoriolis();  --------------------------------------------------
* This is the cuda version of the calcCoriolis routine from the HYDRO_CORE module
*/
__device__ void cudaDevice_normalDistnRAD(float* dists, float* sphere, float* xPos_d, float* yPos_d, float* zPos_d);
	
__device__ void cudaDevice_activatenRAD(float* forces, float* dists, float* sphere, float* rho,
                                        float* u, float* v, float* Frhs_u, float* Frhs_v, float* Frhs_w);

#endif // _NRAD_CUDADEV_CU_H
