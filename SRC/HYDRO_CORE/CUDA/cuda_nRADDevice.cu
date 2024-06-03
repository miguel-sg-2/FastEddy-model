/* FastEddy®: SRC/HYDRO_CORE/CUDA/cuda_coriolisDevice.cu 
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
/*---nRAD*/ 
__constant__ int nRADSelector_d;   /* nRAD selector, (0 = none, 1 = nRAD activated)*/
__constant__ float D_turb_d; /*Turbine rotor diameter in nRAD model*/
__constant__ float D_hub_d; /*Turbine hub diameter in nRAD model*/
__constant__ float z_hh_d; /*Turbine hub height*/
__constant__ float x_turb_d; /*Turbine x-location in domain*/
__constant__ float y_turb_d; /*Turbine y-location in domain*/
__constant__ int nForcesnRAD_d; /*Number of components in force field from the nRAD*/

/*Forces and additional parameters for nRAD parameterization*/
float *forces_nRAD_d; /*Forces acting on the flow in nRAD model*/
float *sphere_nRAD_d; /*Flag showing possible turbine location in grid (turbine yaws, so this is a maybe)*/
float *dist_nRAD_d; /* Distance perpendicular to nRAD */

/*#################------------ nRAD submodule function definitions ------------------#############*/
/*----->>>>> int cuda_nRADDeviceSetup();       ---------------------------------------------------------
 * Used to cudaMalloc and cudaMemcpy parameters and coordinate arrays, and for the NRAD_CUDA submodule.
*/
extern "C" int cuda_nRADDeviceSetup(){
   int errorCode = CUDA_NRAD_SUCCESS;
   int Nelems;

   cudaMemcpyToSymbol(nRADSelector_d, &nRADSelector, sizeof(int));
   cudaMemcpyToSymbol(D_turb_d, &D_turb, sizeof(float));
   cudaMemcpyToSymbol(D_hub_d, &D_hub, sizeof(float));
   cudaMemcpyToSymbol(z_hh_d, &z_hh, sizeof(float));
   cudaMemcpyToSymbol(x_turb_d, &x_turb, sizeof(float));
   cudaMemcpyToSymbol(y_turb_d, &y_turb, sizeof(float));
   cudaMemcpyToSymbol(nForcesnRAD_d, &nForcesnRAD, sizeof(int));

   Nelems = (Nxp+2*Nh)*(Nyp+2*Nh)*(Nzp+2*Nh);
   fecuda_DeviceMalloc(Nelems*nForcesnRAD*sizeof(float), &forces_nRAD_d);
   fecuda_DeviceMalloc(Nelems*sizeof(float), &sphere_nRAD_d);
   fecuda_DeviceMalloc(Nelems*sizeof(float), &dist_nRAD_d);

   cudaMemcpy(sphere_nRAD_d, sphere_nRAD, Nelems*sizeof(float), cudaMemcpyHostToDevice);

   return(errorCode);
} //end cuda_nRADDeviceSetup()

/*----->>>>> extern "C" int cuda_nRADDeviceCleanup();  -----------------------------------------------------------
Used to free all malloced memory by the NRAD submodule.
*/

extern "C" int cuda_nRADDeviceCleanup(){
   int errorCode = CUDA_NRAD_SUCCESS;

   /* Free any NRAD submodule arrays */
   cudaFree(forces_nRAD_d);
   cudaFree(sphere_nRAD_d);
   cudaFree(dist_nRAD_d);

   return(errorCode);

}//end cuda_nRADDeviceCleanup()


/*----->>>>> __device__ void  cudaDevice_normalDistnRAD();  --------------------------------------------------
* This is the cuda version of the normalDistnRAD routine from the HYDRO_CORE module
*/
__device__ void cudaDevice_normalDistnRAD(float* dists, float* sphere, float* xPos_d, float* yPos_d,float* zPos_d){
	
  int i,j,k,ijk;
  int iStride,jStride,kStride;
  float x_hat[3];
  float dr[3];
  float pi = acosf(-1.0);
  float tilt = 0.0;
  float theta_yaw = 0.0;

  /* Unit vectors that define rotor plane*/
  x_hat[0] = cos(tilt*pi/180)*cos(theta_yaw*pi/180);
  x_hat[1] = cos(tilt*pi/180)*sin(theta_yaw*pi/180);
  x_hat[2] = -1*sin(tilt*pi/180);

  i = (blockIdx.x)*blockDim.x + threadIdx.x;
  j = (blockIdx.y)*blockDim.y + threadIdx.y;
  k = (blockIdx.z)*blockDim.z + threadIdx.z;
  iStride = (Ny_d+2*Nh_d)*(Nz_d+2*Nh_d);
  jStride = (Nz_d+2*Nh_d);
  kStride = 1;
  ijk = i*iStride + j*jStride + k*kStride;

  /*Determine if grid point (at mass-center) is contained in sphere surrounding nRAD*/
  if((i >= iMin_d)&&(i < iMax_d) && (j >= jMin_d)&&(j < jMax_d) && (k >= kMin_d)&&(k < kMax_d)){
    /* Vector defined as the difference between the current grid point and the nacelle */
    dr[0] = xPos_d[ijk]-x_turb_d;
    dr[1] = yPos_d[ijk]-y_turb_d;
    dr[2] = zPos_d[ijk]-z_hh_d;

    /* Calculate distance between the current grid point and the nacelle along x_hat (i.e., the axis of rotation of the turbine)*/
    if(sphere[ijk]>0){
      dists[ijk] = dr[0]*x_hat[0] + dr[1]*x_hat[1] + dr[2]*x_hat[2];
    }
  }
}

/*----->>>>> __device__ void  cudaDevice_activatenRAD();  --------------------------------------------------
* This is the cuda version of the activatenRAD routine from the HYDRO_CORE module
*/
__device__ void cudaDevice_activatenRAD(float* forces, float* dists, float* sphere, float* rho, 
		                        float* u, float* v, float* Frhs_u, float* Frhs_v, float* Frhs_w){

  int i,j,k,ijk,iFld;
  int iStride,jStride,kStride,fldStride;
  float x_hat[3];
  float y_hat[3];
  float z_hat[3];
  float pi = acosf(-1.0);
  float tilt = 0.0;
  float theta_yaw = 0.0;
  float thrust;
  float f_xyz[3];
  float Ct = 0.5;
  float U_local = 0.0;
  float proj_dX = abs(dX_d*cos(theta_yaw*pi/180) + dY_d*sin(theta_yaw*pi/180));
  float a = proj_dX*pow(2*pi,0.5);
  float b = 2*pow(proj_dX,2);
  float distribute;
  float A_turb = pi*pow(0.5*D_turb_d,2);
  int numGridCellsAway = 4;

  /* Unit vectors that define rotor plane*/
  x_hat[0] = cos(tilt*pi/180)*cos(theta_yaw*pi/180);
  x_hat[1] = cos(tilt*pi/180)*sin(theta_yaw*pi/180);
  x_hat[2] = -1*sin(tilt*pi/180);
  y_hat[0] = -1*sin(theta_yaw*pi/180);
  y_hat[1] = cos(theta_yaw*pi/180);
  y_hat[2] = 0.0;
  z_hat[0] = cos(theta_yaw*pi/180)*sin(tilt*pi/180);
  z_hat[1] = sin(theta_yaw*pi/180)*sin(tilt*pi/180);
  z_hat[2] = cos(tilt*pi/180);

  i = (blockIdx.x)*blockDim.x + threadIdx.x;
  j = (blockIdx.y)*blockDim.y + threadIdx.y;
  k = (blockIdx.z)*blockDim.z + threadIdx.z;
  iStride = (Ny_d+2*Nh_d)*(Nz_d+2*Nh_d);
  jStride = (Nz_d+2*Nh_d);
  kStride = 1;
  ijk = i*iStride + j*jStride + k*kStride;
  fldStride = (Nx_d+2*Nh_d)*(Ny_d+2*Nh_d)*(Nz_d+2*Nh_d);

  /*Determine if grid point (at mass-center) is contained in sphere surrounding nRAD*/
  if((i >= iMin_d)&&(i < iMax_d) && (j >= jMin_d)&&(j < jMax_d) && (k >= kMin_d)&&(k < kMax_d)){
    /* Activate nRAD if close to the turbine */  
    if((sphere[ijk]>0) && (abs(dists[ijk]) < numGridCellsAway*proj_dX)){
      /* Thrust force from nRAD model*/
      U_local = pow(u[ijk]*u[ijk] + v[ijk]*v[ijk],0.5);
      thrust = 0.5*Ct*rho[ijk]*U_local*U_local; // [N/m2]
      /* Make sure force corresponds to differential element */
      thrust = thrust*proj_dX*dZ_d/A_turb;  // [N/m2]
      /* Project thrust force onto Cartesian grid*/
      for(iFld=0; iFld < nForcesnRAD_d; iFld++){
        f_xyz[iFld] = -1*thrust*x_hat[iFld];
      }
      /* Distribute force over multiple grid cells*/
      distribute = (1/a)*exp(-1*pow(dists[ijk],2)/b); // [1/m]
      /* Forces excerted by the turbine on the flow are equal in magnitude and opposite in direction */
      for(iFld=0; iFld < nForcesnRAD_d; iFld++){
        forces[fldStride*iFld + ijk] = -1*f_xyz[iFld]*distribute; // [N/m3]
      }
      /* Add to tendency terms */
      Frhs_u[ijk] = Frhs_u[ijk] + forces[fldStride*0 + ijk]; // [N/m3]
      Frhs_v[ijk] = Frhs_v[ijk] + forces[fldStride*1 + ijk]; // [N/m3]
      Frhs_w[ijk] = Frhs_w[ijk] + forces[fldStride*2 + ijk]; // [N/m3]
    } //end if to activate nRAD
  } //end if non-halo cells
}








