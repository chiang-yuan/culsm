#include "thermo.h"

// host vectors

double* h_parsum_pe;
double* h_parsum_ke;
double* h_parsum_stress;

// device vectors

double* d_parsum_pe; 		// partial sum of potential energy
double* d_parsum_ke; 		// partial sum of kinetic energy
double* d_parsum_stress;

__global__ void reduce_pe(
	float* ue, float* fa,
	// double* x, float* k, float* r0, int* atom_i, int* atom_j,
	double* parsum_pe, // partial sum of potential energy
	int nbonds)
{
	extern __shared__ double cache[];

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	double temp = 0;

	while (i < nbonds) {

		// int a3_i = atom_i[i] * 3;
        // int a3_j = atom_j[i] * 3;
        
        // double rij = sqrt(pow(x[a3_j] - x[a3_i], 2) + 
		// 	pow(x[a3_j + 1] - x[a3_i + 1], 2) + 
		// 	pow(x[a3_j + 2] - x[a3_i + 2], 2));

		// temp += 1.0/2.0*k[i]*pow(rij - r0[i], 2);
		temp += 1.0/2.0*ue[i]*fa[i];
		
		i += gridDim.x * blockDim.x;
	}
	cache[threadIdx.x] = temp;

	__syncthreads();

	// parallel reduction

	int ihalf = blockDim.x/2;
	while (ihalf != 0) {
		if (threadIdx.x < ihalf) cache[threadIdx.x] += cache[threadIdx.x + ihalf];
		__syncthreads();
		ihalf /= 2;
	}
	if (threadIdx.x == 0) parsum_pe[blockIdx.x] = cache[threadIdx.x];
	__syncthreads();
}

__global__ void reduce_ke(
	float* m, double* v, 
	double* parsum_ke, // partial sum of potential energy
	int natoms)
{
	extern __shared__ double cache[];

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	double temp = 0;
	double v2 = 0;

	while (i < natoms) {

		int a3 = i * 3;
		
		v2 = pow(v[a3], 2) + 
			pow(v[a3 + 1], 2) +
			pow(v[a3 + 2], 2);

		temp += 1.0/2.0*m[i]*v2;
		
		i += gridDim.x * blockDim.x;
	}
	cache[threadIdx.x] = temp;

	__syncthreads();

	// parallel reduction

	int ihalf = blockDim.x/2;
	while (ihalf != 0) {
		if (threadIdx.x < ihalf) cache[threadIdx.x] += cache[threadIdx.x + ihalf];
		__syncthreads();
		ihalf /= 2;
	}
	if (threadIdx.x == 0) parsum_ke[blockIdx.x] = cache[threadIdx.x];
	__syncthreads();
}

__global__ void reduce_stress_comp(
	float* stress, int comp,
	double* parsum_stress, // partial sum of stress
	int natoms)
{
	extern __shared__ double cache[];

	int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	double temp = 0;

	while (i < natoms) {

		int a6 = i * 6;
		
		temp += stress[a6 + comp];
		
		i += gridDim.x * blockDim.x;
	}
	cache[threadIdx.x] = temp;

	__syncthreads();

	// parallel reduction

	int ihalf = blockDim.x/2;
	while (ihalf != 0) {
		if (threadIdx.x < ihalf) cache[threadIdx.x] += cache[threadIdx.x + ihalf];
		__syncthreads();
		ihalf /= 2;
	}
	if (threadIdx.x == 0) parsum_stress[blockIdx.x] = cache[threadIdx.x];
	__syncthreads();
}

Thermo::Thermo(Error *error_)
{
	error = error_;
}

Thermo::~Thermo()
{
}

int Thermo::set(int nthermo_)
{
    nthermo = nthermo_;
	return 0;
}

double Thermo::pe(System & sys)
{
	int threadsPerBlockforBonds = BLOCK_SIZE;
    int blocksPerGridforBonds = (sys.nbonds + threadsPerBlockforBonds - 1)/threadsPerBlockforBonds;

	dim3 dimBlockforBond(threadsPerBlockforBonds, 1, 1);
    dim3 dimGridforBond(blocksPerGridforBonds, 1, 1);

	int bond_double_sm_size = blocksPerGridforBonds*sizeof(double);
	h_parsum_pe = (double*)calloc(blocksPerGridforBonds, sizeof(double));

	cudaMalloc((void**)&d_parsum_pe, bond_double_sm_size);
	cudaMemset(d_parsum_pe, 0, bond_double_sm_size);

	// cudaMemcpy(d_parsum_pe, h_parsum_pe, bond_double_sm_size, cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();

	int sm = threadsPerBlockforBonds*sizeof(double);
	reduce_pe<<<dimGridforBond, dimBlockforBond, sm>>>(
		d_ue, d_fa,
		// d_x, d_ke, d_r0, d_atom_i, d_atom_j,
		d_parsum_pe, // partial sum of potential energy
		sys.nbonds
    );

	cudaDeviceSynchronize();

	cudaMemcpy(h_parsum_pe, d_parsum_pe, bond_double_sm_size, cudaMemcpyDeviceToHost);

	cudaDeviceSynchronize();

	// Add the partial sum of all blocks
	double sum_pe = 0;
	for (int i = 0; i < blocksPerGridforBonds; i++) {
		sum_pe += h_parsum_pe[i];
	}

	// for (b = 0; b < sys.nbonds; b++) {
	// 	btype = sys.bond_type[b];
	// 	// k = sys.bondTypes[itype].coeff[0];
	// 	k = h_k[b];
	// 	r0 = sys.bondTypes[btype - 1].coeff[1];
		
	// 	ai = sys.atom_i[b];
	// 	aj = sys.atom_j[b];
		
	// 	rij = sqrt(pow(sys.x[ai*3] - sys.x[aj*3],2) + 
	// 		pow(sys.x[ai*3 + 1] - sys.x[aj*3 + 1],2) +
	// 		pow(sys.x[ai*3 + 2] - sys.x[aj*3 + 2],2));

	// 	energy += 1.0/2.0*k*pow(rij - r0, 2);
	// }

	cudaFree(d_parsum_pe);

	return sum_pe;
}

double Thermo::ke(System & sys)
{
	int threadsPerBlockforAtoms = BLOCK_SIZE;
    int blocksPerGridforAtoms = (sys.natoms + threadsPerBlockforAtoms - 1)/threadsPerBlockforAtoms;

	dim3 dimBlockforAtom(threadsPerBlockforAtoms, 1, 1);
    dim3 dimGridforAtom(blocksPerGridforAtoms, 1, 1);

	int atom_double_sm_size = blocksPerGridforAtoms*sizeof(double);
	h_parsum_ke = (double*)calloc(blocksPerGridforAtoms, sizeof(double));

	cudaMalloc((void**)&d_parsum_ke, atom_double_sm_size);
	cudaMemset(d_parsum_ke, 0, atom_double_sm_size);
	// cudaMemcpy(d_parsum_ke, h_parsum_ke, atom_double_sm_size, cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();

	int sm = threadsPerBlockforAtoms*sizeof(double);
	reduce_ke<<<dimGridforAtom, dimBlockforAtom, sm>>>(
		d_m, d_v,
		d_parsum_ke, // partial sum of potential energy
		sys.natoms
    );

	cudaDeviceSynchronize();

	cudaMemcpy(h_parsum_ke, d_parsum_ke, atom_double_sm_size, cudaMemcpyDeviceToHost);

	cudaDeviceSynchronize();

	// Add the partial sum of all blocks
	double sum_ke = 0;
	for (int i = 0; i < blocksPerGridforAtoms; i++) {
		sum_ke += h_parsum_ke[i];
	}
	
	// // TODO: parallel reduction
	// double energy = 0;

	// bigint a;
	// int type;
	// double mass;
	// double v2;

	// for (a = 0; a < sys.natoms; a++) {
	// 	type = sys.type[a];
	// 	mass = sys.atomTypes[type - 1].mass;
		
	// 	v2 = pow(sys.v[a*3], 2) + 
	// 		pow(sys.v[a*3 + 1], 2) +
	// 		pow(sys.v[a*3 + 2], 2);

	// 	energy += 1.0/2.0*mass*v2;
	// }

	cudaFree(d_parsum_ke);

	return sum_ke;
}

double Thermo::stress(System & sys, int comp)
{
	int threadsPerBlockforAtoms = BLOCK_SIZE;
    int blocksPerGridforAtoms = (sys.natoms + threadsPerBlockforAtoms - 1)/threadsPerBlockforAtoms;

	dim3 dimBlockforAtom(threadsPerBlockforAtoms, 1, 1);
    dim3 dimGridforAtom(blocksPerGridforAtoms, 1, 1);

	int atom_double_sm_size = blocksPerGridforAtoms*sizeof(double);
	h_parsum_stress = (double*)calloc(blocksPerGridforAtoms, sizeof(double));

	cudaMalloc((void**)&d_parsum_stress, atom_double_sm_size);
	cudaMemset(d_parsum_stress, 0, atom_double_sm_size);

	// cudaMemcpy(d_parsum_stress, h_parsum_stress, atom_double_sm_size, cudaMemcpyHostToDevice);

	cudaDeviceSynchronize();

	int sm = threadsPerBlockforAtoms*sizeof(double);
	reduce_stress_comp<<<dimGridforAtom, dimBlockforAtom, sm>>>(
		d_stress, comp,
		d_parsum_stress, // partial sum of potential energy
		sys.natoms
    );

	cudaDeviceSynchronize();

	cudaMemcpy(h_parsum_stress, d_parsum_stress, atom_double_sm_size, cudaMemcpyDeviceToHost);

	cudaDeviceSynchronize();

	// Add the partial sum of all blocks
	double sum_stress = 0;
	for (int i = 0; i < blocksPerGridforAtoms; i++) {
		sum_stress += h_parsum_stress[i];
	}

	cudaFree(d_parsum_stress);

	return sum_stress;
}

int Thermo::write_thermo(int timestep, System & sys)
{
	printf("%10d \t %16.9e \t %16.9e \t %16.9e\n", 
		timestep, pe(sys), ke(sys), stress(sys, 0));

	return 0;
}


