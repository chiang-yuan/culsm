#include "run.h"
#include "thermo.h"

// host vectors

float* h_m;         // atom masses
float* h_k;         // bond stiffness
float* h_r0;
float* h_rc; 
float* h_stress;   // atom stresses

// double* h_parsum_pe;

// device vectors

int* d_type;        // atom types
float* d_m;         // atom masses
double* d_x;        // atom coordinates
double* d_v;        // atom velocities
double* d_f;        
double* d_a;
float* d_stress;   // atom stresses
// double* d_a_old;    // atom old accerlerations
// double* d_a_new;    // atom new accerlerations

float* d_k;         // bond stiffness
float* d_r0;
float* d_rc;
int* d_atom_i;      // bonded atom i
int* d_atom_j;      // bonded atom j

// double* d_parsum_pe; 		// partial sum of potential energy


__global__ void fix_rigidmove(int* type, double* x, double* v, double* a,
                              int fix_type, double vx, double vy, double vz,
                              float dt, int n) 
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i < n && type[i/3] == fix_type) {
        
        if (i%3 == 0) x[i] = x[i] + vx * dt;
        else if (i%3 == 1) x[i] = x[i] + vy * dt;
        else if (i%3 == 2) x[i] = x[i] + vz * dt;
        
        v[i] = 0.0;
        a[i] = 0.0;
        // a_old[i] = 0.0;
        // a_new[i] = 0.0;
    }

    __syncthreads();
}

__global__ void verlet_update_pos(double* x, double* v,
                                  float dt, int n) 
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n) {
        x[i] = x[i] + v[i] * dt;
    }

    __syncthreads();
}

__global__ void verlet_update_vel(double* v, double* a,
                                  float dt, int n) 
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n) {
        v[i] = v[i] + a[i] * dt;
    }

    __syncthreads();
}

__global__ void calculate_force(double* x, double* f,
                                float* k, float* r0, float* rc, int* atom_i, int* atom_j,
                                int natoms, int nbonds)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i < natoms) {
        f[i*3] = 0;
        f[i*3 + 1] = 0;
        f[i*3 + 2] = 0;
    }
    __syncthreads();

    if (i < nbonds) {

        int ai = atom_i[i];
        int aj = atom_j[i];
        int index_i = ai * 3;
        int index_j = aj * 3;
        
        double r_ij = sqrt(pow(x[index_j] - x[index_i], 2) + 
                          pow(x[index_j + 1] - x[index_i + 1], 2) + 
                          pow(x[index_j + 2] - x[index_i + 2], 2));

        if (r_ij > rc[i]) k[i] = 0;

        if (r_ij < EPS) r_ij += EPS;

        double fix = -k[i] * (r_ij - r0[i]) * (x[index_i] - x[index_j])/r_ij;
        double fiy = -k[i] * (r_ij - r0[i]) * (x[index_i + 1] - x[index_j + 1])/r_ij;
        double fiz = -k[i] * (r_ij - r0[i]) * (x[index_i + 2] - x[index_j + 2])/r_ij;

        atomicAdd(&f[index_i], fix);
        atomicAdd(&f[index_i + 1], fiy);
        atomicAdd(&f[index_i + 2], fiz);
        
        atomicAdd(&f[index_j], -fix);
        atomicAdd(&f[index_j + 1], -fiy);
        atomicAdd(&f[index_j + 2], -fiz);
    }

    __syncthreads();
}

__global__ void verlet_update_acc(double* a, float* m, double* f, 
                                  int natoms)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < natoms) {
        int i3 = i * 3;
        a[i3] = f[i3]/m[i];
        a[i3 + 1] = f[i3 + 1]/m[i];
        a[i3 + 2] = f[i3 + 2]/m[i];
    }
    
    __syncthreads();
}

__global__ void stress_virial_kin(
    float* stress, float* m, double* v,
    int natoms)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i < natoms) {
        int i6 = i * 6;
        int i3 = i * 3;

        stress[i6] = 0;        // xx
        stress[i6 + 1] = 0;    // yy
        stress[i6 + 2] = 0;    // zz
        stress[i6 + 3] = 0;    // xy
        stress[i6 + 4] = 0;    // yz
        stress[i6 + 5] = 0;    // zx

        stress[i6] -= m[i] * v[i3] * v[i3];
        stress[i6 + 1] -= m[i] * v[i3 + 1] * v[i3 + 1];
        stress[i6 + 2] -= m[i] * v[i3 + 2] * v[i3 + 2];
        stress[i6 + 3] -= m[i] * v[i3] * v[i3 + 1];
        stress[i6 + 4] -= m[i] * v[i3 + 1] * v[i3 + 2];
        stress[i6 + 5] -= m[i] * v[i3 + 2] * v[i3];

    }

    __syncthreads();
}


__global__ void stress_virial_pot(
    float* stress, 
    float* m, double* x, double* v,
    float* k, float* r0, float* rc, 
    int* atom_i, int* atom_j,
    int nbonds)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < nbonds) {

        int i3_i = atom_i[i] * 3;
        int i3_j = atom_j[i] * 3;

        int i6_i = atom_i[i] * 6;
        int i6_j = atom_j[i] * 6;
        
        double r_ij = sqrt(pow(x[i3_j] - x[i3_i], 2) + 
                          pow(x[i3_j + 1] - x[i3_i + 1], 2) + 
                          pow(x[i3_j + 2] - x[i3_i + 2], 2));

        // if (r_ij > rc[i]) k[i] = 0;

        if (r_ij < EPS) r_ij += EPS;

        double fix = -k[i] * (r_ij - r0[i]) * (x[i3_i] - x[i3_j])/r_ij;
        double fiy = -k[i] * (r_ij - r0[i]) * (x[i3_i + 1] - x[i3_j + 1])/r_ij;
        double fiz = -k[i] * (r_ij - r0[i]) * (x[i3_i + 2] - x[i3_j + 2])/r_ij;

        atomicAdd(&stress[i6_i], 1.0 / 2.0 * (x[i3_j] - x[i3_i]) * fix);   // sxx
        atomicAdd(&stress[i6_i + 1], 1.0 / 2.0 * (x[i3_j + 1] - x[i3_i + 1]) * fiy);   // syy
        atomicAdd(&stress[i6_i + 2], 1.0 / 2.0 * (x[i3_j + 2] - x[i3_i + 2]) * fiz);   // szz
        atomicAdd(&stress[i6_i + 3], 1.0 / 2.0 * (x[i3_j] - x[i3_i]) * fiy);   // sxy
        atomicAdd(&stress[i6_i + 4], 1.0 / 2.0 * (x[i3_j + 1] - x[i3_i + 1]) * fiz);   // syz
        atomicAdd(&stress[i6_i + 5], 1.0 / 2.0 * (x[i3_j + 2] - x[i3_i + 2]) * fix);   // szx

        atomicAdd(&stress[i6_j], 1.0 / 2.0 * (x[i3_j] - x[i3_i]) * fix);   // sxx
        atomicAdd(&stress[i6_j + 1], 1.0 / 2.0 * (x[i3_j + 1] - x[i3_i + 1]) * fiy);   // syy
        atomicAdd(&stress[i6_j + 2], 1.0 / 2.0 * (x[i3_j + 2] - x[i3_i + 2]) * fiz);   // szz
        atomicAdd(&stress[i6_j + 3], 1.0 / 2.0 * (x[i3_j] - x[i3_i]) * fiy);   // sxy
        atomicAdd(&stress[i6_j + 4], 1.0 / 2.0 * (x[i3_j + 1] - x[i3_i + 1]) * fiz);   // syz
        atomicAdd(&stress[i6_j + 5], 1.0 / 2.0 * (x[i3_j + 2] - x[i3_i + 2]) * fix);   // szx
    }

    __syncthreads();
}


Run::Run(Error *error_)
{
    error = error_;
}

Run::~Run()
{
    
}


int Run::verlet(
    float dt, int timesteps, System &sys, Dump &dump, Thermo &thermo, std::vector<Fix> &fixes) 
{

    printf("Velocity Verlet run for %d steps with timestep size %f...\n", timesteps, dt);

    // TODO: device id
    cudaError_t err = cudaSuccess;
    int gid = 0;
    err = cudaSetDevice(gid);
    if (err != cudaSuccess) {
        char* buffer = new char[MAX_STRING];
        sprintf(buffer, "Cannot select GPU with device ID = %d\n", gid);
        error->message(buffer, -1);
    } 
    else {
        printf("Set GPU with device ID = %d\n", gid);
    }

    // TODO: format camel name of variable
    int threadsPerBlockforAtoms = BLOCK_SIZE;
    int blocksPerGridforAtoms = (sys.natoms*3 + threadsPerBlockforAtoms - 1)/threadsPerBlockforAtoms;
    int blocksPerGridforAtoms6 = (sys.natoms*6 + threadsPerBlockforAtoms - 1)/threadsPerBlockforAtoms;

    int threadsPerBlockforBonds = BLOCK_SIZE;
    int blocksPerGridforBonds = (sys.nbonds + threadsPerBlockforBonds - 1)/threadsPerBlockforBonds;

    dim3 dimBlockforAtom(threadsPerBlockforAtoms, 1, 1);
    dim3 dimGridforAtom(blocksPerGridforAtoms, 1, 1);
    dim3 dimGridforAtom6(blocksPerGridforAtoms6, 1, 1);
    printf("\tthread size for atoms: (%d, %d, %d)\n", dimBlockforAtom.x, dimBlockforAtom.y, dimBlockforAtom.z);
    printf("\tgrid size for atoms: (%d, %d, %d)\n", dimGridforAtom.x, dimGridforAtom.y, dimGridforAtom.z);

    dim3 dimBlockforBond(threadsPerBlockforBonds, 1, 1);
    dim3 dimGridforBond(blocksPerGridforBonds, 1, 1);
    printf("\tthread size for bonds: (%d, %d, %d)\n", dimBlockforBond.x, dimBlockforBond.y, dimBlockforBond.z);
    printf("\tgrid size for bonds: (%d, %d, %d)\n", dimGridforBond.x, dimGridforBond.y, dimGridforBond.z);

    // memory size

    int atom_double3_size = sys.natoms*3*sizeof(double);
    int atom_float6_size = sys.natoms*6*sizeof(float);
    int atom_float_size = sys.natoms*sizeof(float);
    int atom_int_size = sys.natoms*sizeof(int);
    int bond_float_size = sys.nbonds*sizeof(float);
    int bond_int_size = sys.nbonds*sizeof(int);

    int sm_size = threadsPerBlockforAtoms*3*sizeof(float);

    // initialize atom mass m

    h_m = (float *)malloc(atom_float_size);

    for (int a = 0; a < sys.natoms; a++) {
        for (int at = 0; at < sys.no_atom_types; at++) {
            if (sys.type[a] == at + 1) {
                h_m[a] = sys.atomTypes[at].mass;
                break;
            }
        }
    }

    // allocate and zero-initialize atom stresses

    h_stress = (float *) calloc(sys.natoms*6, sizeof(float));

    // initialize bond stiffness k

    h_k = (float *)malloc(bond_float_size);
    h_r0 = (float *)malloc(bond_float_size);
    h_rc = (float *)malloc(bond_float_size);

    for (int b = 0; b < sys.nbonds; b++) {
        for (int bt = 0; bt < sys.no_bond_types; bt++) {
            if (sys.bond_type[b] == bt + 1) {
                h_k[b] = sys.bondTypes[bt].coeff[0];
                h_r0[b] = sys.bondTypes[bt].coeff[1];
                h_rc[b] = sys.bondTypes[bt].coeff[2];
                break;
            }
        }
    }

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start,0);

    // allocate vectors in device memory

    cudaMalloc((void**)&d_type, atom_int_size);
    cudaMalloc((void**)&d_m, atom_float_size);
    cudaMalloc((void**)&d_x, atom_double3_size);
    cudaMalloc((void**)&d_v, atom_double3_size);
    cudaMalloc((void**)&d_f, atom_double3_size);
    cudaMalloc((void**)&d_a, atom_double3_size);
    cudaMalloc((void**)&d_stress, atom_float6_size);

    // cudaMalloc((void**)&d_a_old, atom_double3_size);
    // cudaMalloc((void**)&d_a_new, atom_double3_size);

    cudaMalloc((void**)&d_k, bond_float_size);
    cudaMalloc((void**)&d_r0, bond_float_size);
    cudaMalloc((void**)&d_rc, bond_float_size);
    cudaMalloc((void**)&d_atom_i, bond_int_size);
    cudaMalloc((void**)&d_atom_j, bond_int_size);

    // copy vectors from host memory to device memory

    cudaMemcpy(d_type, sys.type, atom_int_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_m, h_m, atom_float_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, sys.x, atom_double3_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_v, sys.v, atom_double3_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_f, 0, atom_double3_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_a, sys.a, atom_double3_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_stress, h_stress, atom_float6_size, cudaMemcpyHostToDevice);
    
    // cudaMemcpy(d_a_old, sys.a, atom_double3_size, cudaMemcpyHostToDevice);
    // cudaMemcpy(d_a_new, sys.a, atom_double3_size, cudaMemcpyHostToDevice);

    cudaMemcpy(d_k, h_k, bond_float_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_r0, h_r0, bond_float_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rc, h_rc, bond_float_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_atom_i, sys.atom_i, bond_int_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_atom_j, sys.atom_j, bond_int_size, cudaMemcpyHostToDevice);

    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);

    float cpytime;
    cudaEventElapsedTime( &cpytime, start, stop);
    printf("Copy data from host to device... %f (ms) \n", cpytime);

    // volatile bool flag = true;
    for (int ti = 0; ti < timesteps; ti++) {

        cudaEventRecord(start,0);

        for (std::vector<Fix>::iterator ifix = fixes.begin(); ifix != fixes.end(); ifix++) {
            fix_rigidmove<<<dimGridforAtom, dimBlockforAtom>>>(
                d_type, d_x, d_v, d_a,
                ifix->type, ifix->dispx, ifix->dispy, ifix->dispz,
                dt, sys.natoms*3);
        }
        
        cudaDeviceSynchronize();

        // calculate velocity at the next half timestep v(t+dt/2) = v(t) + 1/2*a(t)*dt
        verlet_update_vel<<<dimGridforAtom, dimBlockforAtom>>>(d_v, d_a, dt/2.0, sys.natoms*3);

        cudaDeviceSynchronize();

        // calculate position at the next timestep x(t+dt) = x(t) + v(t+dt/2)*dt 
        verlet_update_pos<<<dimGridforAtom, dimBlockforAtom>>>(d_x, d_v, dt, sys.natoms*3);

        cudaDeviceSynchronize();

        // calculate accerlation at the next timestep a(t+dt) by x(t+dt)
        calculate_force<<<dimGridforBond, dimBlockforBond>>>(
            d_x, d_f,
            d_k, d_r0, d_rc, d_atom_i, d_atom_j, 
            sys.natoms, sys.nbonds
        );

        cudaDeviceSynchronize();

        verlet_update_acc<<<dimGridforAtom, dimBlockforAtom>>>(d_a, d_m, d_f, sys.natoms);

        cudaDeviceSynchronize();

        // calculate velocity at the next full timestep v(t+dt) = v(t+dt/2) + 1/2*a(t+dt)*dt
        verlet_update_vel<<<dimGridforAtom, dimBlockforAtom>>>(d_v, d_a, dt/2.0, sys.natoms*3);

        cudaDeviceSynchronize();

        // TODO: dump setting
        if (ti % dump.ndump == 0 || ti % thermo.nthermo == 0) {
            cudaMemcpy(sys.x, d_x, atom_double3_size, cudaMemcpyDeviceToHost);
            cudaMemcpy(sys.v, d_v, atom_double3_size, cudaMemcpyDeviceToHost);
            cudaMemcpy(h_k, d_k, bond_float_size, cudaMemcpyDeviceToHost);

            if (ti % dump.ndump == 0) {

                stress_virial_kin<<<dimGridforAtom6, dimBlockforAtom>>>(
                    d_stress, d_m, d_v, 
                    sys.natoms);

                cudaDeviceSynchronize();

                stress_virial_pot<<<dimGridforBond, dimBlockforBond>>>(
                    d_stress, d_m, d_x, d_v,
                    d_k, d_r0, d_rc,
                    d_atom_i, d_atom_j,
                    sys.nbonds
                );

                cudaDeviceSynchronize();

                cudaMemcpy(h_stress, d_stress, atom_float6_size, cudaMemcpyDeviceToHost);

                dump.write_lmpdump(ti, sys);
            }
            if (ti % thermo.nthermo == 0) thermo.write_thermo(ti, sys);
        }

        cudaEventRecord(stop,0);
        cudaEventSynchronize(stop);
    }
    
    
    // Release device memory
    cudaFree(d_type);
    cudaFree(d_m);
    cudaFree(d_x);
    cudaFree(d_v);
    cudaFree(d_a);
    cudaFree(d_stress);
    // cudaFree(d_a_old);
    // cudaFree(d_a_new);

    cudaFree(d_k);
    cudaFree(d_r0);
    cudaFree(d_rc);
    cudaFree(d_atom_i);
    cudaFree(d_atom_j);
 
    // Release host memory
    free(h_m);
    free(h_k);
    free(h_r0);
    free(h_rc);
    free(h_stress);

    cudaDeviceReset();    
    return 0;
}




