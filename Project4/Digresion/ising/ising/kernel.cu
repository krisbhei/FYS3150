
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <iostream>

#include <curand.h>
#include <curand_kernel.h>

#include <stdlib.h>
#include <random>

#include "../../codes/common/cpu_anim.h"
#include "../../codes/common/book.h"


//This code is highly inspired by the codes found in 
//chapter 8 from "CUDA by example - an introduction to general purpose GPU programming"
//written by Jason Sanders and Edward Kandrot.
//The necessary headers, cpu_anim.h and book.h are used in the book and therefore used here as well.

texture<float, 2>  texIn;

struct DataBlock //To store necessary data
{
	int dim;
	double temp;
	unsigned char   *output_bitmap;
	float           *dev_inSrc;
	CPUAnimBitmap  *bitmap;

};

__global__ void kernel_spins(float *dst, int rand_x, int rand_y, double rand_energy,double temp) 
{
	//Performing a MC-cycle at one thread
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	if (x == rand_x && y == rand_y)
	{
		double  spin_top, spin_left, spin_this, spin_right, spin_bottom;
		double w[17];
		for (int i = 0; i < 17; i++) w[i] = 0;
		for (int i = -8; i < 9; i += 4) w[i + 8] = exp(-((double)i) / temp);

		
		spin_top = tex2D(texIn, x, y - 1);
		spin_left = tex2D(texIn, x - 1, y);
		spin_this = tex2D(texIn, x, y);
		spin_right = tex2D(texIn, x + 1, y);
		spin_bottom = tex2D(texIn, x, y + 1);
		

		int deltaEnergy = 2 * spin_this*
			(spin_left +
			spin_right +
			spin_top +
			spin_bottom);
		if (rand_energy <= w[deltaEnergy+8])
		{
			int offset = x + y * blockDim.x * gridDim.x;
			dst[offset] = spin_this == 1 ? 0 : 1;

		}
		
	}
	
} //End: kernel_spins


void anim_gpu(DataBlock *d,int ticks) 
{
	int numSpins = d->dim;
	dim3    blocks(numSpins / 16, numSpins / 16);
	dim3    threads(16, 16);
	CPUAnimBitmap  *bitmap = d->bitmap;

	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_real_distribution<double> distr(0.0, 1.0);

	int cycle = 2000 + ticks * 2000;
	double temp = d->temp;
	for (int i = 0; i < 2000; i++)
	{
		float * out = d->dev_inSrc;
		int rand_x = (int)(distr(gen)*(double)numSpins);
		int rand_y = (int)(distr(gen)*(double)numSpins);
		double rand_energy = distr(gen);
		kernel_spins << <blocks, threads >> >(out, rand_x, rand_y, rand_energy, temp);
	}
	
	float_to_color << <blocks, threads >> >(d->output_bitmap,
		d->dev_inSrc);

	//Render the results after the MC-cycles
	HANDLE_ERROR(cudaMemcpy(bitmap->get_ptr(),
		d->output_bitmap,
		bitmap->image_size(),
		cudaMemcpyDeviceToHost));

	std::cout << "Cycle: " << cycle << std::endl;

} //End: anim_gpu

void anim_exit(DataBlock *d) 
{
	cudaUnbindTexture(texIn);
	HANDLE_ERROR(cudaFree(d->dev_inSrc));
} //End: anim_exit


int main(void) 
{
	DataBlock   data;

	int numSpins;
	std::cout << "Number of spins in each dimension: ";
	std::cin >> numSpins;
	std::cout << "Temperature: ";
	double temp;
	std::cin >> temp;

	CPUAnimBitmap bitmap(numSpins, numSpins, &data);
	data.temp = temp;
	data.dim = numSpins;
	data.bitmap = &bitmap;

	int imageSize = bitmap.image_size();

	HANDLE_ERROR(cudaMalloc((void**)&data.output_bitmap,
		imageSize));

	HANDLE_ERROR(cudaMalloc((void**)&data.dev_inSrc,
		imageSize));

	cudaChannelFormatDesc desc = cudaCreateChannelDesc<float>();

	HANDLE_ERROR(cudaBindTexture2D(NULL, texIn,
		data.dev_inSrc,
		desc, numSpins, numSpins,
		sizeof(float) * numSpins));


	//Initilizing at a ground state
	float *spins = (float*)malloc(imageSize);

	for (int i = 0; i<numSpins*numSpins; i++)
	{
		spins[i] = 1;
	}

	//Send over the ground state to the GPU
	HANDLE_ERROR(cudaMemcpy(data.dev_inSrc, spins,
		imageSize,
		cudaMemcpyHostToDevice));
	free(spins);

	//Render and show
	bitmap.anim_and_exit((void(*)(void*, int))anim_gpu,
		(void(*)(void*))anim_exit);
	return 0;

}//End: main