//Compile with -lm

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <complex.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "vor.h"

#define PI 3.14159265359
#define NUM_THREADS	4

//Global Vars
char *filepath = "./VOR-diadr-0.raw";	//File Source
int Fs = 262144;  				//Sample rate
int entries = 262144; 			//Complex entries
int fir_size = 45;
int *am_max_peaks, *fm_max_peaks, *am_min_peaks, *fm_min_peaks;
float *am, *fm;

void *PeakThread(void *threadid)
{
	long tid;
	tid = (long)threadid;

	if (tid == 0) am_max_peaks = find_max_peaks(am,4000);
	else if (tid == 1) fm_max_peaks = find_max_peaks(fm,4000);
	else if (tid == 2) am_min_peaks = find_min_peaks(am,4000);
	else fm_min_peaks = find_min_peaks(fm,4000);

	pthread_exit(NULL);
}

//MAIN
int main(){

	struct timespec start, finish;
	double elapsed;
	clock_gettime(CLOCK_MONOTONIC, &start);

	system("clear");

	float complex* z;
	float *z_real, *z_imag;	

	float fm_period_samp, am_period_samp;
	float fm_freq, am_freq;
	float am_mod_dpth ,fm_mod_dpth;

	int phase;
	float phase_deg;

	z_real = (float*) malloc(entries * sizeof(float));
	z_imag = (float*) malloc(entries * sizeof(float));

	//while (1){
	
		system("rm plot.dat; touch plot.dat");
		system("rm am.dat; touch am.dat");
		system("rm fm.dat; touch fm.dat");
		//system("rm capture.dat; rtl_sdr -f 110800000 -s 250000 -n 250000 capture.dat");

		z = getData();

		for (int i=0; i<entries; i++){
			z_real[i] = creal(z[i]);
			z_imag[i] = cimag(z[i]);
		}	

		plot(z_real, entries, "plot.dat");

		am = amdemod(z_real, z_imag);
		plot(am, entries, "am.dat");

		fm = fmdemod(z_real, z_imag);
		plot(fm, entries, "fm.dat");

		pthread_t threads[NUM_THREADS];
	   	int rc;
	   	long t;

		for(t=0;t<NUM_THREADS;t++){
		 	rc = pthread_create(&threads[t], NULL, PeakThread, (void *)t);
		 	if (rc){
		   		printf("ERROR; return code from pthread_create() is %d\n", rc);
		   		exit(-1);
		   	}
		}
		for (t=0; t<NUM_THREADS; t++)
       		pthread_join(threads[t], NULL);

		am_period_samp = find_freq(am_max_peaks);
		fm_period_samp = find_freq(fm_max_peaks);
		am_freq = Fs/find_freq(am_max_peaks);
		fm_freq = Fs/find_freq(fm_max_peaks);
		am_mod_dpth = find_am_mod_depth(am, am_max_peaks, am_min_peaks);
		fm_mod_dpth = find_am_mod_depth(fm, fm_max_peaks, fm_min_peaks);

		printf("AM Period in Samples: %f \t FM Period in Samples: %f\n", 
			am_period_samp, fm_period_samp);
		printf("\nAM Frequency: %fHz \t\t FM Frequency: %fHz\n", 
			am_freq, fm_freq);
		printf("\nAM Modulation Depth: %f%% \t FM Modulation Depth: %f%%\n", 
			am_mod_dpth, fm_mod_dpth);

		phase = find_phase(am_max_peaks, fm_max_peaks);
		phase_deg = ((float)phase/(float)find_freq(am_max_peaks))*360.0;
		printf("\nPhase diff: %d samples or %f deg\n", phase, phase_deg);

		printf("\n");
	//}

	clock_gettime(CLOCK_MONOTONIC, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	printf("Program execution time: %g seconds\n", elapsed);

	pthread_exit(NULL);
	return 0;

}

//Read data from file
float complex* getData(){
/*
	unsigned char data[entries*2];
	float real[entries], imag[entries];
	float complex *com;

	com = (float complex*) malloc(entries * sizeof(float complex));

	FILE *fp = fopen(filepath, "r");
	if (fp == NULL){
		perror("fopen failed error");
		return NULL;
	}

	fread(data, sizeof(unsigned char), entries*2, fp);
	
	for (int i=0; i<entries; i++){
		real[i] = 0.008 * (((float)data[i*2]) - 127.0);         //*10?     -127?
		imag[i] = 0.008 * (((float)data[i*2+1]) - 127.0);
		com[i] = real[i] + imag[i] * I;
	}

	fclose(fp);
	return com;
*/
	float data[entries*2];
	float complex *com;

	com = (float complex*) malloc(entries * sizeof(float complex));

	FILE *fp = fopen(filepath, "r");
	if (fp == NULL){
		perror("fopen failed error");
		return NULL;
	}

	fread(data, sizeof(float), entries*2, fp);
	
	for (int i=0; i<entries; i++){
		com[i] = data[i*2] + data[i*2+1] * I;
	}

	fclose(fp);
	return com;

}

//Plot a signal - Make file for GNUPlot
//~$ gnuplot
//~$ plot "plot.dat"
void plot(float* vector, int len, char* filename){
	
	FILE *plot;
	plot = fopen(filename, "rb+");
	for (int i=0; i<len; i++) fprintf(plot, "%d %g\n", i, vector[i]);

	fflush(plot);

}

/*get filter from files   

<url: www.arc.id.au/FilterDesign.html>
Fb = 30kHz
M(odd) length = 45
Fs = 262144Hz
Att 45dB
1- lpf3k.filter - low pass filter 3kHz
2- hpf5k.filter - high pass filter 5kHz
3- lpf30k.filter - low pass filter 30kHz
	
*/	
										
float* getFilter(int type){

	char* filterpath;
	int size = fir_size;
	float* data;

	data = (float*) malloc(size * sizeof(float));

	if (type == 1)
		filterpath = "./lpf3k.filter";
	else if (type == 2)
		filterpath = "./hpf5k.filter";
	else
		filterpath = "./lpf30k.filter";

	FILE *fp = fopen(filterpath, "r");
	if (fp == NULL){
		perror("fopen failed error");
		return NULL;
	}

	for (int i=0; i<size; i++){
		float temp;
		fscanf(fp, "%f", &temp);
		data[i] = temp;
	}

	fclose(fp);
	return data;

}

float* convolve(float* X, float* h){

	float* Y;
	Y = (float*) malloc(entries*sizeof(float));
	float sum;
	int N = fir_size;
	
	for (int n=0; n<entries; n++){
		sum = 0;

		for (int k=0; k<N-1; k++){
			if (n>=k) sum += h[k] * X[n-k];
		}

		Y[n] = sum;
	}

	return Y;
}

//AM Demodulation
float* amdemod(float* z_real, float* z_imag){

	//apply lpf30k
	float *filtered_r = lpf(z_real, 30000);
	float *filtered_i = lpf(z_imag, 30000);

	//complex to mag^2
	float* mag2;
	float* am;
	mag2 = (float*) malloc(entries*sizeof(float));
	for (int i=0; i<entries; i++){
		mag2[i] = pow(z_real[i],2) + pow(z_imag[i],2);
		//multiply x60
		mag2[i] = 60 * mag2[i];
	}

	am = moving_avg(2000, mag2);
	
	return am;
}

float* moving_avg(int frameSize, float* data){

	float sum=0;
	float* avgPoints;
	avgPoints = (float*) malloc(entries*sizeof(float));
	avgPoints[0]=0;
	data[0] = 0;

	for (int i=1; i<frameSize; i++){
		avgPoints[i] = avgPoints[i-1] + data[i];
	}

	for (int i=frameSize; i<entries; i++){
		avgPoints[i] = avgPoints[i-1] + data[i] - data[i - frameSize];
	}
	
	return avgPoints;
} 
 /*
int* find_max_peaks(float* am, int win_size){

	int* peaks;
	int k=0;
	peaks = (int*) malloc(entries/1000*sizeof(int));

	for(int i=0; i<entries-win_size; i++){
		float max = 0;
		for(int j=0; j<win_size; j++){
			if (am[j+i]>max){
				max = am[j+i];
			}
		}
    	if (max == am[i+win_size/2]){
			peaks[k] = i+win_size/2;
			k++;
		}
    }

	peaks[k]=0;

	return peaks;
			
}
*/

int* find_max_peaks(float* am, int win_size){

	int* peaks;
	int k=0;
	int prev_max_pos=0;
	peaks = (int*) malloc(entries/1000*sizeof(int));

	for(int i=1; i<entries-win_size; i++){
		float max = 0;
		if (i==1 || prev_max_pos<i){
			for(int j=0; j<win_size; j++){
				if (am[j+i]>max){
					max = am[j+i];
					prev_max_pos = j+i;
				}
			}
		}
		else {
			for(int j=prev_max_pos-i; j<win_size; j++){
				if (am[j+i]>max){
					max = am[j+i];
					prev_max_pos = j+i;
				}
			}
		}

    	if (max == am[i+win_size/2]){
			peaks[k] = i+win_size/2;
			k++;
		}
    }

	peaks[k]=0;

	return peaks;
			
}

int* find_min_peaks(float* am, int win_size){

	int* peaks;
	int k=0;
	peaks = (int*) malloc(entries/1000*sizeof(int));

	for(int i=0; i<entries-win_size; i++){
		float min = 3.4028e+38;
		for(int j=0; j<win_size; j++){
			if (am[j+i]<min){
				min = am[j+i];
			}
		}
    	if (min == am[i+win_size/2]){
			peaks[k] = i+win_size/2;
			k++;
			i=i+win_size/2-1;
		}
    }

	peaks[k]=0;

	return peaks;
			
}

float find_freq(int *peaks){

	int k = 0,sum = 0;
	for(int i=1; peaks[i]>0; i++){
		sum += peaks[i]-peaks[i-1];
        k++;
	}
	return sum/(float)k;

}

float find_am_mod_depth(float* am, int* max_peaks, int* min_peaks){

	float sum = 0.0;
	int k=0;

	for(int i=0; max_peaks[i]>0 && min_peaks[i]>0; i++){
		//printf("%d: pos %d val %f\n", i, peaks[i], am[peaks[i]]);
		//printf("%d: pos %d val %f\n", i, min_peaks[i], am[min_peaks[i]]);
		sum += ( am[max_peaks[i]] - am[min_peaks[i]] ) / ( am[max_peaks[i]] + am[min_peaks[i]] );
		k++;
	}

	return (100*sum/k);

}

//FM Demodulation
float* fmdemod(float* z_real, float* z_imag){

	float *fm_real, *fm_imag, *fm;
	fm_real = (float*) malloc(entries*sizeof(float));
	fm_imag = (float*) malloc(entries*sizeof(float));
	fm = (float*) malloc(entries*sizeof(float));

	//apply hpf5k
	float *filtered_r = hpf(z_real, 5000);
	float *filtered_i = hpf(z_imag, 5000);

	//make complex cosine signal
	float *cos_real, *cos_imag;
	cos_real = (float*) malloc(entries*sizeof(float));
	cos_imag = (float*) malloc(entries*sizeof(float));

	//multiply
	for (int i=0; i<entries; i++){
		cos_real[i] = cos(2*PI*9960*i/Fs);
		cos_imag[i] = -sin(2*PI*9960*i/Fs);
		fm_real[i] = filtered_r[i]*cos_real[i] - filtered_i[i]*cos_imag[i];
		fm_imag[i] = filtered_r[i]*cos_imag[i] + filtered_i[i]*cos_real[i];
	}
	
	//apply lpf3k
	filtered_r = lpf(fm_real, 3000);
	filtered_i = lpf(fm_imag, 3000);
	
	//split
	float *fm_real_delay, *fm_imag_delay;
	fm_real_delay = (float*) malloc(entries*sizeof(float));
	fm_imag_delay = (float*) malloc(entries*sizeof(float));

	//delay
	for (int i=1; i<entries; i++){
		fm_real_delay[i] = filtered_r[i-1];
		fm_imag_delay[i] = filtered_i[i-1];
	}

	//conjurate
	for (int i=0; i<entries; i++) filtered_i[i] = -filtered_i[i];

	//multiply
	for (int i=0; i<entries; i++){
		fm_real[i] = filtered_r[i]*fm_real_delay[i] - filtered_i[i]*fm_imag_delay[i];
		fm_imag[i] = filtered_r[i]*fm_imag_delay[i] + filtered_i[i]*fm_real_delay[i];
	}

	//complex to float + devide
	for (int i=0; i<entries; i++){
		fm[i] = atan(fm_imag[i] / fm_real[i]);		
	}

	//moving average
	float* fm_final;
	fm_final = (float*)malloc(entries*sizeof(float));
	fm_final = moving_avg(4000, fm);
	
	return fm_final;

}

int find_phase(int* am_peaks, int* fm_peaks){

	int counter=0, sum=0;

	for (int i=1; i<31; i++){
		if (fm_peaks[i] <= am_peaks[i+1] && fm_peaks[i] >= am_peaks[i]){
			sum += fm_peaks[i] - am_peaks[i];
			counter++;	
		}
	}

	return sum/counter;

}

float* hpf(float* signal, int CUTOFF){

    float RC = 1.0/(CUTOFF*2*3.14);
    float dt = 1.0/Fs;
    float alpha = RC/(RC + dt);
    float *filteredArray;
	filteredArray = (float*)malloc(entries*sizeof(float));

    filteredArray[0] = signal[0];
    for (int i = 1; i<entries; i++){
        filteredArray[i] = alpha * (filteredArray[i-1] + signal[i] - signal[i-1]);
    }
    return filteredArray;

}

float* lpf(float* signal, int CUTOFF){

	float RC = 1.0/(CUTOFF*2*3.14);
    float dt = 1.0/Fs;
    float alpha = dt/(RC+dt);
    float *filteredArray;
	filteredArray = (float*)malloc(entries*sizeof(float));

    filteredArray[0] = signal[0];
    for(int i=1; i<entries; i++){
        filteredArray[i] = filteredArray[i-1] + (alpha*(signal[i] - filteredArray[i-1]));
    }
    return filteredArray;

}
