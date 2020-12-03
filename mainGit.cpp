#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "Wave.h"

using namespace std;

#define DO3  261.63
#define RE3  293.66
#define MI3  329.63
#define FA3  349.23
#define SOL3 392.00
#define LA3  440.00
#define SI3  493.88

#define NB_NOTES 7

double NOTES[NB_NOTES] = {DO3, RE3, MI3, FA3, SOL3, LA3, SI3};

void add_note(vector<double>& signal, double freq_note, double duration, double freq_ech, double amplitude = 1.0);

// signal in [-1,-1] -> [0, 255]
unsigned char normalize(double signal_ech);

vector<unsigned char> to_data8(vector<double>& signal);

void write_signal(string& filepath, vector<double>& signal, int nb_channels, int sampling_freq);

vector<double> read_signal(string filepath, int nb_channels, int sampling_freq);

/*
	This FFT has been proposed by Paul Bourke
	http://paulbourke.net/miscellaneous/dft/
	This computes an in-place complex-to-complex FFT
	x and y are the real and imaginary arrays of 2^m points.
	dir =  1 gives forward transform
	dir = -1 gives reverse transform
	You MUST compute first the value m such that
	2^(m-1) < n (size of your signal) <= 2^m
	allocate a new signal of nm=2^m values
	then fill the n first values of this new signal
	with your signal and fill the rest with 0
	WARNING : you must pass m, not nm !!!
*/
int DFT(int dir,int m,double *x1,double *y1);

vector<double> DFT_visualize(vector<double>& signal_real, vector<double>& signal_imaginary);

void normalisation(double* signal, int N);


/* Renvoie m tel que 2^m >= n  */
int next_pow2(int n) {

    int next_pow2 = 1;
    int m = 0;

    while (next_pow2 < n) {
        next_pow2 *= 2;
        m++;
    }

    return m;
}

void add_gamme(vector<double>& signal, vector<double>& gamme, int sampling_freq){
    double max_freq = sampling_freq / 2.0;
    double max_index = signal.size() / 2.0;

    for (size_t i = 0; i < gamme.size(); ++i){
        size_t index = (max_index * gamme[i]) / max_freq;
        size_t index_sym = max_index + (max_index - index);
        signal[index] = 1;
        signal[index_sym] = 1;
    }
}

void filter_low_pass(vector<double>& signal, double);


int main(int argc, const char *argv[]) {

    cout << "Initialisation du signal...\n";

    double sample_time = 2; // en secondes
    string filename = "sounds/la.wav";
    string filename_FFT = "sounds/la_FFT.wav";
    int sampling_freq = 44100;
    int nb_channels = 1;

    vector<double> signal = read_signal("sounds/La.wav", nb_channels, sampling_freq);

    //add_note(signal, 0, 5, sampling_freq);

    /* les signaux doivent avoir une taille en puissance de deux pour que la FFT
       puisse fonctionnée convenablement */

    int FFT_m = next_pow2(signal.size());
    long long FFT_size = pow(2, FFT_m);

    cout << FFT_m << " " << FFT_size << "\n";

    vector<double> signal_real(FFT_size);
    vector<double> signal_imaginary(FFT_size);

    copy(signal.begin(), signal.end(), signal_real.begin());

    cout << "Transformation de fourrier en cours..." << endl;

    DFT(1, FFT_m, signal_real.data(), signal_imaginary.data());

    vector<double> gamme;

    //gamme.push_back(DO3);
    //gamme.push_back(MI3);
    //gamme.push_back(SOL3);

    //add_gamme(signal_real, gamme, sampling_freq);
    vector<double> visuel =DFT_visualize(signal_real,signal_imaginary);
    write_signal(filename_FFT, signal_real, nb_channels, sampling_freq);

    DFT(-1, FFT_m, signal_real.data(), signal_imaginary.data());

    cout << "Ecriture du signal dans " << filename << "\n";

    copy(signal_real.begin(), signal_real.begin() + signal.size(), signal.begin());

    normalisation(signal.data(), signal.size());

    write_signal(filename, , nb_channels, sampling_freq);

    return 0;
}

/*
S(t) = sin(2 * PI * freq * t)
S(n) = S(n * tau) = sin(2 * Pi * freq * n * tau) = sin(2 * PI * freq * 1 / (freq * tau) * n
*/

void add_note(vector<double>& signal, double freq_note, double duration, double freq_ech, double amplitude){

    long long nb_ech = (long long)(duration * freq_ech);

    for (long long i = 0; i < nb_ech; ++i)
        signal.push_back(amplitude * sin(2 * M_PI * freq_note * (i / freq_ech)));
}

// signal in [-1,-1] -> [0, 255]
unsigned char normalize(double signal_ech){
    return (signal_ech + 1.0) * 127.5;
}

vector<unsigned char> to_data8(vector<double>& signal) {

    vector<unsigned char> data8(signal.size());

    for (size_t i = 0; i < signal.size(); ++i)
        data8[i] = normalize(signal[i]);

    return move(data8);
}




vector<double> read_signal(string filepath, int nb_channels, int sampling_freq) {

    unsigned char* data8 = NULL;
    int size = 0;
    Wave wave(data8, size, nb_channels, sampling_freq);
    wave.read((char*)filepath.c_str());

    wave.getData8(&data8, &size);

    vector<double> signal(size);

    for (int i = 0; i < size; ++i)
        signal[i] = ((data8[i] * 2.0) / 255.0) - 1;

    delete[] data8;

    return move(signal);
}

void normalisation(double* signal, int N) {

    double max, min;
    double *pt, *fin;

    pt = signal;
    fin = pt + N;
    min = (*pt);
    max = (*pt);

    while (pt < fin) {
        min = (*pt) < min ? (*pt) : min;
        max = (*pt) > max ? (*pt) : max;
        pt++;
    }

    if(max > min + 1e-8)
        max = 2.0 / (max - min);
    else
        max = 1e-8;

    pt = signal;

    while (pt < fin) {
        (*pt) = (((*pt) - min) * max) - 1.0;
        pt++;
    }
}

vector<double> DFT_visualize(vector<double>& signal_real, vector<double>& signal_imaginary) {

    double max = 0;
    vector<double> norme(signal_real.size());

    for (size_t i = 0; i < signal_real.size(); ++i) {
        norme[i] = sqrt(signal_real[i] * signal_real[i] + signal_imaginary[i] * signal_imaginary[i]);
        max = max < norme[i] ? norme[i] : max;
    }

    /* normalize val -> [-1, 1] */

    for (size_t i = 0; i < signal_real.size(); ++i) {
        norme[i] = (2.0 * norme[i] / max) - 1.0;
    }

    return move(norme);
}

int DFT(int dir,int m,double *x1,double *y1)
{
    long i,k;
    double arg;
    double cosarg,sinarg;
    double *x2=NULL,*y2=NULL;

    x2 = static_cast<double *>(malloc(m * sizeof(double)));
    y2 = static_cast<double *>(malloc(m * sizeof(double)));
    if (x2 == NULL || y2 == NULL)
        return(false);

    for (i=0;i<m;i++) {
        x2[i] = 0;
        y2[i] = 0;
        arg = - dir * 2.0 * 3.141592654 * (double)i / (double)m;
        for (k=0;k<m;k++) {
            cosarg = cos(k * arg);
            sinarg = sin(k * arg);
            x2[i] += (x1[k] * cosarg - y1[k] * sinarg);
            y2[i] += (x1[k] * sinarg + y1[k] * cosarg);
        }
    }

    /* Copy the data back */
    if (dir == 1) {
        for (i=0;i<m;i++) {
            x1[i] = x2[i] / (double)m;
            y1[i] = y2[i] / (double)m;
        }
    } else {
        for (i=0;i<m;i++) {
            x1[i] = x2[i];
            y1[i] = y2[i];
        }
    }

    free(x2);
    free(y2);
    return(true);
}