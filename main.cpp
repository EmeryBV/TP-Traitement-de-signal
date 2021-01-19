#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <cstring>
#include <vector>
#include "Wave.h"

using namespace std;
#define PI  3.14159265358979323846
#define DO 261.63
#define RE 293.66
#define MI 329.63
#define FA 349.23
#define SOL 392.00
#define LA 440
#define SI 493.88
#define sampling_freq 44100
#define duree 0.3

void DFT(double *signal, double *partie_reelle, double *partie_imaginaire, int N);
void IDFT(double * signal, double * partie_reelle, double * partie_imaginaire, int N);
int FFT(int dir,int m,double *partie_reelle,double *partie_imaginaire);
void normalisation(double* signal, int N);
void filtre_passe_haut(double *signal, int N, double frequence, int i);
unsigned char normaliser_double_to_unsiCharn(double signal) {
    return (signal +  1.0) * 127.5;
}

double normaliser_unsiChar_to_double(unsigned char signal) {
    return (signal / 127.5) - 1;
}

float rapport(float sampling_freg, float test_frequence) {
    return (1 / test_frequence) / (1 / sampling_freg);
}
template<int N>
vector<double> creerMusique(vector<double> &partition, unsigned char (&data8)[N]) {
    cout << partition.size() << "\n";
    vector<double> signal;
    for (int n = 0; n < partition.size(); n++) {
        float rapportNote = rapport(sampling_freq, (float) partition[n]);
        for (int i = 0; i < sampling_freq * duree / partition.size(); i++) {
            signal.push_back(normaliser_double_to_unsiCharn(sin((float) (i) / (rapportNote / (2 * PI)))));
            data8[i + n * sizeof(data8) / partition.size()] = normaliser_double_to_unsiCharn(
                    sin((float) (i) / (rapportNote / (2 * PI))));
        }
    }
    return signal;
}

void DFT_visualize(vector<double> signal_real, vector<double> signal_imaginary,vector<double> &norme) {

    double max = 0;


    for (size_t i = 0; i < signal_real.size(); ++i) {
        norme[i] = sqrt(signal_real[i] * signal_real[i] + signal_imaginary[i] * signal_imaginary[i]);
        max = max < norme[i] ? norme[i] : max;
    }

    /* normalize val -> [-1, 1] */

    for (size_t i = 0; i < signal_real.size(); ++i) {
        norme[i] = (2.0 * norme[i] / max) - 1.0;
    }

}

vector<unsigned char> to_data8(vector<double> &signal) {

    vector<unsigned char> data8(signal.size());

    for (size_t i = 0; i < signal.size(); ++i)
        data8[i] = normaliser_double_to_unsiCharn(signal[i]);

    return move(data8);
}

vector<double> read_signal(string filepath, int nb_channels) {

    unsigned char *data8 = NULL;
    int size = 0;
    Wave wave(data8, size, nb_channels, sampling_freq);
    wave.read((char *) filepath.c_str());

    wave.getData8(&data8, &size);

    vector<double> signal(size);

    for (int i = 0; i < size; ++i)
        signal[i] = ((data8[i] * 2.0) / 255.0) - 1;

    delete[] data8;

    return move(signal);
}

void write_signal(string filename, vector<double> &signal, int nb_channels) {

    vector<unsigned char> data8(to_data8(signal));

    Wave wave(data8.data(), data8.size(), nb_channels, sampling_freq);

// convertie en char* car la fonction write de modifie pas la constante de toute façon
    wave.write((char *) filename.c_str());
}

/* Renvoie m tel que 2^m >= n  */
int puissance2(int n) {

    int puissance= 1;
    int m = 0;

    while (puissance< n) {
        puissance *= 2;
        m++;
    }

    return m;
}


int main() {
    // en secondes
    char *filename = "sounds/La.wav";
    char *filenameDFT = "sounds/LaDFT.wav";
    char *filenameIDFT = "sounds/LaIDFT.wav";
    char *filenameFFT = "sounds/LaFFT.wav";
    char *filenameIFFT = "sounds/LaIFFT.wav";
    short nb_channels = 1;

    vector<double> partition;
//   partition.push_back(DO);
   partition.push_back(LA);
//   partition.push_back(MI);
//   partition.push_back(FA);
//   partition.push_back(SOL);
//   partition.push_back(LA);
//   partition.push_back(SI);


    unsigned char data8[static_cast<int>(sampling_freq * duree)];
    creerMusique(partition, data8);
    Wave wave = Wave(data8, sampling_freq * duree, nb_channels, sampling_freq);
    wave.write(filename);
    vector<double> signal = read_signal(filename, 1);


    int nbr_puissance = puissance2(signal.size());
    long long longueurMin = pow(2, nbr_puissance);

    vector<double> partie_reelle(longueurMin);
    vector<double> partie_imaginaire(longueurMin);
    copy(signal.begin(), signal.end(), partie_reelle.begin());
    cout << "Transformation de fourrier en cours..." << endl;


    //DFT
    copy(signal.begin(), signal.end(), partie_reelle.begin());
    DFT(signal.data(), partie_reelle.data(), partie_imaginaire.data(), signal.size());
    copy(partie_reelle.begin(), partie_reelle.begin() + signal.size(), signal.begin());
    vector<double> norme(partie_reelle.size());

    DFT_visualize(partie_reelle,partie_imaginaire, norme);

    write_signal(filenameDFT,norme, 1);

    //IDFT
    IDFT(norme.data(),partie_reelle.data(),partie_imaginaire.data(),signal.size());
    copy(partie_reelle.begin(), partie_reelle.begin() + signal.size(), signal.begin());
    write_signal(filenameIDFT,norme, 1);



//    FFT(1, nbr_puissance, partie_reelle.data(),partie_imaginaire.data());
//    vector<double> norme2(partie_reelle.size());
////    DFT_visualize(partie_reelle,partie_imaginaire, norme2);
//    write_signal(filenameFFT,signal, 1);
//
//
//    //IFFT
//    FFT(-1,nbr_puissance,partie_reelle.data(),partie_imaginaire.data());
//
//    DFT_visualize(partie_reelle, partie_imaginaire,norme2);
//    write_signal(filenameIFFT,norme2, 1);


    return 0;

}
void filtre_passe_haut(double *signal, int N, double frequence, int i) {
    for (int i = 0; i < N; ++i) {
    if (signal[i]<frequence){
        signal[i]= 0;
        }
    }
}
void DFT(double *signal, double *partie_reelle, double *partie_imaginaire, int N) {
    int n, k;
    double deux_pi_N = 2.0 * PI / (double) N;
    double *pti, *ptr, *pt;
    double omega, theta;
    ptr = partie_reelle;
    pti = partie_imaginaire;

    for (k = 0; k < N; k++) {
        omega = deux_pi_N * (double) k;
        (*ptr) = 0.0;
        (*pti) = 0.0;
        pt = signal;
        for (n = 0; n < N; n++) {
            theta = omega * (double) n;
            (*ptr) += (*pt) * cos(theta);
            (*pti) -= (*pt++) * sin(theta);
        }
        ptr++;
        pti++;
    }
}

void IDFT(double * signal, double * partie_reelle, double * partie_imaginaire, int N){
    for(int i=0;i<N;i++){
        double sum = 0;
        for(int j=0;j<N;j++){
            sum += partie_reelle[j]*cos(2*PI*j*i/N) - partie_imaginaire[j]*sin(2*PI*j*i/N);
        }
        signal[i] = sum/N;
    }
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


/*
	This FFT has been proposed by Paul Bourke
	http://paulbourke.net/miscellaneous/dft/
	This computes an in-place complex-to-complex FFT
	x and partie_imaginaire are the real and imaginary arrays of 2^m points.
	dir =  1 gives forward transform
	dir = -1 gives reverse transform
	You MUST compute first the value m such that
	2^(m-1) < n (size of your signal) <= 2^m
	allocate a new signal of nm=2^m values
	then fill the n first values of this new signal
 with your signal and fill the rest with 0
	WARNING : you must pass m, not nm !!!
	*/

int FFT(int dir, int m, double *x, double *y) {

    int n,i,i1,j,k,i2,l,l1,l2;
    double c1,c2,tx,ty,t1,t2,u1,u2,z;

    /* Calculate the number of points */

    n = 1;

    for (i = 0 ; i < m ; i++)
        n *= 2;

    /* Do the bit reversal */

    j = 0;
    i2 = n >> 1;

    for (i = 0 ; i < n - 1 ; i++) {

        if (i < j) {
            tx = x[i];
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }

        k = i2;

        while (k <= j) {
            j -= k;
            k >>= 1;
        }

        j += k;
    }

    /* Compute the FFT */

    c1 = -1.0;
    c2 = 0.0;
    l2 = 1;

    for (l = 0 ; l < m ; l++) {

        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;

        for (j = 0 ; j < l1 ; j++) {

            for (i = j ; i < n ; i += l2) {
                i1 = i + l1;
                t1 = u1 * x[i1] - u2 * y[i1];
                t2 = u1 * y[i1] + u2 * x[i1];
                x[i1] = x[i] - t1;
                y[i1] = y[i] - t2;
                x[i] += t1;
                y[i] += t2;
            }

            z =  u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }

        c2 = sqrt((1.0 - c1) / 2.0);

        if (dir == 1)
            c2 = -c2;

        c1 = sqrt((1.0 + c1) / 2.0);
    }

    /* Scaling for forward transform */

    if (dir == 1) {
        for (i = 0 ; i < n ; i++) {
            x[i] /= n;
            y[i] /= n;
        }
    }

    return 1;
}