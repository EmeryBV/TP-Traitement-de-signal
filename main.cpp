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
#define duree 0.2

void DFT(double *signal, double *partie_reelle, double *partie_imaginaire, int N);
void IDFT(double *signal, double *partie_reelle, double *partie_imaginaire, int N);
unsigned char normaliser_double_to_unsiCharn(double signal) {
    return (signal + 1.0) * 127.5;
}

double normaliser_unsiChar_to_double(unsigned char signal) {
    return (signal / 127.5) - 1;
}

void Normalisation(double *signal, int N);

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

vector<double> DFT_visualize(vector<double> &signal_real, vector<double> &signal_imaginary) {

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


int main() {
    // en secondes
    char *filename = "sounds/La.wav";
    char *filenameDFT = "sounds/LaDFT.wav";
    char *filenameBis = "sounds/Labis.wav";
    short nb_channels = 1;

    vector<double> partition;
    partition.push_back(LA);
      partition.push_back(RE);
      partition.push_back(MI);
//      partition.push_back(FA);
//      partition.push_back(SOL);
//      partition.push_back(LA);
//      partition.push_back(SI);
//      partition.push_back(DO);


    unsigned char data8[static_cast<int>(sampling_freq * duree)];
    creerMusique(partition, data8);
    Wave wave = Wave(data8, sampling_freq * duree, nb_channels, sampling_freq);
    wave.write(filename);
    vector<double> signal = read_signal(filename, 1);
    vector<unsigned char> data8bis(to_data8(signal));
    Wave waveBis = Wave(data8bis.data(), sampling_freq * duree, nb_channels, sampling_freq);
    waveBis.write(filenameBis);

    vector<double> partie_reelle(signal.size());
    vector<double> partie_imaginaire(signal.size());
    cout << "Transformation de fourrier en cours..." << endl;

    DFT(signal.data(), partie_reelle.data(), partie_imaginaire.data(), signal.size());

    double norme = sqrt( pow(partie_reelle.size(),2) +pow(partie_imaginaire.size(),2));
    cout<<"la norme est "<< norme;
    copy(partie_reelle.begin(), partie_reelle.begin() + signal.size(), signal.begin());
    vector<double> test = DFT_visualize(partie_reelle,partie_imaginaire);
//    Normalisation(partie_reelle.data(), partie_reelle.size());
    write_signal(filenameDFT,test , 1);

//    vector<unsigned char> data8DFT(to_data8(partie_reelle));
//    Wave waveDFT = Wave(data8DFT.data(), data8DFT.size(), nb_channels, sampling_freq);
//    waveDFT.write(filenameDFT);

    return 0;

}

void DFT(double *signal, double *partie_reelle, double *partie_imaginaire, int N) {
    int n, k;
    double deux_pi_N = 2.0 * M_PI / (double) N;
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

void IDFT(double *signal, double *partie_reelle, double *partie_imaginaire, int N) {
    int n, k;
    double deux_pi_N = 2.0 * M_PI / (double) N;
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


void Normalisation(double *signal, int N) {
    double mini, maxi, delta;
    mini = signal[0];
    maxi = mini;

    for (int i = 0; i < N; i++) {
        cout<<"Max = "<<signal[i]<<"\n";
        maxi = maxi >  (signal[i]) ? (signal[i]) : maxi;
        mini = mini >  (signal[i]) ? (signal[i]) : mini;

    }
    cout<<"Delta = "<<maxi - mini<<"\n";
    cout<<"Max = "<<maxi<<"\n";
    cout<<"Min= "<<mini<<"\n";
    delta = maxi - mini;
    if(delta>1e-16) delta =1.0/delta;
    for(int i = 0;i<N;i++){
//        cout<<signal[i]<<"\n";
//        cout<<delta<<"\n";
        signal[i] = (signal[i]-mini)*delta;
//        cout<<signal[i]<<"\n";
    }

}