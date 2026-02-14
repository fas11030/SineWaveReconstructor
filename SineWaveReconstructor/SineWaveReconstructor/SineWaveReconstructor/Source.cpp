/*

#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Simple FFT Implementation ---
void fft(std::vector<std::complex<double>>& x) {
    size_t n = x.size();
    if (n <= 1) return;
    std::vector<std::complex<double>> even(n / 2), odd(n / 2);
    for (size_t i = 0; i < n / 2; ++i) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }
    fft(even);
    fft(odd);
    for (size_t k = 0; k < n / 2; ++k) {
        std::complex<double> t = std::polar(1.0, -2.0 * M_PI * k / (double)n) * odd[k];
        x[k] = even[k] + t;
        x[k + n / 2] = even[k] - t;
    }
}

// --- Helper Functions ---
size_t pow2_floor(size_t n) {
    size_t p = 1;
    while (p * 2 <= n) p *= 2;
    return p;
}

double interpolate_peak(double m1, double m2, double m3) {
    double den = (m1 + m3) - (2.0 * m2);
    return (std::abs(den) < 1e-9) ? 0.0 : (m1 - m3) / (2.0 * den);
}

// --- The Core Logic ---
void process_audio(const std::vector<double>& input, std::vector<double>& output,
    double sampleRate, int numSines, int density) {

    size_t selSamps = input.size();
    size_t winlen = std::min((size_t)16384, pow2_floor(std::max((size_t)1024, selSamps)));
    size_t fftlen = winlen;
    size_t halfLen = fftlen / 2;
    double binWidth = sampleRate / (double)fftlen;

    std::vector<double> totalMags(halfLen, 0.0);
    double skipDistSamps = (double)selSamps / (double)density;

    // ANALYSIS
    for (int s = 0; s < density; ++s) {
        size_t startIdx = (size_t)(s * skipDistSamps);
        if (startIdx + winlen > selSamps) break;

        std::vector<std::complex<double>> fftData(fftlen, 0.0);
        for (size_t i = 0; i < winlen; ++i) {
            double window = 0.5 * (1.0 - cos((2.0 * M_PI * i) / (double)(winlen - 1)));
            fftData[i] = input[startIdx + i] * window;
        }

        fft(fftData);

        for (size_t i = 0; i < halfLen; ++i) {
            totalMags[i] += std::abs(fftData[i]);
        }
    }

    // RESYNTHESIS
    output.assign(selSamps, 0.0);
    for (int k = 0; k < numSines; ++k) {
        double maxVal = -1.0;
        size_t maxIdx = 0;

        for (size_t j = 2; j < halfLen - 2; ++j) {
            if (totalMags[j] > maxVal && totalMags[j] > totalMags[j - 1] && totalMags[j] > totalMags[j + 1]) {
                maxVal = totalMags[j];
                maxIdx = j;
            }
        }

        if (maxVal > 0) {
            double delta = interpolate_peak(totalMags[maxIdx - 1], totalMags[maxIdx], totalMags[maxIdx + 1]);
            double freq = ((double)maxIdx + delta) * binWidth;
            double phaseInc = (2.0 * M_PI * freq) / sampleRate;
            double phase = 0.0;

            for (size_t i = 0; i < selSamps; ++i) {
                output[i] += maxVal * sin(phase);
                phase += phaseInc;
            }
            // Wipe spectral peak to find the next one
            totalMags[maxIdx] = totalMags[maxIdx - 1] = totalMags[maxIdx + 1] = 0.0;
        }
    }

    // NORMALIZE
    double maxOut = 0.0;
    for (double s : output) maxOut = std::max(maxOut, std::abs(s));
    if (maxOut > 1e-9) {
        double scale = 0.8 / maxOut;
        for (double& s : output) s *= scale;
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cout << "Usage: sine_recon input.wav output.wav" << std::endl;
        return 1;
    }

    unsigned int channels;
    unsigned int sampleRate;
    drwav_uint64 totalFrameCount;

    // Using f32 version for better compatibility
    float* pSampleData = drwav_open_file_and_read_pcm_frames_f32(argv[1], &channels, &sampleRate, &totalFrameCount, NULL);

    if (pSampleData == NULL) {
        std::cerr << "Error: Could not open or read input file." << std::endl;
        return 1;
    }

    std::vector<double> input((size_t)totalFrameCount);
    for (size_t i = 0; i < (size_t)totalFrameCount; ++i) {
        input[i] = (double)pSampleData[i * channels]; // Take first channel
    }

    std::vector<double> output;
    process_audio(input, output, (double)sampleRate, 10, 20);

    drwav_data_format format;
    format.container = drwav_container_riff;
    format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    format.channels = 1;
    format.sampleRate = sampleRate;
    format.bitsPerSample = 32;

    drwav wav;
    if (!drwav_init_file_write(&wav, argv[2], &format, NULL)) {
        std::cerr << "Error: Could not open output file for writing." << std::endl;
        return 1;
    }

    // Convert output to float for writing
    std::vector<float> outputF32(output.begin(), output.end());
    drwav_write_pcm_frames(&wav, totalFrameCount, outputF32.data());
    drwav_uninit(&wav);

    drwav_free(pSampleData, NULL);
    std::cout << "Success! Created: " << argv[2] << std::endl;
    return 0;
}
*/

#define DR_WAV_IMPLEMENTATION
#include "dr_wav.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Core Math ---
void fft(std::vector<std::complex<double>>& x) {
    size_t n = x.size();
    if (n <= 1) return;
    std::vector<std::complex<double>> even(n / 2), odd(n / 2);
    for (size_t i = 0; i < n / 2; ++i) {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }
    fft(even); fft(odd);
    for (size_t k = 0; k < n / 2; ++k) {
        std::complex<double> t = std::polar(1.0, -2.0 * M_PI * k / (double)n) * odd[k];
        x[k] = even[k] + t;
        x[k + n / 2] = even[k] - t;
    }
}

double get_peak(const std::vector<double>& v) {
    double p = 0.0;
    for (double s : v) p = std::max(p, std::abs(s));
    return p;
}

// Processing function for a single channel
void process_channel(const std::vector<double>& input, std::vector<double>& output,
    double sr, int count, int density, bool preserve) {

    size_t selSamps = input.size();
    size_t winlen = std::min((size_t)16384, (size_t)pow(2, floor(log2(std::max((double)1024, (double)selSamps)))));
    size_t fftlen = winlen;
    size_t halfLen = fftlen / 2;
    double binWidth = sr / (double)fftlen;
    std::vector<double> totalMags(halfLen, 0.0);

    // Analysis
    double skip = (double)selSamps / density;
    for (int s = 0; s < density; ++s) {
        size_t start = (size_t)(s * skip);
        if (start + winlen > selSamps) break;
        std::vector<std::complex<double>> fftData(fftlen, 0.0);
        for (size_t i = 0; i < winlen; ++i) {
            double w = 0.5 * (1.0 - cos((2.0 * M_PI * i) / (winlen - 1)));
            fftData[i] = input[start + i] * w;
        }
        fft(fftData);
        for (size_t i = 0; i < halfLen; ++i) totalMags[i] += std::abs(fftData[i]);
    }

    // Resynthesis
    output.assign(selSamps, 0.0);
    for (int k = 0; k < count; ++k) {
        double maxV = -1.0; size_t maxI = 0;
        for (size_t j = 2; j < halfLen - 2; ++j) {
            if (totalMags[j] > maxV && totalMags[j] > totalMags[j - 1] && totalMags[j] > totalMags[j + 1]) {
                maxV = totalMags[j]; maxI = j;
            }
        }
        if (maxV > 0) {
            double den = (totalMags[maxI - 1] + totalMags[maxI + 1]) - (2.0 * totalMags[maxI]);
            double delta = (den == 0) ? 0.0 : (totalMags[maxI - 1] - totalMags[maxI + 1]) / (2.0 * den);
            double freq = (maxI + delta) * binWidth;
            double phaseInc = (2.0 * M_PI * freq) / sr;
            double phase = 0.0;
            for (size_t i = 0; i < selSamps; ++i) {
                output[i] += totalMags[maxI] * sin(phase);
                phase += phaseInc;
            }
            totalMags[maxI] = totalMags[maxI - 1] = totalMags[maxI + 1] = 0.0;
        }
    }

    // Scaling
    double inPeak = get_peak(input);
    double outPeak = get_peak(output);
    if (outPeak > 0) {
        double factor = preserve ? (inPeak / outPeak) : (0.8 / outPeak);
        for (double& s : output) s *= factor;
    }
}

int main() {
    std::string inPath, outPath;
    int count, density, choice;
    bool preserve;

    // --- Simple "Prompt" GUI ---
    std::cout << "--- Sine Peak Reconstructor ---\n";
    std::cout << "Input WAV path: "; std::getline(std::cin, inPath);
    if (inPath.front() == '\"') inPath = inPath.substr(1, inPath.size() - 2); // Clean quotes

    std::cout << "Output WAV path: "; std::getline(std::cin, outPath);
    std::cout << "Number of Sines (1-100): "; std::cin >> count;
    std::cout << "Density (Snapshots): "; std::cin >> density;
    std::cout << "Preserve Original Amplitude? (1=Yes, 0=No): "; std::cin >> choice;
    preserve = (choice == 1);

    // Read Audio
    unsigned int channels;
    unsigned int sampleRate;
    drwav_uint64 totalFrames;
    float* pData = drwav_open_file_and_read_pcm_frames_f32(inPath.c_str(), &channels, &sampleRate, &totalFrames, NULL);
    if (!pData) return 1;

    std::vector<std::vector<double>> inputChans(channels, std::vector<double>(totalFrames));
    std::vector<std::vector<double>> outputChans(channels, std::vector<double>(totalFrames));

    // De-interleave
    for (size_t i = 0; i < totalFrames; ++i) {
        for (size_t c = 0; c < channels; ++c) {
            inputChans[c][i] = pData[i * channels + c];
        }
    }

    // Process each channel
    for (size_t c = 0; c < channels; ++c) {
        std::cout << "Processing Channel " << c + 1 << "...\n";
        process_channel(inputChans[c], outputChans[c], sampleRate, count, density, preserve);
    }

    // Interleave and Write
    std::vector<float> finalOut(totalFrames * channels);
    for (size_t i = 0; i < totalFrames; ++i) {
        for (size_t c = 0; c < channels; ++c) {
            finalOut[i * channels + c] = (float)outputChans[c][i];
        }
    }

    drwav_data_format format = { drwav_container_riff, DR_WAVE_FORMAT_IEEE_FLOAT, channels, sampleRate, 32 };
    drwav wav;
    drwav_init_file_write(&wav, outPath.c_str(), &format, NULL);
    drwav_write_pcm_frames(&wav, totalFrames, finalOut.data());
    drwav_uninit(&wav);

    drwav_free(pData, NULL);
    std::cout << "Done! Press Enter to exit.";
    std::cin.ignore(); std::cin.get();
    return 0;
}