#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#define Q31_1_BASE (1LL<<30)
#define Q31_1_MAX  ((1<<30)-1)
#define Q31_1_MIN  (-1<<30)
#define SIZE (128)
#define PI (3.14159265358979323846)

int32_t flt2fixd(double x)
{
	if (x >= 1)
		return Q31_1_MAX;
	else if (x < -1)
		return Q31_1_MIN;

	int32_t res = x * (double)Q31_1_BASE;
	return res;
}

float fixd2flt(int32_t x)
{
	float res = (float)(x) / ((float)Q31_1_BASE);
	return res;
	return res;
}

struct HEADER
{
	uint8_t ChunkID[4];
	uint32_t ChunkSize;
	uint8_t Format[4];
	uint8_t Subchunk1ID[4];
	uint32_t Subchunk1Size;
	uint16_t AudioFormat;
	uint16_t NumChannels;
	uint32_t SampleRate;
	uint32_t ByteRate;
	uint16_t BlockAlign;
	uint16_t BitsPreSample;
};

struct CHUNK
{
	uint8_t ID[4];
	uint32_t size;
};

void coeffscalc(int32_t* coeffs, double* coeffs_double, float filterfreq, int32_t samplerate, float Q)
{
	double a0, a1, a2, b1, b2, norm;

	double K = tan(PI * filterfreq / samplerate);
	double K2 = K * K;

	norm = 1.0 / (1.0 + K / Q + K2);
	a0 = K2 * norm;
	a1 = 2 * a0;
	a2 = a0;
	b1 = 2 * (K2 - 1) * norm;
	b2 = (1 - K / Q + K2) * norm;

	coeffs_double[0] = a0;
	coeffs_double[1] = a1;
	coeffs_double[2] = a2;
	coeffs_double[3] = b1;
	coeffs_double[4] = b2;

	coeffs[0] = flt2fixd(a0);
	coeffs[1] = flt2fixd(a1);
	coeffs[2] = flt2fixd(a2);
	coeffs[3] = flt2fixd(b1);
	coeffs[4] = flt2fixd(b2);
}

int64_t acc = 0;

int32_t IIR(int32_t* buffer, int32_t* coeffs,int16_t sample)
{
	buffer[0] = buffer[1];
	buffer[1] = buffer[2];
	buffer[2] = (int32_t)sample;
	buffer[3] = buffer[4];
	buffer[4] = buffer[5];

	acc += (int64_t)buffer[2] * coeffs[0] + (int64_t)buffer[1] * coeffs[1] + (int64_t)buffer[0] * coeffs[2] - (int64_t)buffer[4] * coeffs[3] - (int64_t)buffer[3] * coeffs[4];

	buffer[5] = (acc >> 30);
	acc = acc & 0x3fffffff;
	return buffer[5];
}

int main() {

	FILE* file_in;
	FILE* file_out;
	struct HEADER header;
	struct CHUNK chunk;

	int8_t filename[30];
	puts("Enter filename:");
	gets(filename);

	fopen_s(&file_in, filename, "rb");
	fopen_s(&file_out, "filteredwave.wav", "wb");

	fread(&header, sizeof(header), 1, file_in);

	printf("NumChannels=%i\n", header.NumChannels);
	printf("SampleRate=%i\n", header.SampleRate);
	printf("BitsPreSample=%i\n", header.BitsPreSample);

	if (header.BitsPreSample != 16) {
		printf("Audiofile isn't 16bit\n");
		system("pause");
		return 0;
	}

	fread(&chunk, sizeof(chunk), 1, file_in);
	
	uint16_t sample_size = header.BitsPreSample / 8;
	uint32_t samples_count = chunk.size / sample_size;
	printf("Samples=%d\n", samples_count);

	fwrite(&header, sizeof(header), 1, file_out);
	fwrite(&chunk, sizeof(chunk), 1, file_out);

	int32_t sample_buffer_L[6] = { 0, 0, 0, 0, 0, 0 };
	int32_t sample_buffer_R[6] = { 0, 0, 0, 0, 0, 0 };
	int16_t s_buf[2*SIZE];
	int16_t outputL;
	int16_t outputR;
	int16_t signalL;
	int16_t signal;

	double a0, a1, a2, b1, b2;
	int32_t coeffs[5] = { 0, 0, 0, 0, 0 };

	a0 = 0.11205483175086794;
	a1 = 0.2241096635017359;
	a2 = 0.11205483175086794;
	b1 = -0.8559866467934273;
	b2 = 0.3042059737968993;

	coeffs[0] = flt2fixd(a0);
	coeffs[1] = flt2fixd(a1);
	coeffs[2] = flt2fixd(a2);
	coeffs[3] = flt2fixd(b1);
	coeffs[4] = flt2fixd(b2);

	while (1)
	{
		size_t read_size = fread(s_buf, 2*sample_size, 128, file_in);

		if (!read_size)
			break;

		for (int n = 0; n < read_size; n++) {

			outputL =  IIR(sample_buffer_L, coeffs,s_buf[2*n]);
			outputR = IIR(sample_buffer_R, coeffs, s_buf[2 * n + 1]);

			s_buf[2*n] = outputL;
			s_buf[2 * n + 1] = outputR;
		}
		fwrite(s_buf, 2*sample_size, read_size, file_out);
	}

	fclose(file_in);
	fclose(file_out);
	printf("\n WAV file has been filtered.\n");
	system("pause");
	return 0;
}