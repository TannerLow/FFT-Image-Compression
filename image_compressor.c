#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#define PI 3.14159265359

#define ROWS 2048
#define COLS 2048
#define INPUT_IMAGE "sweet_p_grayscale.bmp"
#define OUTPUT_IMAGE "sweet_p_compressed.bmp"
#define NORMALIZE_PIXELS 1
#define THRESHOLD_PERCENT 0.98f

void fft(float* rex, float* imx, int sampleCount);
void ifft(float* rex, float* imx, int sampleCount);
void printFFTData(FILE* stream, float* rex, float* imx, int sampleCount);

int isLittleEndian();
int compareFloats(const void *a, const void *b);

int main() {
    FILE* outputLogFile = fopen("logs.txt", "w");
    if(outputLogFile == NULL) {
        perror("Error opening file");
        return 1;
    }

    printf("isLittleEndian: %d\n", isLittleEndian());

    FILE *file;
    unsigned char header[54]; // BMP header is 54 bytes
    unsigned int dataPos;     // Position where actual pixel data begins

    // Open the BMP file in binary mode for reading
    file = fopen(INPUT_IMAGE, "rb");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    // Read the BMP header (first 54 bytes)
    fread(header, 1, 54, file);

    // BMP files start with "BM" (ASCII code 66 and 77)
    if (header[0] != 'B' || header[1] != 'M') {
        printf("Not a BMP file\n");
        fclose(file);
        fclose(outputLogFile);
        return 1;
    }

    // Extract the image dimensions and data offset from the header
    dataPos = *(int*)&(header[0x0A]);
    printf("Data Offset: %u\n", dataPos);

    int rows = ROWS;
    int cols = COLS;
    int totalPixelCount = rows * cols;
    
    unsigned char* pixelData;
    pixelData = (unsigned char*)malloc(totalPixelCount);
    if(pixelData == NULL) {
        printf("Failed to allocate memory\n");
        fclose(file);
        fclose(outputLogFile);
        return 1;
    }

    // skip to the pixel data portion of file
    if(fseek(file, dataPos, SEEK_SET) != 0) {
        printf("failed to seek to data position in input file\n");
        fclose(file);
        fclose(outputLogFile);
        return 1;
    }
    
    int pixelsRead = fread(pixelData, 1, totalPixelCount, file);
    if(pixelsRead != totalPixelCount) {
        printf("Failed to read from file. Bytes read did not match expected\n");
        fclose(file);
        fclose(outputLogFile);
        return 1;
    }
    printf("Data check: %d %d %d\n", (int)pixelData[0], (int)pixelData[1], (int)pixelData[totalPixelCount-1]);


    // Perform compression algorithm
    float* reAllCoefficients = (float*)malloc(totalPixelCount * sizeof(float));
    float* imAllCoefficients = (float*)malloc(totalPixelCount * sizeof(float));
    memset(reAllCoefficients, 0, totalPixelCount * sizeof(float));
    memset(imAllCoefficients, 0, totalPixelCount * sizeof(float));

    float* reRowCoefficients = (float*)malloc(cols * sizeof(float));
    float* imRowCoefficients = (float*)malloc(cols * sizeof(float));

    if(reAllCoefficients == NULL || imAllCoefficients == NULL || reRowCoefficients == NULL || imRowCoefficients == NULL) {
        printf("Failed to allocate memory for coefficients\n");
        fclose(outputLogFile);
        fclose(file);
        free(reAllCoefficients);
        free(imAllCoefficients);
        free(reRowCoefficients);
        free(imRowCoefficients);
        return 1;    
    }

    // load in normalized samples and perform FFT
    printf("Beginning FFT of rows\n");

    for(int row = 0; row < rows; row++) {
        memset(reRowCoefficients, 0, cols * sizeof(float));
        memset(imRowCoefficients, 0, cols * sizeof(float));
        for(int col = 0; col < cols; col++) {
            if(NORMALIZE_PIXELS) {
                reRowCoefficients[col] = pixelData[row * cols + col] / 256.f;
            }
            else {
                reRowCoefficients[col] = pixelData[row * cols + col];
            }
        }

        fft(reRowCoefficients, imRowCoefficients, cols);

        memcpy(reAllCoefficients + row * cols, reRowCoefficients, cols * sizeof(float));
        memcpy(imAllCoefficients + row * cols, imRowCoefficients, cols * sizeof(float));

        // debugging purposes
        if(reAllCoefficients[50 + row * cols] != reRowCoefficients[50]) {
            printf("%.02f ", reAllCoefficients[50 + row * cols]);
            printf("%.02f \n", reRowCoefficients[50]);
            printf("row: %d\n", row);
            printf("panic!!!\n");

            fclose(outputLogFile);
            fclose(file);
            free(reAllCoefficients);
            free(imAllCoefficients);
            free(reRowCoefficients);
            free(imRowCoefficients);
            return 1;  
        }
    }

    free(reRowCoefficients);
    free(imRowCoefficients);
    reRowCoefficients = NULL;
    imRowCoefficients = NULL;

    printf("Beginning FFT of columns\n");

    // FFT the columns
    float* reColCoefficients = (float*)malloc(rows * sizeof(float));
    float* imColCoefficients = (float*)malloc(rows * sizeof(float));
    if(reColCoefficients == NULL || imColCoefficients == NULL) {
        printf("Failed to allocate memory for column coefficients\n");
        fclose(outputLogFile);
        fclose(file);
        free(reAllCoefficients);
        free(imAllCoefficients);
        free(reColCoefficients);
        free(imColCoefficients);
        return 1;
    }

    for(int col = 0; col < cols; col++) {
        memset(reColCoefficients, 0, rows * sizeof(float));
        memset(imColCoefficients, 0, rows * sizeof(float));
        for(int row = 0; row < rows; row++) {
            reColCoefficients[row] = reAllCoefficients[row * cols + col];
            imColCoefficients[row] = imAllCoefficients[row * cols + col];
        }

        fft(reColCoefficients, imColCoefficients, rows);

        // logging values for 1 example to showcase data transformation
        if(col == 0) {
            fprintf(outputLogFile, "FFT coefficients after initial column-wise FFT:\n");
            printFFTData(outputLogFile, reColCoefficients, imColCoefficients, rows);
        }

        // save coefficients to set of all coefficients
        for(int row = 0; row < rows; row++) {
            reAllCoefficients[row * cols + col] = reColCoefficients[row];
            imAllCoefficients[row * cols + col] = imColCoefficients[row];
        }

        // debugging purposes
        if(reAllCoefficients[50 * cols + col] != reColCoefficients[50]) {
            printf("%.02f ", reAllCoefficients[50 * cols + col]);
            printf("%.02f \n", reColCoefficients[50]);
            printf("row: %d\n", col);
            printf("panic!!!\n");

            fclose(outputLogFile);
            fclose(file);
            free(reAllCoefficients);
            free(imAllCoefficients);
            free(reColCoefficients);
            free(imColCoefficients);
            return 1;  
        }
    }

    // zero the lowest N% of frequencies by amplitude
    float* sortedByMagnitudes;
    sortedByMagnitudes = (float*)malloc(totalPixelCount * sizeof(float));
    if(sortedByMagnitudes == NULL) {
        printf("Failed to allocate memory when zeroing data\n");
        fclose(outputLogFile);
        fclose(file);
        free(reAllCoefficients);
        free(imAllCoefficients);
        free(reColCoefficients);
        free(imColCoefficients);
        return 1;
    }

    for(int row = 0; row < rows; row++) {
        for(int col = 0; col < cols; col++) {
            int index = row * cols + col;
            sortedByMagnitudes[index]  = reAllCoefficients[index] * reAllCoefficients[index];
            sortedByMagnitudes[index] += imAllCoefficients[index] * imAllCoefficients[index];
            //sortedByMagnitudes[row]  = sqrt(sortedByMagnitudes[row]);
        }
    }
    qsort(sortedByMagnitudes, totalPixelCount, sizeof(float), compareFloats);
    int thresholdIndex = totalPixelCount * THRESHOLD_PERCENT;
    float threshold = sortedByMagnitudes[thresholdIndex];
    
    // zero low thresholds
    for(int row = 0; row < rows; row++) {
        for(int col = 0; col < cols; col++) {
            int index = row * cols + col;
            float magnitude;
            magnitude  = reAllCoefficients[index] * reAllCoefficients[index];
            magnitude += imAllCoefficients[index] * imAllCoefficients[index];
            //magnitude  = sqrt(magnitude);
            if(magnitude < threshold) {
                reAllCoefficients[index] = 0;
                imAllCoefficients[index] = 0;
            }
        }

        // logging values for 1 example to showcase data transformation
        if(row == 256) {
            fprintf(outputLogFile, "FFT coefficients after zeroing %.02f of low frequencies:\nREX: ", THRESHOLD_PERCENT);
            for(int col = 0; col < cols; col++) {
                fprintf(outputLogFile, "%.02f ", reAllCoefficients[row * cols + col]);
            }
            fprintf(outputLogFile, "\n\n");
        }
    }



    // Image is now ready to be compressed

    // Pretend we did compress and save the image
    // To actually compress the file, you could use a simple algorithm like run-length compression
    // Because of all the consecutive zeros, it would compress each contigous set of zeros to just ~8 bytes each
    // Calculation to see how much data would have been saved
    printf("Pixel data byte count with no compression: %d\n", totalPixelCount);
    int hypotheticalByteCount = 0;
    bool reConnectedToZero = false;
    bool imConnectedToZero = false;
    for(int i = 0; i < totalPixelCount; i++) {
        float reCoefficient = reAllCoefficients[i];
        if(reCoefficient == 0 && !reConnectedToZero) {
            hypotheticalByteCount += sizeof(float) + sizeof(int);
            reConnectedToZero = true;
        }
        else if (reCoefficient != 0){
            hypotheticalByteCount += sizeof(float);
            reConnectedToZero = false;
        }

        float imCoefficient = imAllCoefficients[i];
        if(imCoefficient == 0 && !imConnectedToZero) {
            hypotheticalByteCount += sizeof(float) + sizeof(int);
            imConnectedToZero = true;
        }
        else if (imCoefficient != 0){
            hypotheticalByteCount += sizeof(float);
            imConnectedToZero = false;
        }
    }
    printf("Pixel data byte count with some basic compression: %d\n", hypotheticalByteCount);

    // Now to recreate image from compressed data
    // If we had done run-length compression, this is the stage we would decompress at
    printf("Beginning IFFT of columns\n");

    for(int col = 0; col < cols; col++) {
        memset(reColCoefficients, 0, rows * sizeof(float));
        memset(imColCoefficients, 0, rows * sizeof(float));
        for(int row = 0; row < rows; row++) {
            reColCoefficients[row] = reAllCoefficients[row * cols + col];
            imColCoefficients[row] = imAllCoefficients[row * cols + col];
        }

        ifft(reColCoefficients, imColCoefficients, rows);

        // save coefficients to set of all coefficients
        for(int row = 0; row < rows; row++) {
            reAllCoefficients[row * cols + col] = reColCoefficients[row];
            imAllCoefficients[row * cols + col] = imColCoefficients[row];
        }

        // debugging purposes
        if(reAllCoefficients[50 * cols + col] != reColCoefficients[50]) {
            printf("%.02f ", reAllCoefficients[50 * cols + col]);
            printf("%.02f \n", reColCoefficients[50]);
            printf("row: %d\n", col);
            printf("panic!!!\n");

            fclose(outputLogFile);
            fclose(file);
            free(reAllCoefficients);
            free(imAllCoefficients);
            free(reColCoefficients);
            free(imColCoefficients);
            return 1;  
        }
    }

    free(reColCoefficients);
    free(imColCoefficients);
    reColCoefficients = NULL;
    imColCoefficients = NULL;


    printf("Beginning IFFT of rows\n");

    reRowCoefficients = (float*)malloc(cols * sizeof(float));
    imRowCoefficients = (float*)malloc(cols * sizeof(float));
    if(reAllCoefficients == NULL || imAllCoefficients == NULL || reRowCoefficients == NULL || imRowCoefficients == NULL) {
        printf("Failed to allocate memory for coefficients\n");
        fclose(outputLogFile);
        fclose(file);
        free(reAllCoefficients);
        free(imAllCoefficients);
        free(reRowCoefficients);
        free(imRowCoefficients);
        return 1;
    }

    for(int row = 0; row < rows; row++) {
        memset(reRowCoefficients, 0, cols * sizeof(float));
        memset(imRowCoefficients, 0, cols * sizeof(float));
        for(int col = 0; col < cols; col++) {
            reRowCoefficients[col] = reAllCoefficients[row * cols + col];
            imRowCoefficients[col] = imAllCoefficients[row * cols + col];
        }

        ifft(reRowCoefficients, imRowCoefficients, cols);

        memcpy(reAllCoefficients + row * cols, reRowCoefficients, cols * sizeof(float));
        memcpy(imAllCoefficients + row * cols, imRowCoefficients, cols * sizeof(float));

        // debugging purposes
        if(reAllCoefficients[50 + row * cols] != reRowCoefficients[50]) {
            printf("%.02f ", reAllCoefficients[50 + row * cols]);
            printf("%.02f \n", reRowCoefficients[50]);
            printf("row: %d\n", row);
            printf("panic!!!\n");

            fclose(outputLogFile);
            fclose(file);
            free(reAllCoefficients);
            free(imAllCoefficients);
            free(reRowCoefficients);
            free(imRowCoefficients);
            return 1;  
        }
    }
    free(reRowCoefficients);
    free(imRowCoefficients);
    free(imAllCoefficients);
    reRowCoefficients = NULL;
    imRowCoefficients = NULL;
    imAllCoefficients = NULL;


    // save image to file
    if(fseek(file, 0, SEEK_SET) != 0) {
        printf("failed to seek to beginning input file\n");
        fclose(file);
        free(reAllCoefficients);
        free(pixelData);
        fclose(outputLogFile);
        return 1;
    }

    char* prePixelData = (char*)malloc(dataPos);
    int bytesRead = fread(prePixelData, sizeof(char), dataPos, file);
    fclose(file);
    if(bytesRead != dataPos) {
        printf("failed to read pre data portion of input file\n");
        free(reAllCoefficients);
        free(pixelData);
        fclose(outputLogFile);
        return 1;
    }
    
    for(int i = 0; i < totalPixelCount; i++) {
        unsigned char pixelValue;
        if(NORMALIZE_PIXELS) {
            pixelValue = reAllCoefficients[i] * 256;
        }
        else {
            pixelValue = reAllCoefficients[i];
        }
        pixelData[i] = pixelValue;
    }

    FILE* outputFile = fopen(OUTPUT_IMAGE, "wb");
    if(outputFile != NULL) {
        if((fwrite(prePixelData, sizeof(char), dataPos, outputFile) == dataPos) &&
           (fwrite(pixelData, sizeof(char), totalPixelCount, outputFile) == totalPixelCount)) {
            printf("Image compressed successfully!\n");
        }
    }
    else {
        printf("Failed to open output file\n");
    }

    
    // clean up the mess
    free(reAllCoefficients);
    free(pixelData);
    fclose(outputLogFile);
    fclose(outputFile);
    return 0;
}

void fft(float* rex, float* imx, int sampleCount) {
    // constants
    int last = sampleCount - 1;
    int halfCount = sampleCount / 2;
    int log = ilogb(sampleCount);
    int j = halfCount;

    // bit reversal ordering
    for(int i = 1; i < last; i++) {
        if(i < j) {
            float tempReal = rex[j];
            float tempImag = imx[j];
            rex[j] = rex[i];
            imx[j] = imx[i];
            rex[i] = tempReal;
            imx[i] = tempImag;
        }

        int k = halfCount;
        while(k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }

    // Iterative FFT algorithm
    for(int stage = 1; stage <= log; stage++) {
        int stageSampleCount = exp2(stage);
        halfCount = stageSampleCount / 2;
        float ur = 1;
        float ui = 0;
        float realSinusoid =  cos(PI / halfCount);
        float imagSinusoid = -sin(PI / halfCount);
        for(int j = 0; j < halfCount; j++) {
            for(int i = j; i < sampleCount; i += stageSampleCount) {
                int ptr = i + halfCount;
                float tempReal = rex[ptr] * ur - imx[ptr] * ui;
                float tempImag = rex[ptr] * ui + imx[ptr] * ur;
                rex[ptr] = rex[i] - tempReal;
                imx[ptr] = imx[i] - tempImag;
                rex[i] += tempReal;
                imx[i] += tempImag;
            }
            float tempUr = ur;
            ur = ur * realSinusoid - ui * imagSinusoid;
            ui = tempUr * imagSinusoid + ui * realSinusoid;
        }
    }
}

void ifft(float* rex, float* imx, int sampleCount) {
    for(int i = 0; i < sampleCount; i++) {
        imx[i] = -imx[i];
    }

    fft(rex, imx, sampleCount);

    for(int i = 0; i < sampleCount; i++) {
        rex[i] /= sampleCount;
        imx[i] /= -sampleCount;
    }
}

void printFFTData(FILE* stream, float* rex, float* imx, int sampleCount) {
    fprintf(stream, "REX: ");
    for(int i = 0; i < sampleCount; i++) {
        fprintf(stream, "%.02f ", rex[i]);
    }
    fprintf(stream, "\n");

    fprintf(stream, "IMX: ");
    for(int i = 0; i < sampleCount; i++) {
        fprintf(stream, "%.02f ", imx[i]);
    }
    fprintf(stream, "\n\n");
}

int isLittleEndian() {
    unsigned int num = 1;
    // Cast the integer to a byte pointer
    unsigned char *ptr = (unsigned char *)&num;
    // If the first byte is the least significant byte (little endian), return 1
    // Otherwise, return 0 (big endian)
    return (*ptr == 1);
}

int compareFloats(const void *a, const void *b) {
    return (*(float*)a - *(float*)b);
}