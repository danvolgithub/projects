//******************************************************************************
//  Demo Application03 for MSP430/ADS1293 Interface Code Library v1.0
//  Stream Read ADC Data with interrupt
//
//                MSP430F5529
//             -----------------
//         /|\|              XIN|-
//          | |                 |
//          --|RST          XOUT|-
//            |                 |
//            |    P4.4/UCA1SIMO|--> SDI
//            |    P4.5/UCA1SOMI|<-- SDO
//            |     P4.0/UCA1CLK|--> CLK
//            |             P2.0|--> CSB
//            |             P2.4|<-- DRDYB
//            |                 |
//
//   Vishy Natarajan
//   Texas Instruments Inc.
//   April 2013
//   Built with IAR Embedded Workbench Version:  5.5x
//******************************************************************************
/*  Copyright 2011-2012 Texas Instruments Incorporated. All rights reserved.

  IMPORTANT: Your use of this Software is limited to those specific rights
  granted under the terms of a software license agreement between the user who
  downloaded the software, his/her employer (which must be your employer) and
  Texas Instruments Incorporated (the "License"). You may not use this Software
  unless you agree to abide by the terms of the License. The License limits your
  use, and you acknowledge, that the Software may not be modified, copied or
  distributed unless embedded on a Texas Instruments microcontroller which is
  integrated into your product. Other than for the foregoing purpose, you may
  not use, reproduce, copy, prepare derivative works of, modify, distribute,
  perform, display or sell this Software and/or its documentation for any
  purpose.

  YOU FURTHER ACKNOWLEDGE AND AGREE THAT THE SOFTWARE AND DOCUMENTATION ARE
  PROVIDED “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESS OR IMPLIED,
  INCLUDING WITHOUT LIMITATION, ANY WARRANTY OF MERCHANTABILITY, TITLE,
  NON-INFRINGEMENT AND FITNESS FOR A PARTICULAR PURPOSE. IN NO EVENT SHALL TEXAS
  INSTRUMENTS OR ITS LICENSORS BE LIABLE OR OBLIGATED UNDER CONTRACT,
  NEGLIGENCE, STRICT LIABILITY, CONTRIBUTION, BREACH OF WARRANTY, OR OTHER LEGAL
  EQUITABLE THEORY ANY DIRECT OR INDIRECT DAMAGES OR EXPENSES INCLUDING BUT NOT
  LIMITED TO ANY INCIDENTAL, SPECIAL, INDIRECT, PUNITIVE OR CONSEQUENTIAL
  DAMAGES, LOST PROFITS OR LOST DATA, COST OF PROCUREMENT OF SUBSTITUTE GOODS,
  TECHNOLOGY, SERVICES, OR ANY CLAIMS BY THIRD PARTIES (INCLUDING BUT NOT
  LIMITED TO ANY DEFENSE THEREOF), OR OTHER SIMILAR COSTS.

  Should you have any questions regarding your right to use this Software,
  contact Texas Instruments Incorporated at www.TI.com.
*******************************************************************************/

#include <stdint.h>
#include "TI_ADS1293.h"
#include "TI_ADS1293_register_settings.h"
#include "TI_MSP430.h"
#include "TI_MSP430_hardware_board.h"
#include "TI_MSP430_spi.h"

#define WINDOWSIZE 20   // Integrator window size, in samples. The article recommends 150ms. So, FS*0.15.
                        // However, you should check empirically if the waveform looks ok.
#define NOSAMPLE -32000 // An indicator that there are no more samples to read. Use an impossible value for a sample.
#define FS 360          // Sampling frequency.
#define BUFFSIZE 600    // The size of the buffers (in samples). Must fit more than 1.66 times an RR interval, which
                        // typically could be around 1 second.

#define DELAY 22        // Delay introduced by the filters. Filter only output samples after this one.
                        // Set to 0 if you want to keep the delay. Fixing the delay results in DELAY less samples
                        // in the final end result.

#include "panTompkins.h"
#include <stdio.h>      // Remove if not using the standard file functions.

void output();
void TI_ADS1293_WriteRegSettings(void);                                        // Configure ADS1293 registers
void process_adc_output(uint32_t *);                                           // dummy adc data process function: toggles LED
volatile uint8_t ADS1293_ADCDataReady = 0;                                     // set to 1 in DRDYB interrupt service routine
#define SAMPLE_ARRAY_SIZE 10                                                   // Store last 10 ecg samples for processing
#define CH_DATA_SIZE 6                                                         // 6 bytes: ch1 ecg + ch2 ecg
uint32_t adc_sample_array[SAMPLE_ARRAY_SIZE];
//******************************************************************************
void main(void)
{
  uint8_t count, i;
  uint8_t read_buf[CH_DATA_SIZE];
  uint32_t adc_data;


  WDTCTL = WDTPW+WDTHOLD;                                                      // Stop WDT

  TI_ADS1293_LED_PxOUT |= TI_ADS1293_LED_PIN;                                  // Set LED ON
  TI_ADS1293_LED_PxDIR |= TI_ADS1293_LED_PIN;                                  // Set pin direction is output

  // configure Port Pin to handle Data Ready Bar Output (DRDYB) from ADS1293
  TI_ADS1293_DRDYB_PxDIR &= ~TI_ADS1293_DRDYB_PIN;                             // Set up port pin for DRDYB
  TI_ADS1293_DRDYB_PxIES |= TI_ADS1293_DRDYB_PIN;                              // Interrupt Edge Select
  TI_ADS1293_DRDYB_PxIFG &= ~TI_ADS1293_DRDYB_PIN;                             // Clear Interrupt Flag
  TI_ADS1293_DRDYB_PxIE |= TI_ADS1293_DRDYB_PIN;                               // Enable Port interrupt

  TI_ADS1293_SPISetup();                                                       // Initilaize MSP430 SPI Block

  TI_ADS1293_WriteRegSettings();                                               // Set up ADS1293 for Channel Scan

  count = CH_DATA_SIZE;                                                        // bytes to read: ADC_DOUT2 - ADCDOUT0
  i = 0;
  TI_ADS1293_SPIWriteReg(TI_ADS1293_CONFIG_REG,
                          TI_ADS1293_CONFIG_REG_VALUE | ADS1293_START_CONV);   // start conversion
  while (1)
  {
    if (ADS1293_ADCDataReady)
    {
      ADS1293_ADCDataReady = 0;                                                // clear flag
      TI_ADS1293_SPIStreamReadReg(read_buf, count);                            // read adc output into read_buf

      adc_data = ((uint32_t) read_buf[0] << 16)
                 | ((uint16_t) read_buf[1] << 8) | read_buf[2];                // form raw adc output data
      adc_sample_array[i] = adc_data;
      if (++i == SAMPLE_ARRAY_SIZE)                                            // sample array is full
      {
        process_adc_output(adc_sample_array);                                  // dummy app function: no error toggles led
        i = 0;
      }
    }
    __bis_SR_register(LPM0_bits + GIE);                                        // Enter LPM0, enable interrupts
    __no_operation();                                                          // For debugger
  }
}
//******************************************************************************
//  void TI_ADS1293_WriteRegSettings(void)
//
//  DESCRIPTION:
//  ADS1293 registers are configured to the values defined ADS1293_register_settings.h
//  These register settings can easily be obtained from the "Register configuration file"
//  saved from Sensor AFE Software
//
//  ARGUMENTS:
//
//******************************************************************************

void TI_ADS1293_WriteRegSettings(void)
{

  TI_ADS1293_SPIWriteReg(TI_ADS1293_CONFIG_REG,
                          TI_ADS1293_CONFIG_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_FLEX_CH1_CN_REG,
                          TI_ADS1293_FLEX_CH1_CN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_FLEX_CH2_CN_REG,
                          TI_ADS1293_FLEX_CH2_CN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_FLEX_CH3_CN_REG,
                          TI_ADS1293_FLEX_CH3_CN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_FLEX_PACE_CN_REG,
                          TI_ADS1293_FLEX_PACE_CN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_FLEX_VBAT_CN_REG,
                          TI_ADS1293_FLEX_VBAT_CN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_LOD_CN_REG,
                          TI_ADS1293_LOD_CN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_LOD_EN_REG,
                          TI_ADS1293_LOD_EN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_LOD_CURRENT_REG,
                          TI_ADS1293_LOD_CURRENT_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_LOD_AC_CN_REG,
                          TI_ADS1293_LOD_AC_CN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_CMDET_EN_REG,
                          TI_ADS1293_CMDET_EN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_CMDET_CN_REG,
                          TI_ADS1293_CMDET_CN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_RLD_CN_REG,
                          TI_ADS1293_RLD_CN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_WILSON_EN1_REG,
                          TI_ADS1293_WILSON_EN1_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_WILSON_EN2_REG,
                          TI_ADS1293_WILSON_EN2_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_WILSON_EN3_REG,
                          TI_ADS1293_WILSON_EN3_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_WILSON_CN_REG,
                          TI_ADS1293_WILSON_CN_REG_VALUE);

  TI_ADS1293_SPIWriteReg(TI_ADS1293_REF_CN_REG,
                          TI_ADS1293_REF_CN_REG_VALUE);

  TI_ADS1293_SPIWriteReg(TI_ADS1293_OSC_CN_REG,
                          TI_ADS1293_OSC_CN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_AFE_RES_REG,
                          TI_ADS1293_AFE_RES_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_AFE_SHDN_CN_REG,
                          TI_ADS1293_AFE_SHDN_CN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_AFE_FAULT_CN_REG,
                          TI_ADS1293_AFE_FAULT_CN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_AFE_DITHER_EN_REG,
                          TI_ADS1293_AFE_DITHER_EN_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_AFE_PACE_CN_REG,
                          TI_ADS1293_AFE_PACE_CN_REG_VALUE);

  TI_ADS1293_SPIWriteReg(TI_ADS1293_R2_RATE_REG,
                          TI_ADS1293_R2_RATE_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_R3_RATE1_REG,
                          TI_ADS1293_R3_RATE1_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_R3_RATE2_REG,
                          TI_ADS1293_R3_RATE2_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_R3_RATE3_REG,
                          TI_ADS1293_R3_RATE3_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_P_DRATE_REG,
                          TI_ADS1293_P_DRATE_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_DIS_EFILTER_REG,
                          TI_ADS1293_DIS_EFILTER_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_DRDYB_SRC_REG,
                          TI_ADS1293_DRDYB_SRC_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_SYNCOUTB_SRC_REG,
                          TI_ADS1293_SYNCOUTB_SRC_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_MASK_DRDYB_REG,
                          TI_ADS1293_MASK_DRDYB_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_MASK_ERR_REG,
                          TI_ADS1293_MASK_ERR_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_ALARM_FILTER_REG,
                          TI_ADS1293_ALARM_FILTER_REG_VALUE);
  TI_ADS1293_SPIWriteReg(TI_ADS1293_CH_CNFG_REG,
                          TI_ADS1293_CH_CNFG_REG_VALUE);

}
//******************************************************************************
//  void process_adc_output(uint32_t *)
//
//  DESCRIPTION:
//  Dummy ADC data process function: Toggles LED
//******************************************************************************
void process_adc_output(uint32_t *data)
{
    // The signal array is where the most recent samples are kept. The other arrays are the outputs of each
    // filtering module: DC Block, low pass, high pass, integral etc.
    // The output is a buffer where we can change a previous result (using a back search) before outputting.
    uint32_t signal[BUFFSIZE], dcblock[BUFFSIZE], lowpass[BUFFSIZE], highpass[BUFFSIZE], derivative[BUFFSIZE], squared[BUFFSIZE], integral[BUFFSIZE], outputSignal[BUFFSIZE];

    // rr1 holds the last 8 RR intervals. rr2 holds the last 8 RR intervals between rrlow and rrhigh.
    // rravg1 is the rr1 average, rr2 is the rravg2. rrlow = 0.92*rravg2, rrhigh = 1.08*rravg2 and rrmiss = 1.16*rravg2.
    // rrlow is the lowest RR-interval considered normal for the current heart beat, while rrhigh is the highest.
    // rrmiss is the longest that it would be expected until a new QRS is detected. If none is detected for such
    // a long interval, the thresholds must be adjusted.
    int rr1[8], rr2[8], rravg1, rravg2, rrlow = 0, rrhigh = 0, rrmiss = 0;

    // i and j are iterators for loops.
    // sample counts how many samples have been read so far.
    // lastQRS stores which was the last sample read when the last R sample was triggered.
    // lastSlope stores the value of the squared slope when the last R sample was triggered.
    // currentSlope helps calculate the max. square slope for the present sample.
    // These are all long unsigned int so that very long signals can be read without messing the count.
    long unsigned int i, j, sample = 0, lastQRS = 0, lastSlope = 0, currentSlope = 0;

    // This variable is used as an index to work with the signal buffers. If the buffers still aren't
    // completely filled, it shows the last filled position. Once the buffers are full, it'll always
    // show the last position, and new samples will make the buffers shift, discarding the oldest
    // sample and storing the newest one on the last position.
    int current;

    // There are the variables from the original Pan-Tompkins algorithm.
    // The ones ending in _i correspond to values from the integrator.
    // The ones ending in _f correspond to values from the DC-block/low-pass/high-pass filtered signal.
    // The peak variables are peak candidates: signal values above the thresholds.
    // The threshold 1 variables are the threshold variables. If a signal sample is higher than this threshold, it's a peak.
    // The threshold 2 variables are half the threshold 1 ones. They're used for a back search when no peak is detected for too long.
    // The spk and npk variables are, respectively, running estimates of signal and noise peaks.
    dataType peak_i = 0, peak_f = 0, threshold_i1 = 0, threshold_i2 = 0, threshold_f1 = 0, threshold_f2 = 0, spk_i = 0, spk_f = 0, npk_i = 0, npk_f = 0;

    // qrs tells whether there was a detection or not.
    // regular tells whether the heart pace is regular or not.
    // prevRegular tells whether the heart beat was regular before the newest RR-interval was calculated.
    bool qrs, regular = true, prevRegular;

    // Initializing the RR averages
    for (i = 0; i < 8; i++)
    {
        rr1[i] = 0;
        rr2[i] = 0;
    }

    // The main loop where everything proposed in the paper happens. Ends when there are no more signal samples.
    do{
        // Test if the buffers are full.
        // If they are, shift them, discarding the oldest sample and adding the new one at the end.
        // Else, just put the newest sample in the next free position.
        // Update 'current' so that the program knows where's the newest sample.
        if (sample >= BUFFSIZE)
        {
            for (i = 0; i < BUFFSIZE - 1; i++)
            {
                signal[i] = signal[i+1];
                dcblock[i] = dcblock[i+1];
                lowpass[i] = lowpass[i+1];
                highpass[i] = highpass[i+1];
                derivative[i] = derivative[i+1];
                squared[i] = squared[i+1];
                integral[i] = integral[i+1];
                outputSignal[i] = outputSignal[i+1];
            }
            current = BUFFSIZE - 1;
        }
        else
        {
            current = sample;
        }
        signal[current] = (uint32_t)data;

        // If no sample was read, stop processing!
        if (signal[current] == (unsigned int)NOSAMPLE)
            break;
        sample++; // Update sample counter

        // DC Block filter
        // This was not proposed on the original paper.
        // It is not necessary and can be removed if your sensor or database has no DC noise.
        if (current >= 1)
            dcblock[current] = signal[current] - signal[current-1] + 0.995*dcblock[current-1];
        else
            dcblock[current] = 0;

        // Low Pass filter
        // Implemented as proposed by the original paper.
        // y(nT) = 2y(nT - T) - y(nT - 2T) + x(nT) - 2x(nT - 6T) + x(nT - 12T)
        // Can be removed if your signal was previously filtered, or replaced by a different filter.
        lowpass[current] = dcblock[current];
        if (current >= 1)
            lowpass[current] += 2*lowpass[current-1];
        if (current >= 2)
            lowpass[current] -= lowpass[current-2];
        if (current >= 6)
            lowpass[current] -= 2*dcblock[current-6];
        if (current >= 12)
            lowpass[current] += dcblock[current-12];

        // High Pass filter
        // Implemented as proposed by the original paper.
        // y(nT) = 32x(nT - 16T) - [y(nT - T) + x(nT) - x(nT - 32T)]
        // Can be removed if your signal was previously filtered, or replaced by a different filter.
        highpass[current] = -lowpass[current];
        if (current >= 1)
            highpass[current] -= highpass[current-1];
        if (current >= 16)
            highpass[current] += 32*lowpass[current-16];
        if (current >= 32)
            highpass[current] += lowpass[current-32];

        // Derivative filter
        // This is an alternative implementation, the central difference method.
        // f'(a) = [f(a+h) - f(a-h)]/2h
        // The original formula used by Pan-Tompkins was:
        // y(nT) = (1/8T)[-x(nT - 2T) - 2x(nT - T) + 2x(nT + T) + x(nT + 2T)]
        derivative[current] = highpass[current];
        if (current > 0)
            derivative[current] -= highpass[current-1];

        // This just squares the derivative, to get rid of negative values and emphasize high frequencies.
        // y(nT) = [x(nT)]^2.
        squared[current] = derivative[current]*derivative[current];

        // Moving-Window Integration
        // Implemented as proposed by the original paper.
        // y(nT) = (1/N)[x(nT - (N - 1)T) + x(nT - (N - 2)T) + ... x(nT)]
        // WINDOWSIZE, in samples, must be defined so that the window is ~150ms.

        integral[current] = 0;
        for (i = 0; i < WINDOWSIZE; i++)
        {
            if (current >= (dataType)i)
                integral[current] += highpass[current - i];
            else
                break;
        }
        integral[current] /= (dataType)i;

        qrs = false;

        // If the current signal is above one of the thresholds (integral or filtered signal), it's a peak candidate.
        if (integral[current] >= threshold_i1 || highpass[current] >= threshold_f1)
        {
            peak_i = integral[current];
            peak_f = highpass[current];
        }

        // If both the integral and the signal are above their thresholds, they're probably signal peaks.
        if ((integral[current] >= threshold_i1) && (highpass[current] >= threshold_f1))
        {
            // There's a 200ms latency. If the new peak respects this condition, we can keep testing.
            if (sample > lastQRS + FS/5)
            {
                // If it respects the 200ms latency, but it doesn't respect the 360ms latency, we check the slope.
                if (sample <= lastQRS + (long unsigned int)(0.36*FS))
                {
                    // The squared slope is "M" shaped. So we have to check nearby samples to make sure we're really looking
                    // at its peak value, rather than a low one.
                    currentSlope = 0;
                    for (j = current - 10; j <= current; j++)
                        if (squared[j] > currentSlope)
                            currentSlope = squared[j];

                    if (currentSlope <= (dataType)(lastSlope/2))
                    {
                        qrs = false;
                    }

                    else
                    {
                        spk_i = 0.125*peak_i + 0.875*spk_i;
                        threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
                        threshold_i2 = 0.5*threshold_i1;

                        spk_f = 0.125*peak_f + 0.875*spk_f;
                        threshold_f1 = npk_f + 0.25*(spk_f - npk_f);
                        threshold_f2 = 0.5*threshold_f1;

                        lastSlope = currentSlope;
                        qrs = true;
                    }
                }
                // If it was above both thresholds and respects both latency periods, it certainly is a R peak.
                else
                {
                    currentSlope = 0;
                    for (j = current - 10; j <= current; j++)
                        if (squared[j] > currentSlope)
                            currentSlope = squared[j];

                    spk_i = 0.125*peak_i + 0.875*spk_i;
                    threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
                    threshold_i2 = 0.5*threshold_i1;

                    spk_f = 0.125*peak_f + 0.875*spk_f;
                    threshold_f1 = npk_f + 0.25*(spk_f - npk_f);
                    threshold_f2 = 0.5*threshold_f1;

                    lastSlope = currentSlope;
                    qrs = true;
                }
            }
            // If the new peak doesn't respect the 200ms latency, it's noise. Update thresholds and move on to the next sample.
            else
            {
                peak_i = integral[current];
                npk_i = 0.125*peak_i + 0.875*npk_i;
                threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
                threshold_i2 = 0.5*threshold_i1;
                peak_f = highpass[current];
                npk_f = 0.125*peak_f + 0.875*npk_f;
                threshold_f1 = npk_f + 0.25*(spk_f - npk_f);
                threshold_f2 = 0.5*threshold_f1;
                qrs = false;
                outputSignal[current] = qrs;
                if (sample > DELAY + BUFFSIZE)
                    output(outputSignal[0]);
                continue;
            }

        }

        // If a R-peak was detected, the RR-averages must be updated.
        if (qrs)
        {
            // Add the newest RR-interval to the buffer and get the new average.
            rravg1 = 0;
            for (i = 0; i < 7; i++)
            {
                rr1[i] = rr1[i+1];
                rravg1 += rr1[i];
            }
            rr1[7] = sample - lastQRS;
            lastQRS = sample;
            rravg1 += rr1[7];
            rravg1 *= 0.125;

            // If the newly-discovered RR-average is normal, add it to the "normal" buffer and get the new "normal" average.
            // Update the "normal" beat parameters.
            if ( (rr1[7] >= rrlow) && (rr1[7] <= rrhigh) )
            {
                rravg2 = 0;
                for (i = 0; i < 7; i++)
                {
                    rr2[i] = rr2[i+1];
                    rravg2 += rr2[i];
                }
                rr2[7] = rr1[7];
                rravg2 += rr2[7];
                rravg2 *= 0.125;
                rrlow = 0.92*rravg2;
                rrhigh = 1.16*rravg2;
                rrmiss = 1.66*rravg2;
            }

            prevRegular = regular;
            if (rravg1 == rravg2)
            {
                regular = true;
            }
            // If the beat had been normal but turned odd, change the thresholds.
            else
            {
                regular = false;
                if (prevRegular)
                {
                    threshold_i1 /= 2;
                    threshold_f1 /= 2;
                }
            }
        }
        // If no R-peak was detected, it's important to check how long it's been since the last detection.
        else
        {
            // If no R-peak was detected for too long, use the lighter thresholds and do a back search.
            // However, the back search must respect the 200ms limit and the 360ms one (check the slope).
            if ((sample - lastQRS > (long unsigned int)rrmiss) && (sample > lastQRS + FS/5))
            {
                for (i = current - (sample - lastQRS) + FS/5; i < (long unsigned int)current; i++)
                {
                    if ( (integral[i] > threshold_i2) && (highpass[i] > threshold_f2))
                    {
                        currentSlope = 0;
                        for (j = i - 10; j <= i; j++)
                            if (squared[j] > currentSlope)
                                currentSlope = squared[j];

                        if ((currentSlope < (dataType)(lastSlope/2)) && (i + sample) < lastQRS + 0.36*lastQRS)
                        {
                            qrs = false;
                        }
                        else
                        {
                            peak_i = integral[i];
                            peak_f = highpass[i];
                            spk_i = 0.25*peak_i+ 0.75*spk_i;
                            spk_f = 0.25*peak_f + 0.75*spk_f;
                            threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
                            threshold_i2 = 0.5*threshold_i1;
                            lastSlope = currentSlope;
                            threshold_f1 = npk_f + 0.25*(spk_f - npk_f);
                            threshold_f2 = 0.5*threshold_f1;
                            // If a signal peak was detected on the back search, the RR attributes must be updated.
                            // This is the same thing done when a peak is detected on the first try.
                            //RR Average 1
                            rravg1 = 0;
                            for (j = 0; j < 7; j++)
                            {
                                rr1[j] = rr1[j+1];
                                rravg1 += rr1[j];
                            }
                            rr1[7] = sample - (current - i) - lastQRS;
                            qrs = true;
                            lastQRS = sample - (current - i);
                            rravg1 += rr1[7];
                            rravg1 *= 0.125;

                            //RR Average 2
                            if ( (rr1[7] >= rrlow) && (rr1[7] <= rrhigh) )
                            {
                                rravg2 = 0;
                                for (i = 0; i < 7; i++)
                                {
                                    rr2[i] = rr2[i+1];
                                    rravg2 += rr2[i];
                                }
                                rr2[7] = rr1[7];
                                rravg2 += rr2[7];
                                rravg2 *= 0.125;
                                rrlow = 0.92*rravg2;
                                rrhigh = 1.16*rravg2;
                                rrmiss = 1.66*rravg2;
                            }

                            prevRegular = regular;
                            if (rravg1 == rravg2)
                            {
                                regular = true;
                            }
                            else
                            {
                                regular = false;
                                if (prevRegular)
                                {
                                    threshold_i1 /= 2;
                                    threshold_f1 /= 2;
                                }
                            }

                            break;
                        }
                    }
                }

                if (qrs)
                {
                    outputSignal[current] = false;
                    outputSignal[i] = true;
                    if (sample > DELAY + BUFFSIZE)
                        output(outputSignal[0]);
                    continue;
                }
            }

            // Definitely no signal peak was detected.
            if (!qrs)
            {
                // If some kind of peak had been detected, then it's certainly a noise peak. Thresholds must be updated accordinly.
                if ((integral[current] >= threshold_i1) || (highpass[current] >= threshold_f1))
                {
                    peak_i = integral[current];
                    npk_i = 0.125*peak_i + 0.875*npk_i;
                    threshold_i1 = npk_i + 0.25*(spk_i - npk_i);
                    threshold_i2 = 0.5*threshold_i1;
                    peak_f = highpass[current];
                    npk_f = 0.125*peak_f + 0.875*npk_f;
                    threshold_f1 = npk_f + 0.25*(spk_f - npk_f);
                    threshold_f2 = 0.5*threshold_f1;
                }
            }
        }
        // The current implementation outputs '0' for every sample where no peak was detected,
        // and '1' for every sample where a peak was detected. It should be changed to fit
        // the desired application.
        // The 'if' accounts for the delay introduced by the filters: we only start outputting after the delay.
        // However, it updates a few samples back from the buffer. The reason is that if we update the detection
        // for the current sample, we might miss a peak that could've been found later by backsearching using
        // lighter thresholds. The final waveform output does match the original signal, though.
        outputSignal[current] = qrs;
        if (sample > DELAY + BUFFSIZE)
            output(outputSignal[0]);
    } while (signal[current] != (unsigned int)NOSAMPLE);

    // Output the last remaining samples on the buffer
    for (i = 1; i < BUFFSIZE; i++)
        output(outputSignal[i]);

}

void output(int out)
{
    printf("%d\n", out);
}
//******************************************************************************
// TI_ADS1293_SPI_DRDYB_PIN interrupt service routine
#pragma vector=TI_ADS1293_DRDYB_VECTOR
__interrupt void TI_ADS1293_DRDY_PORTx(void)
{
  TI_ADS1293_DRDYB_PxIFG &= ~TI_ADS1293_DRDYB_PIN;                             //  IFG cleared
  ADS1293_ADCDataReady = 1;                                                    // set flag
  __bic_SR_register_on_exit(LPM0_bits);                                        // Exit active CPU
}
//******************************************************************************
//EOF
