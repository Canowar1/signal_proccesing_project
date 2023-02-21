# signal_proccesing_project
It is the codes and report of the matlab project, which is the final assignment of the Signal and Systems course I took at the undergraduate level and includes the observations on the processing and modification of signals in the systems.
 Name Surname: Caner Budak Academic ID: 280070846121 A.M.: el22606
The report I prepared contains the explanations needed for the questions related to the Matlab assignment.
1.1 Lowpass Filters (Moving Average)
a) According to calculations in Matlab our coefficient can be found as; for N=2;
a=1, b=0.046221498602406 for N=4;
a=1, b=0.907557002795188 for N=10;
a=1, b=0.046221498602406
b) What can be observed by “freqz() command and their plots;
The code in matlab will plot the magnitude response and phase response of the three moving average filters with lengths N = 2, 4, and 10. The magnitude response shows the amount of attenuation or amplification of each frequency component of the input signal, while the phase response shows the phase shift introduced by the filter at each frequency.
Here's what we can observe from the plots:
Magnitude Response: As the length of the moving average filter increases, the magnitude response becomes flatter in the passband and exhibits a steeper roll-off in the stopband. This means that the filter is better able to remove high-frequency noise as the length increases. The cutoff frequency of the filter is equal to 1/N times the Nyquist frequency, so the cutoff frequency also decreases as the length increases.
Phase Response: The phase response of a moving average filter is linear, which means that the filter introduces a constant delay for all frequencies. The delay is equal to half the filter length in samples, which means that the delay increases as the length of the filter increases.
c) What can be observed in question three;
Firstly, I would like to mention I couldn’t find zpd.m. When the code runs, it will plot the pole-zero diagram of the moving average filters with lengths N = 2, 4, and 10. The zplane() function is used to plot the poles and zeros in the complex plane. The tf2zp() function is used to convert the filters to zero-pole representation, which is a way of expressing a filter's transfer function in terms of its zeros and poles.

 Here's what you can observe from the pole-zero diagrams:
The moving average filter with N = 2 has one zero at the origin and two poles at -1. This means that the filter has a linear phase response and introduces a delay of one sample.
The moving average filter with N = 4 has four zeros at the origin and four poles at -1. This means that the filter has a linear phase response and introduces a delay of two samples.
The moving average filter with N = 10 has ten zeros at the origin and ten poles at -1. This means that the filter has a linear phase response and introduces a delay of five samples.
In all cases, the filters have a simple pole-zero structure, with the zeros and poles arranged symmetrically about the origin. The phase response is linear, and the delay introduced by the filter is equal to half the filter length in samples.
1.2 Bandpass Filters
i) According to our calculations, transfer function coefficients are; (1,1), (-1,3600,-0600) and (0.7225,-0.7200)
ii) From the amplitude response plot, we can see that the filter has a low-pass characteristic, with a maximum magnitude of approximately 1 at 0 Hz (i.e., DC). The magnitude decreases rapidly as the frequency increases, reaching a cutoff frequency of around 0.5 times the sampling rate (or pi radians/sample in normalized frequency). This confirms the expectation from the pole-zero diagram that this filter is a low-pass filter.
iii)
iv) We can easily observe that the characteristics of filters changed a bit due to the change of pole locations.
v) We observe that when the input pulse has a period of T=50s, the system's response exhibits a decaying oscillatory behavior, which is expected given the location of the poles in the complex plane. As the period of the pulse is reduced to T=5s, the oscillatory behavior becomes more pronounced and the system takes longer to settle to its steady state. This is also expected, since a shorter period means that the input signal contains more high frequency components, which excite the system's natural oscillations more strongly.

 vi) We observe that the new filter has a wider pass band than the original filter, since the magnitude response does not drop off as quickly as before. This means that the new filter will allow more high frequency components to pass through, while attenuating the low frequency components to a lesser extent.
2.1 Analysis of Musical Signals
i) According to frequencies;
cello_note: 44.1 kHz
clarinet_note: 88.2 kHz
flute_note: 88.2 kHz
So the instrument which have a 44.1 kHz frequency is Cello
ii) By observing the waveforms, we can see that the signals are indeed periodic, which is expected for musical notes. To estimate the period, we can count the number of cycles in the plotted segment and divide by the segment duration.
iii) These frequencies are consistent with the expected pitch of the notes. We can also see that the relationship between the fundamental frequency and the period holds, since the period is the inverse of the fundamental frequency.
Looking at the spectrum shapes, we can see that all three signals have a similar pattern, with a strong peak at the fundamental frequency and weaker peaks at harmonics (integer multiples of the fundamental frequency). However, the relative magnitudes of the harmonics are different for each signal, which is what gives each note its unique timbre.
For the flute note, we can see that there are strong harmonics up to the 6th or 7th order, which contributes to its bright and piercing sound. For the clarinet note, we can see that there are fewer harmonics and they are weaker in magnitude compared to the fundamental frequency, which contributes to its smooth and mellow sound. For the cello note, we can see that there are many harmonics but they are weaker in magnitude and decay quickly, which contributes to its rich and warm sound.
Overall, the spectrum of each note gives us insight into its acoustic characteristics, and we can connect these characteristics to the number and magnitudes of the harmonics in the spectrum.
iv) From the plots, we can observe that the energy of the signal varies over time, reflecting the variations in the amplitude of the signal. We can also see that the energy plot is smoother than the original signal plot, since the energy is averaged over a window. The sharp peaks and dips in the original signal plot are smoothed out in the energy plot.

 v) From the spectrum plot, we can observe that the spectrum of the noisy signal has additional peaks or noise around the original harmonics, which are not present in the spectrum of the original cello note. This is due to the added Gaussian noise with zero mean and a standard deviation of 0.05, which introduces random variations in the signal. The noise also reduces the overall magnitude of the harmonics, making them harder to distinguish from the added noise.
vi) From the spectrum plot, we can see that the filter has effectively removed the added noise and restored the original harmonics of the cello note. The filter has also smoothed out some of the sharp peaks in the original spectrum, which is a typical behavior of low pass filters. The filtered signal sounds much cleaner and more similar to the original cello note, with a reduction in the hissing and buzzing sounds introduced by the added noise.
vii) From the spectrum plots, we can see that the band-pass filters have effectively isolated the desired harmonic frequencies, and attenuated the other frequency components. The filtered signals in the time domain show a clear oscillation at the frequency of the isolated harmonic, with minimal interference from other frequency components.
However, note that there are still some residual frequency components in the filtered signals, especially at the edges of the passbands. This is due to the transition band of the band-pass filter, and can be reduced by increasing the order of the filter or adjusting the filter parameters. Also, the filtered signals have lower energy and shorter duration than the original cello note, which is expected as we have removed most of the other frequency components.
2.2 Reconstruction of Musical Signals as Sums of Sinusoids
Since at this stage only calculations are included for most parts, I did not add it to the report step by step. Relevant data can be returned by running codes.
