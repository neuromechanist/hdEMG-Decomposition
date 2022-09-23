# hdEMG Decomposition pipeline

Do you ever wonder why EMG signals look so random?
You are not alone, many have looked at EMG signals as a random mess, and tried to filter it to get to linear envelopes that just show the overall muscle activity.

However, if you look closely, EMG signals recorded on the surface are a mixture of electrical activity (i.e., [action potentials](https://en.wikipedia.org/wiki/Action_potential))  from muscle [motor units](https://en.wikipedia.org/wiki/Motor_unit) that are transformed (or got convoluted) during their journey to the skin.

## Can we decompose EMG to read the motor unit activities?

Well, yes and no. **No**, because the widely available *bipolar* EMG sensors may not be enough to decompose the EMG signals to the motor unit activity. **Yes**, because if you have “enough” independent recordings from a muscle, you might be able to decompose those recordings to the motor unit activity using specialized bind-source separation (BSS) techniques. 

The “enough” recordings are usually collected using so-called **high-density** EMG arrays. These arrays have as low as 4 to as high as 64 electrodes that would record EMG signals in *monopolar* fashion, and the electrode distance is usually <1cm (hence it is called high-density). There are also some [research](https://ieeexplore.ieee.org/document/5409593) suggesting that the surface EMG decomposition can indeed pick up the motor units revealed by the needle EMG, which is the gold-standard to motor-unit recording.

## Awesome! Is there any open-source BSS available?

Well, the answer to this question is also yes and no 😅. There are several papers that offer pseudocode for the BSS methods, but almost none of them share the code itself.  The published BSS methods are mostly based on special cases of independent component analysis (ICA). Here is a nonexhaustive list of the methods:

1. Convolution Kernel Compensation (CKC), [Holobar and Zazula 2007](https://ieeexplore.ieee.org/document/4291854)
2. Fast ICA,  [ Meng et al 2022](https://doi.org/10.1016/j.bspc.2022.103615) [Chen and Zhou 2016](https://ieeexplore.ieee.org/document/7058391)
3. Convolutive sphering , [Negro et al 2016](https://iopscience.iop.org/article/10.1088/1741-2560/13/2/026027)

I only found the pipeline of one of the methods shared online. The CKC method and Convolute sphering are most likely behind a [paywall](https://demuse.feri.um.si/). The fast ICA approach implemented in *Meng 2022* is available from [physionet](https://physionet.org/content/hd-semg/1.0.0/) from the Hyser dataset ([Jiang et al 2021](https://ieeexplore.ieee.org/document/9438637)), and for convenience, from [GitHub](https://github.com/neuromechanist/fastICA_EMG_decomp) (Let’s call this the open-source code). 

## So, if there is one pipeline available, what is this code?

Good question! This pipeline (in its current form) is ~**not**~ a new method; it is rather a refined implementation of the open-source code. Upon reviewing the code, I found the open-source fast ICA is mostly customized toward the Hyper d

## What is the next step?

There are two routes, namely because there are two main processing platforms:

### Matlab development:

1. Implement the convolutive sphering method
2. Implement the pulse to noise ratio (PNR), [Holobar et al 2014](https://iopscience.iop.org/article/10.1088/1741-2560/11/1/016008/meta)
3. Using other ICA approaches, such as infomax
4. Develop processing pipelines for MU metrics, such as firing rate, conduction velocity, etc.

### Python/PyTorch development:

1. Move the pipeline to a PyTorch implementation. This has a benefit of using the GPU natively
2. Implement items 2 and 3 under Matlab development.

## How to run the code

The code is designed to take the high-density EMG as an array, or as a path pointing to a `mat` file structure as designed by OT Bioelecttronica [BioLab+](https://otbioelettronica.it/files/47/Software/129/OTBiolab-v1580.exe)(version 1.5.8.0 at the time of writing).

There are also two sample datasets that you can use for testing.
If you want to just run the code to get the motor units, their spike trains, and plot of the spike trains from the sample1 dataset, type in `run_decomposition` while you are in the root directory. It should work like a charm 😁. The `run_ICA` function, however, would take quite some time to complete.

If you want to bypass the `run_ICA`, and jump to the rest of the code, you just need to add two arguments to the function: `run_decomposition(‘load_ICA’,1,’save_flag’,0)` . I think the input arguments are self-explanatory.

Feel free to reach the documentation of each function and the input arguments of the `run_decomposition` function under the `initialize`` section to become familiar how to use the input arguments.

## Cool, how I can contribute?

Sure, any and every help is welcome.
**First**, please fork the repository, this helps with keeping track of the development. 

**Second**, check if the bug or improvement you are thinking about is not already in the issue’s list.

**Third**, you can start working on the bug fix based not he fork you made, and then have a pull request, so I can review and hopefully merge it to the repository. A more efficient way is to write an issue first, so we can discuss the approach and then work on the fix and go for a pull request.

## Are there any public datasets available to try?

Yes! To my knowledge, there are two great datasets available for upper limb with high-density EMG recordings:

1. Jiang, Xinyu, Xiangyu Liu, Jiahao Fan, Xinming Ye, Chenyun Dai, Edward A. Clancy, Metin Akay, and Wei Chen. 2021. “Open Access Dataset, Toolbox and Benchmark Processing Results of High-Density Surface Electromyogram Recordings.” *IEEE TNSRE*. https://doi.org/10.1109/TNSRE.2021.3082551.
2. Malešević, Nebojša, Alexander Olsson, Paulina Sager, Elin Andersson, Christian Cipriani, Marco Controzzi, Anders Björkman, and Christian Antfolk. 2021. “A Database of High-Density Surface Electromyogram Signals Comprising 65 Isometric Hand Gestures.” *Scientific Data* 8 (1): 63. http://dx.doi.org/10.1038/s41597-021-00843-9

There is also a great intramuscular EMG recordings from the upper limb:
Malesevic, Nebojsa, Anders Björkman, Gert S. Andersson, Ana Matran-Fernandez, Luca Citi, Christian Cipriani, and Christian Antfolk. 2020. “A Database of Multi-Channel Intramuscular Electromyogram Signals during Isometric Hand Muscles Contractions.” *Scientific Data* 7 (1): 10. http://dx.doi.org/10.1038/s41597-019-0335-8

## Further reading

There are many resources regarding the hd EMG:

There is a great tutorial on the underlying mechanisms and characteristics of the motor unit ([Del Vecchio et al 2022](http://dx.doi.org/10.1016/j.jelekin.2020.102426)).

There is also a paper on the best practices of hd EMG recording ([Gallina et al 2022](http://dx.doi.org/10.1016/j.jelekin.2022.102656)).

Also check the IEEE magazine article for the possible uses of the hdEMG ([Holobar and Farina 202](http://dx.doi.org/10.1109/MSP.2021.3057051))

## How to cite this repository

Please use the following to cite this repository:

Shirazi, S.Y. 2022: hdEMG-Decompostion (github.com/neuromechanist/hdEMG-Decomposition/tag/0.1), GitHub.