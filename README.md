# HADES
simple example of analysis

flowCount.C - is the main executing program. It should be compiled via ROOT with HYDRA on it.

flowHisto.C - is the initialization of all histograms that are needed for flowCount.C macro to prevent bad_alloc().

flowVar.C - is the initialization of all variables that are needed for flowCount.C macro to prevent bad_alloc().

Get* (GetRecenter.C, GetFlatSinCos.C, ...) - macros that generate calibration "*.cc" files for flowCount.C.

GetFourier.C - plot EP angle distribution and its Fourier cosine coefficients (till 6 harmonic).

flow3Vn.C - generates v[TYPE][HARMONIC][CENTRALITY][RAPIDITY][Pt] from obtained data. TYPE: raw,recentered, flattened.

FastFlow3.C - uses v[...] from flow3Vn.C to plot v1,v2 as functions of Pt/rapidity in centrality bins.

FuncPlot.C - file with class that draws histograms of (phi-Psi_EP) distribution.

ReadHistFunc.C - file with class that reads histograms from obtained data.

ReadRead.C - reading and drawing (phi-Psi_EP) distributions using FuncPlot.C and ReadHistFunc.C.

"Files/" directory contains all writed macros and other helpfull files. This is backup, just in case.
