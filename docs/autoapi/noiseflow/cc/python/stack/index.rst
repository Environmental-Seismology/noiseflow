:py:mod:`noiseflow.cc.python.stack`
===================================

.. py:module:: noiseflow.cc.python.stack


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   noiseflow.cc.python.stack.DOST
   noiseflow.cc.python.stack.StackClass_python



Functions
~~~~~~~~~

.. autoapisummary::

   noiseflow.cc.python.stack.adaptive_filter
   noiseflow.cc.python.stack.clusterstack
   noiseflow.cc.python.stack.nroot
   noiseflow.cc.python.stack.power2pad
   noiseflow.cc.python.stack.pws
   noiseflow.cc.python.stack.rms
   noiseflow.cc.python.stack.robust
   noiseflow.cc.python.stack.seisstack
   noiseflow.cc.python.stack.selective
   noiseflow.cc.python.stack.stack
   noiseflow.cc.python.stack.tfpws



.. py:class:: DOST(data)

   .. py:method:: dost(d)

      Discrete Orthonormal Stockwell Transform

      PARAMETERS:
      ---------------------
      d: array of time-series data (numpy.ndarray)

      RETURNS:
      ---------------------
      d_dost: array of DOST coefficients (numpy.ndarray)


   .. py:method:: dostbw(D)

      Calculate size of the DOST bandwidths

      PARAMETERS:
      ---------------------
      D: length of time-series data (int)

      RETURNS:
      ---------------------
      bw: list of DOST bandwidths


   .. py:method:: fourier(d)

      Normalized and centered fft

      PARAMETERS:
      ---------------------
      d: array of time series data (numpy.ndarray)

      RETURNS:
      ---------------------
      fftIn: array of frequency-domain data (numpy.ndarray)


   .. py:method:: idost(d)

      Inverse Discrete Orthonormal Stockwell Transform

      PARAMETERS:
      ---------------------
      d: array of DOST coefficients (numpy.ndarray)

      RETURNS:
      ---------------------
      d_idost: array of time-series data (numpy.ndarray)


   .. py:method:: ifourier(d)

      Normalized and centered ifft

      PARAMETERS:
      ---------------------
      d: array of frequency-domain data (numpy.ndarray)

      RETURNS:
      ---------------------
      ifftIn: array of time series data (numpy.ndarray)


   .. py:method:: pad(data)

      Zero pad data such that its length is a power of 2

      PARAMETERS:
      ---------------------
      data: array of time series data (numpy.ndarray)

      RETURNS:
      ---------------------
      data: zero-padded array of time series data (numpy.ndarray)



.. py:class:: StackClass_python(corr_data, stack_method, par, stack_all, stack_len, stack_step, pick, median_high, median_low, flag, jobs)

   Bases: :py:obj:`object`

   .. py:method:: check_ngood(nindex)


   .. py:method:: pick(pair_id)


   .. py:method:: process_chunk(chunk_start, chunk_end)


   .. py:method:: run()


   .. py:method:: stack(pair_id, nindex)



.. py:function:: adaptive_filter(d, g=1)

   the adaptive covariance filter to enhance coherent signals. Fellows the method of
   Nakata et al., 2015 (Appendix B)

   the filtered signal [x1] is given by x1 = ifft(P*x1(w)) where x1 is the ffted spectra
   and P is the filter. P is constructed by using the temporal covariance matrix.

   PARAMETERS:
   ----------------------
   d: numpy.ndarray contains the 2D traces of daily/hourly cross-correlation functions
   g: a positive number to adjust the filter harshness [default is 1]
   RETURNS:
   ----------------------
   newstack: numpy vector contains the stacked cross correlation function


.. py:function:: clusterstack(d, h=0.75, win=None, axis=0, normalize=True, plot=False)

   Performs stack after clustering. The data will be clustered into two groups.
   If the two centers of the clusters are similar (defined by corrcoef >= "t"), the original
   traces associated with both clusters will be used to produce the final linear stack, weighted by
   normalized SNR (phase clarity) of each cluster. Otherwise, the one with larger phase clarity
   (defined as max(abs(amplitudes))/rms(abs(amplitudes))) will be used to get the final stack.

   PARAMETERS:
   ---------------------
   d: N length array of time series data (numpy.ndarray)
   h: corrcoeff threshold to decide which group/cluster to use. Default 0.75.
   win: [start_index,end_index] used to compute the weight, instead of the entire trace. Default None.
           When None, use the entire trace.
   axis: which axis to stack. default 0.
   normalize: Normalize the traces before clustering. This will only influence the cluster.
           The final stack will be produced using the original data.
   plot: plot clustering results. default False.

   RETURNS:
   ---------------------
   newstack: final stack.


.. py:function:: nroot(d, p=2)

   this is nth-root stacking algorithm translated based on the matlab function
   from https://github.com/xtyangpsp/SeisStack (by Xiaotao Yang; follows the
   reference of Millet, F et al., 2019 JGR)

   Parameters:
   ------------
   d: numpy.ndarray contains the 2D cross correlation matrix
   p: np.int, nth root for the stacking. Default is 2.

   Returns:
   ------------
   newstack: np.ndarray, final stacked waveforms

   Written by Chengxin Jiang @ANU (May2020)


.. py:function:: power2pad(data)

   Zero pad data such that its length is a power of 2


.. py:function:: pws(d, p=2)

   Performs phase-weighted stack on array of time series. Modified on the noise function by Tim Climents.
   Follows methods of Schimmel and Paulssen, 1997.
   If s(t) is time series data (seismogram, or cross-correlation),
   S(t) = s(t) + i*H(s(t)), where H(s(t)) is Hilbert transform of s(t)
   S(t) = s(t) + i*H(s(t)) = A(t)*exp(i*phi(t)), where
   A(t) is envelope of s(t) and phi(t) is phase of s(t)
   Phase-weighted stack, g(t), is then:
   g(t) = 1/N sum j = 1:N s_j(t) * | 1/N sum k = 1:N exp[i * phi_k(t)]|^v
   where N is number of traces used, v is sharpness of phase-weighted stack

   PARAMETERS:
   ---------------------
   d: N length array of time series data (numpy.ndarray)
   p: exponent for phase stack (int). default is 2

   RETURNS:
   ---------------------
   newstack: Phase weighted stack of time series data (numpy.ndarray)


.. py:function:: rms(d)


.. py:function:: robust(d, epsilon=1e-05, maxstep=10, win=None, stat=False, ref=None)

   this is a robust stacking algorithm described in Pavlis and Vernon 2010. Generalized
   by Xiaotao Yang.

   PARAMETERS:
   ----------------------
   d: numpy.ndarray contains the 2D cross correlation matrix
   epsilon: residual threhold to quit the iteration (a small number). Default 1E-5
   maxstep: maximum iterations. default 10.
   win: [start_index,end_index] used to compute the weight, instead of the entire trace. Default None.
           When None, use the entire trace.
   ref: reference stack, with the same length as individual data. Default: None. Use median().
   RETURNS:
   ----------------------
   newstack: numpy vector contains the stacked cross correlation

   Written by Marine Denolle
   Modified by Xiaotao Yang


.. py:function:: seisstack(d, method, par=None)

   This is the same as stack(), to be compatible with old usage.


.. py:function:: selective(d, cc_min, epsilon=1e-05, maxstep=10, win=None, stat=False, ref=None)

   this is a selective stacking algorithm developed by Jared Bryan/Kurama Okubo.

   PARAMETERS:
   ----------------------
   d: numpy.ndarray contains the 2D cross correlation matrix
   epsilon: residual threhold to quit the iteration
   cc_min: numpy.float, threshold of correlation coefficient to be selected
   epsilon: residual threhold to quit the iteration (a small number). Default 1E-5
   maxstep: maximum iterations. default 10.
   win: [start_index,end_index] used to compute the weight, instead of the entire trace. Default None.
           When None, use the entire trace.
   ref: reference stack, with the same length as individual data. Default: None. Use mean().
   RETURNS:
   ----------------------
   newstack: numpy vector contains the stacked cross correlation
   nstep: np.int, total number of iterations for the stacking

   Originally ritten by Marine Denolle
   Modified by Chengxin Jiang @Harvard (Oct2020)


.. py:function:: stack(d, method, par=None)

   this is a wrapper for calling individual stacking functions.
   d: data. 2-d array
   method: stacking method, one of "linear","pws","robust","acf","nroot","selective",
           "cluster"
   par: dictionary containing all parameters for each stacking method. defaults will
       be used if not specified.

   RETURNS:
   ds: stacked data, which may be a list depending on the method.


.. py:function:: tfpws(d, p=2, axis=0)

   Performs time-frequency domain phase-weighted stack on array of time series.

   $C_{ps} = |(\sum{S*e^{i2\pi}/|S|})/M|^p$, where $C_{ps}$ is the phase weight. Then
   $S_{pws} = C_{ps}*S_{ls}$, where $S_{ls}$ is the S transform of the linea stack
   of the whole data.

   Reference:
   Schimmel, M., Stutzmann, E., & Gallart, J. (2011). Using instantaneous phase
   coherence for signal extraction from ambient noise data at a local to a
   global scale. Geophysical Journal International, 184(1), 494–506.
   https://doi.org/10.1111/j.1365-246X.2010.04861.x

   PARAMETERS:
   ---------------------
   d: N length array of time series data (numpy.ndarray)
   p: exponent for phase stack (int). default is 2
   axis: axis to stack, default is 0.

   RETURNS:
   ---------------------
   newstack: Phase weighted stack of time series data (numpy.ndarray)


