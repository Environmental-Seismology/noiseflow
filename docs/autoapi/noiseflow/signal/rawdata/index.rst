:py:mod:`noiseflow.signal.rawdata`
==================================

.. py:module:: noiseflow.signal.rawdata


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   noiseflow.signal.rawdata.RawData_Class




.. py:class:: RawData_Class(data, sampling_rate)

   Bases: :py:obj:`object`

   .. py:method:: bandpass(freqmin, freqmax, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False)


   .. py:method:: bandstop(freqmin, freqmax, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False)


   .. py:method:: decimate(sampling_rate_new, flag=False, flag_gap=None, threads=1, py=False)


   .. py:method:: detrend(type='linear', flag=False, flag_gap=None, threads=1, py=False)


   .. py:method:: highpass(freq, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False)


   .. py:method:: lowpass(freq, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False)


   .. py:method:: plot(channel_num=0, win_num=0, save=False, save_path=None, dpi=100)


   .. py:method:: save(save_path, format='npz', compression=False, h5_compression_format='gzip', h5_compression_opts=3)


   .. py:method:: taper(max_percentage=0.05, type='hann', side='both', flag=False, flag_gap=None, threads=1, py=False)



