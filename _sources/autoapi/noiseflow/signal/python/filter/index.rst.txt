:py:mod:`noiseflow.signal.python.filter`
========================================

.. py:module:: noiseflow.signal.python.filter


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   noiseflow.signal.python.filter.bandpass
   noiseflow.signal.python.filter.bandpass_py
   noiseflow.signal.python.filter.bandstop
   noiseflow.signal.python.filter.bandstop_py
   noiseflow.signal.python.filter.highpass
   noiseflow.signal.python.filter.highpass_py
   noiseflow.signal.python.filter.lowpass
   noiseflow.signal.python.filter.lowpass_py



.. py:function:: bandpass(data, freqmin, freqmax, df, corners=4, zerophase=False)

   Butterworth-Bandpass Filter.
   Filter data from ``freqmin`` to ``freqmax`` using ``corners``
   corners.
   The filter uses :func:`scipy.signal.iirfilter` (for design)
   and :func:`scipy.signal.sosfilt` (for applying the filter).
   :type data: numpy.ndarray
   :param data: Data to filter.
   :param freqmin: Pass band low corner frequency.
   :param freqmax: Pass band high corner frequency.
   :param df: Sampling rate in Hz.
   :param corners: Filter corners / order.
   :param zerophase: If True, apply filter once forwards and once backwards.
       This results in twice the filter order but zero phase shift in
       the resulting filtered trace.
   :return: Filtered data.


.. py:function:: bandpass_py(data, freqmin, freqmax, df, corners=4, zerophase=True, flag=False, flag_gap=None)


.. py:function:: bandstop(data, freqmin, freqmax, df, corners=4, zerophase=False)

   Butterworth-Bandstop Filter.
   Filter data removing data between frequencies ``freqmin`` and ``freqmax``
   using ``corners`` corners.
   The filter uses :func:`scipy.signal.iirfilter` (for design)
   and :func:`scipy.signal.sosfilt` (for applying the filter).
   :type data: numpy.ndarray
   :param data: Data to filter.
   :param freqmin: Stop band low corner frequency.
   :param freqmax: Stop band high corner frequency.
   :param df: Sampling rate in Hz.
   :param corners: Filter corners / order.
   :param zerophase: If True, apply filter once forwards and once backwards.
       This results in twice the number of corners but zero phase shift in
       the resulting filtered trace.
   :return: Filtered data.


.. py:function:: bandstop_py(data, freqmin, freqmax, df, corners=4, zerophase=True, flag=False, flag_gap=None)


.. py:function:: highpass(data, freq, df, corners=4, zerophase=False)

   Butterworth-Highpass Filter.
   Filter data removing data below certain frequency ``freq`` using
   ``corners`` corners.
   The filter uses :func:`scipy.signal.iirfilter` (for design)
   and :func:`scipy.signal.sosfilt` (for applying the filter).
   :type data: numpy.ndarray
   :param data: Data to filter.
   :param freq: Filter corner frequency.
   :param df: Sampling rate in Hz.
   :param corners: Filter corners / order.
   :param zerophase: If True, apply filter once forwards and once backwards.
       This results in twice the number of corners but zero phase shift in
       the resulting filtered trace.
   :return: Filtered data.


.. py:function:: highpass_py(data, freq, df, corners=4, zerophase=True, flag=False, flag_gap=None)


.. py:function:: lowpass(data, freq, df, corners=4, zerophase=False)

   Butterworth-Lowpass Filter.
   Filter data removing data over certain frequency ``freq`` using ``corners``
   corners.
   The filter uses :func:`scipy.signal.iirfilter` (for design)
   and :func:`scipy.signal.sosfilt` (for applying the filter).
   :type data: numpy.ndarray
   :param data: Data to filter.
   :param freq: Filter corner frequency.
   :param df: Sampling rate in Hz.
   :param corners: Filter corners / order.
   :param zerophase: If True, apply filter once forwards and once backwards.
       This results in twice the number of corners but zero phase shift in
       the resulting filtered trace.
   :return: Filtered data.


.. py:function:: lowpass_py(data, freq, df, corners=4, zerophase=True, flag=False, flag_gap=None)


