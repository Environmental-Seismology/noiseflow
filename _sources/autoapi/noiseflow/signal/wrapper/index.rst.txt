:py:mod:`noiseflow.signal.wrapper`
==================================

.. py:module:: noiseflow.signal.wrapper


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   noiseflow.signal.wrapper.bandpass
   noiseflow.signal.wrapper.bandstop
   noiseflow.signal.wrapper.decimate
   noiseflow.signal.wrapper.detrend
   noiseflow.signal.wrapper.highpass
   noiseflow.signal.wrapper.lowpass
   noiseflow.signal.wrapper.taper



Attributes
~~~~~~~~~~

.. autoapisummary::

   noiseflow.signal.wrapper.NOISEFLOW_USE_CPP
   noiseflow.signal.wrapper.config
   noiseflow.signal.wrapper.config_path
   noiseflow.signal.wrapper.env


.. py:function:: bandpass(data, freqmin, freqmax, df, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False)


.. py:function:: bandstop(data, freqmin, freqmax, df, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False)


.. py:function:: decimate(data, df, df_new, flag=False, flag_gap=None, threads=1, py=False)


.. py:function:: detrend(data, type='linear', flag=False, flag_gap=None, threads=1, py=False)


.. py:function:: highpass(data, freq, df, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False)


.. py:function:: lowpass(data, freq, df, corners=4, zerophase=True, flag=False, flag_gap=None, threads=1, py=False)


.. py:function:: taper(data, max_percentage=0.05, type='hann', side='both', flag=False, flag_gap=None, threads=1, py=False)


.. py:data:: NOISEFLOW_USE_CPP
   :value: True

   

.. py:data:: config

   

.. py:data:: config_path

   

.. py:data:: env

   

