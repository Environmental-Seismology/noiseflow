:py:mod:`noiseflow.cc.python.rfft`
==================================

.. py:module:: noiseflow.cc.python.rfft


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   noiseflow.cc.python.rfft.RFFTClass_python



Functions
~~~~~~~~~

.. autoapisummary::

   noiseflow.cc.python.rfft.time_norm_clip
   noiseflow.cc.python.rfft.time_norm_smooth



.. py:class:: RFFTClass_python(raw_data, dt, cc_len, cc_step, time_norm, clip_std, smooth_N, freq_norm, freqmin, freqmax, whiten_npad, smoothspect_N, flag, jobs)

   Bases: :py:obj:`object`

   .. py:method:: freq_norm(rfft_data)


   .. py:method:: process_chunk(chunk_start, chunk_end)


   .. py:method:: rfft(time_norm_data)


   .. py:method:: run()


   .. py:method:: time_norm(data)



.. py:function:: time_norm_clip(data, clip_std)


.. py:function:: time_norm_smooth(data, smooth_N)


