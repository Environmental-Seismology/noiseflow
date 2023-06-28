:py:mod:`noiseflow.client.client`
=================================

.. py:module:: noiseflow.client.client


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   noiseflow.client.client.downloader_https



Functions
~~~~~~~~~

.. autoapisummary::

   noiseflow.client.client.downloader



Attributes
~~~~~~~~~~

.. autoapisummary::

   noiseflow.client.client.cpu_cores


.. py:class:: downloader_https(url, show_info=True, resume=True, filename=None, num_threads=cpu_cores, timeout=10, chunk_size=1024 * 1000, header=None, proxies=None)

   :param url: link address
   :param filename: file name

   .. py:method:: check_url()

      check url support break-point resume and support multi-thread downloading


   .. py:method:: download()


   .. py:method:: download_by_piece(_range)


   .. py:method:: get_range(start=0)

      set download range
      eg: [(0, 1023), (1024, 2047), (2048, 3071) ...]


   .. py:method:: print()



.. py:function:: downloader(url, type='https', show_info=True, resume=True, filename=None, num_threads=cpu_cores, timeout=10, chunk_size=1024 * 1000, header=None, proxies=None)


.. py:data:: cpu_cores

   

