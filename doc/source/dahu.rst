Dahu: online data analysis server
=================================

*Dahu* is a JSON-RPC server operated over Tango able to execute some Python code remotely for online data-analysis with a low latency.

The *dahu* server executes **jobs**:
------------------------------------

* Each job lives in its own thread (yes, thread, not process, it the plugin developper to ensure the work he is doing is GIL-compliant).
* Each job executes one plugin, provided by the plugin developper (i.e. the scientist)
* The job de/serialises JSON strings coming from/returning to Tango
* Jobs are executed asynchronously, the request for calculation is answered instantaneously with a *jobid*.
* The *jobid* can be used to poll the server for the status of the job or for manual synchronization (Tango can time-out!).
* When jobs are finished, the client is notified via Tango events about the status
* Results can be retrieved after the job has finished.

Jobs execute **plugin**:
------------------------

* Plugins are written in Python
* Plugins can be classes or simple function
* The input and output must be JSON-seriablisable as simple dictionnaries
* Plugins are dynamically loaded from Python modules
* Plugins can be profiled for performance analysis

Offline processing
------------------

All jobs can be run offline using the `dahu-reprocess` command line tool.
This tool is not multithreaded and plugins are directly run, it is intended for:

* offline developments
* re-processing some failed online processing.

Dahu is light !
---------------

Dahu is a small project started at ESRF in 2013 with less than 1000 lines of code.
It is used in production since then on a couple of beamlines.
With its FIFO scheduler, dahu is very fast (1µs, 0.3ms from Tango)

