Dahu: Data Analysis RPC server over Tango
=========================================

Software chunks initially developped for **ID02 upgrade program** in 2012-2014.

Dahu is a lightweight plugin based framework...
   ... technically a JSON-RPC server over Tango
* plugins can be class or can be generated from any state-less function written in Python
* a plugin is executed within a job, each job lives in its own thread.
* plugins have empty constructors plus 4 methods (or more)
 - `setup` allows to set the input parameters. Sanitization is performed here.
 - `process` does the work.
 - `teardown` sets the output, the logging and cleans up (if needed).
 - `abort` can be used to stop the processing if a plugin is a daemon.
* the job is responsible for serializing on disk of the plugin's input and output
* jobs can be launched using the tango interface (or other ...)
* plugins have a single input and output, they are simple JSON-serializable dictionaries.

