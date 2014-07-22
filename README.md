UPBL09a
=======

Software chunks for ID02 upgrade program

Dahu is a lightweight plugin based framework...
   ... technically a JSON-RPC server over Tango
* plugin can be class or can be generated from stateless function
* a plugin is executed within a job, each job lives in its own thread.
* plugins have empty constructors plus 4 methods (or more) 
 - setup allows to set the input parameters. It performs sanitization if needed
 - process does the taff
 - teardown sets the output and the logging and cleans up if needed
 - abort can be used to stop the processing if a plugin is a daemon.
* the job is responsible for serializing on disk the plugin input and output  
* jobs can be launched using the tango interface (or other ...)
* plugins have a single input and output.
* input/output are dictionaries where keys are strings and values are JSON-serializable (maybe pickle for numpy arrays?)
 
