#!/usr/bin/env python3

"""
Reprocess a job using the current
"""

import os, json, logging, time
from argparse import ArgumentParser
import dahu.factory
from dahu.job import Job
from dahu.utils import get_workdir, NumpyEncoder
logging.basicConfig()

STATE_UNINITIALIZED = Job.STATE_UNINITIALIZED
STATE_RUNNING = Job.STATE_RUNNING
STATE_SUCCESS = Job.STATE_SUCCESS
STATE_FAILURE = Job.STATE_FAILURE
STATE_ABORTED = Job.STATE_ABORTED
STATE = Job.STATE


def parse():
    """
    Parse the command line a return the parsed arguments
    """
    # TODO
    parser = ArgumentParser(description='reProcess some data using the Dahu')
    parser.add_argument("-d", '--debug', dest='debug', action='store_true',
                        default=False, help='debug mode')
    parser.add_argument(dest='args', nargs='+', help='job inputs to be re-processed')
    args = parser.parse_args()
    if args.debug:
        logging.root.setLevel(logging.DEBUG)
    return args.args


def _run_(plugin, what):
    """
    run setup, process, teardown or abort ...

    @param plugin: plugin instance
    @param what: setup, process or teardown
    @return : execution code
    """

    methods = {"process":  plugin.DEFAULT_PROCESS,
               "setup":    plugin.DEFAULT_SET_UP,
               "teardown": plugin.DEFAULT_TEAR_DOWN,
               "abort":    plugin.DEFAULT_ABORT    }
    assert what in methods
    name = methods.get(what)
    if name in dir(plugin):
        method = plugin.__getattribute__(name)
        if "__call__" in dir(method):
            try:
                method()
            except Exception as error:
                import traceback
                err_msg = [traceback.format_exc(limit=10), ""
                           "Error %s while calling %s.%s" %
                           (error, plugin.__class__.__name__, what)]
                print(os.linesep.join(err_msg))
                return STATE_FAILURE
            else:
                return STATE_RUNNING
    else:
        print("No such method %s in class %s" % (what, plugin.__class__.__name__))
        return STATE_FAILURE


def process(args):
    """Process a set of arguments
    
    :param args: list of files to process
    """
    working_dir = get_workdir()
    for idx, fn in enumerate(args):
        print("Processing #%i: %s" % (idx, fn))
        if os.path.exists(fn):
            with open(fn, "r") as fp:
                dico = json.load(fp)
        else:
            logging.warning("No such file: %s" % fn)
            continue
        plugin_name = dico["plugin_name"]
        dico["job_id"] = idx
        plugin = dahu.factory.plugin_factory(plugin_name)
        #Drop some keys from input specific for reprocess
        for ignore in getattr(plugin, "REPROCESS_IGNORE", []):
            if ignore in dico:
                dico.pop(ignore, None)
        start_time = time.time()
        plugin.input = dico
        state = _run_(plugin, "setup")
        if state == STATE_RUNNING:
            status1 = _run_(plugin, "process")
            status2 = _run_(plugin, "teardown")
            if status1 == STATE_RUNNING and status2 == STATE_RUNNING:
                state = STATE_SUCCESS
            else:
                state = STATE_FAILURE
        print("Finished with state: %s" % state)
        plugin.output["job_runtime"] = time.time() - start_time
        result = json.dumps(plugin.output, indent=4, cls=NumpyEncoder)
        print(result)
        basename = os.path.join(working_dir, "%05i_%s" % (idx, plugin.get_name()))
        with open(basename + ".out", "w") as fp:
            fp.write(result)
        with open(basename + ".inp", "w") as fp:
            json.dump(plugin.input, fp, indent=4, cls=NumpyEncoder)


def main():
    args = parse()
    process(args)


if __name__ == "__main__":
    main()
