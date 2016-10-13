# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import
import logging
import os
from subprocess import Popen, CalledProcessError
import time

# pylint: disable=import-error
from six.moves.queue import Queue


logger = logging.getLogger(__name__)


class AsyncProcess(object):
    """
    A thin wrapper around Popen which starts a process asynchronously,
    suppresses stdout printing, and raises an exception if the return code
    of wait() isn't 0
    """
    def __init__(
            self,
            args,
            suppress_stderr=False,
            redirect_stdout_file=None):
        assert len(args) > 0
        self.cmd = args[0]
        self.args = args
        self.redirect_stdout_file = redirect_stdout_file

        with open(os.devnull, 'w') as devnull:
            stdout = redirect_stdout_file if redirect_stdout_file else devnull
            stderr = devnull if suppress_stderr else None
            self.process = Popen(args, stdout=stdout, stderr=stderr)

    def poll(self):
        """
        Peeks at whether the process is done or not, without
        waiting for it. Leaves exception handling and such to wait().
        """
        return self.process.poll()

    def wait(self):
        ret_code = self.process.wait()
        logger.debug(
            "%s finished with return code %s",
            self.cmd,
            ret_code)
        if ret_code:
            raise CalledProcessError(ret_code, self.cmd)
        return ret_code

def run_command(args, **kwargs):
    """
    Given a list whose first element is a command name, followed by arguments,
    execute it and show timing info.
    """
    assert len(args) > 0
    cmd = args[0]
    start_time = time.time()
    process = AsyncProcess(args, **kwargs)
    process.wait()
    elapsed_time = time.time() - start_time
    logger.info("%s took %0.4f seconds", cmd, elapsed_time)

def run_multiple_commands_redirect_stdout(
        multiple_args_dict,
        print_commands=True,
        process_limit=0,
        polling_freq=1,
        **kwargs):
    """
    Run multiple shell commands in parallel, write each of their
    stdout output to files associated with each command.

    Parameters
    ----------
    multiple_args_dict : dict
        A dictionary whose keys are files and values are args list.
        Run each args list as a subprocess and write stdout to the
        corresponding file.

    print_commands : bool
        Print shell commands before running them.

    process_limit : int
        Limit the number of concurrent processes to this number. 0
        if there is no limit

    polling_freq : int
        Number of seconds between checking for done processes, if
        we have a process limit
    """
    assert len(multiple_args_dict) > 0
    assert all(len(args) > 0 for args in multiple_args_dict.values())
    assert all(hasattr(f, 'name') for f in multiple_args_dict.keys())
    start_time = time.time()
    processes = Queue(maxsize=process_limit)

    def add_to_queue(process):
        if print_commands:
            handler = logging.FileHandler(process.redirect_stdout_file.name)
            handler.setLevel(logging.DEBUG)
            logger.addHandler(handler)
            logger.debug(" ".join(process.args))
            logger.removeHandler(handler)
        processes.put(process)

    for f, args in multiple_args_dict.items():
        p = AsyncProcess(
            args,
            redirect_stdout_file=f,
            **kwargs)
        if not processes.full():
            add_to_queue(p)
        else:
            while processes.full():
                # Are there any done processes?
                to_remove = []
                for possibly_done in processes.queue:
                    if possibly_done.poll() is not None:
                        possibly_done.wait()
                        to_remove.append(possibly_done)
                # Remove them from the queue and stop checking
                if to_remove:
                    for process_to_remove in to_remove:
                        processes.queue.remove(process_to_remove)
                    break
                # Check again in a second if there weren't
                time.sleep(polling_freq)
            add_to_queue(p)

    # Wait for all the rest of the processes
    while not processes.empty():
        processes.get().wait()

    elapsed_time = time.time() - start_time
    logger.info(
        "Ran %d commands in %0.4f seconds",
        len(multiple_args_dict),
        elapsed_time)
