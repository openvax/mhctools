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
from Queue import Queue

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

        with open(os.devnull, 'w') as devnull:
            stdout = redirect_stdout_file if redirect_stdout_file else devnull
            stderr = devnull if suppress_stderr else None
            self.process = Popen(args, stdout=stdout, stderr=stderr)

    def wait(self):
        ret_code = self.process.wait()
        logging.info(
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
    logging.info("%s took %0.4f seconds", cmd, elapsed_time)

def run_multiple_commands(multiple_args_lists, print_commands=True,
                          process_limit=0, **kwargs):
    """
    Run multiple shell commands in parallel.

    Parameters
    ----------

    multiple_args_lists : list of lists
        A collection of args lists to run on the shell.
        For example:
            [["ls", "-als"], ["rm", "-rf", "/dev"]]

    print_commands : bool
        Print shell commands before running them

   process_limit : int
        Limit the number of concurrent processes to this number. 0
        if there is no limit.
    """
    assert len(multiple_args_lists) > 0
    assert all(len(args) > 0 for args in multiple_args_lists)
    start_time = time.time()
    command_names = [args[0] for args in multiple_args_lists]
    processes = Queue(maxsize=process_limit)
    for args in multiple_args_lists:
        if print_commands:
            print(" ".join(args))
        p = AsyncProcess(args, **kwargs)
        if not processes.full():
            processes.put(p)
        else:
            processes.get().wait()
            processes.put(p)

    while not processes.empty():
        processes.get().wait()

    elapsed_time = time.time() - start_time
    logging.info("Ran %d commands (%s) in %0.4f seconds",
        len(multiple_args_lists),
        ",".join(command_names),
        elapsed_time
    )

def run_multiple_commands_redirect_stdout(
        multiple_args_dict,
        print_commands=True,
        process_limit=0,
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
        if there is no limit.
    """
    assert len(multiple_args_dict) > 0
    assert all(len(args) > 0 for args in multiple_args_dict.values())
    assert all(hasattr(f, 'name') for f in multiple_args_dict.keys())
    start_time = time.time()
    processes = Queue(maxsize=process_limit)
    for f, args in multiple_args_dict.iteritems():
        if print_commands:
            print(" ".join(args), ">", f.name)
        p = AsyncProcess(
            args,
            redirect_stdout_file=f,
            **kwargs)
        if not processes.full():
            processes.put(p)
        else:
            processes.get().wait()

    while not processes.empty():
        processes.get().wait()

    elapsed_time = time.time() - start_time
    logging.info("Ran %d commands in %0.4f seconds",
        len(multiple_args_dict),
        elapsed_time)
