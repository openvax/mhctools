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

import errno
import logging
import os
from subprocess import Popen, CalledProcessError
import time
from multiprocessing import cpu_count

logger = logging.getLogger(__name__)

# Retry settings for Popen when the OS is out of process slots.
_POPEN_MAX_RETRIES = 5
_POPEN_BASE_DELAY = 1.0  # seconds, doubled each retry

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
        self.suppress_stderr = suppress_stderr
        self.redirect_stdout_file = redirect_stdout_file
        self.process = None

    def start(self):
        with open(os.devnull, 'w') as devnull:
            stdout = (
                self.redirect_stdout_file if self.redirect_stdout_file
                else devnull)
            stderr = devnull if self.suppress_stderr else None
            delay = _POPEN_BASE_DELAY
            for attempt in range(_POPEN_MAX_RETRIES):
                try:
                    self.process = Popen(
                        self.args, stdout=stdout, stderr=stderr)
                    return
                except OSError as e:
                    if e.errno != errno.EAGAIN and not isinstance(
                            e, BlockingIOError):
                        raise
                    if attempt == _POPEN_MAX_RETRIES - 1:
                        raise
                    logger.warning(
                        "Resource temporarily unavailable starting "
                        "%s (attempt %d/%d), retrying in %.1fs",
                        self.cmd, attempt + 1, _POPEN_MAX_RETRIES,
                        delay)
                    time.sleep(delay)
                    delay *= 2

    def poll(self):
        """
        Peeks at whether the process is done or not, without
        waiting for it. Leaves exception handling and such to wait().
        """
        if self.process is None:
            self.start()
        return self.process.poll()

    def wait(self):
        if self.process is None:
            self.start()
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
    start_time = time.time()
    process = AsyncProcess(args, **kwargs)
    process.wait()
    elapsed_time = time.time() - start_time
    logger.info("%s took %0.4f seconds", args[0], elapsed_time)

def run_multiple_commands_redirect_stdout(
        multiple_args_dict,
        print_commands=True,
        process_limit=-1,
        polling_freq=0.5,
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
        if there is no limit, -1 to use max number of processors

    polling_freq : float
        Number of seconds between checking for done processes when
        at the concurrency limit.
    """
    assert len(multiple_args_dict) > 0
    assert all(len(args) > 0 for args in multiple_args_dict.values())
    assert all(hasattr(f, 'name') for f in multiple_args_dict.keys())
    if process_limit < 0:
        process_limit = cpu_count()
        logger.debug("Using %d processes", process_limit)

    start_time = time.time()
    active = []

    def _reap_finished():
        """Remove and wait() on completed processes."""
        still_running = []
        for p in active:
            if p.poll() is not None:
                p.wait()
            else:
                still_running.append(p)
        active[:] = still_running

    for f, args in multiple_args_dict.items():
        # Wait until we're below the limit before starting a new process
        while process_limit and len(active) >= process_limit:
            _reap_finished()
            if len(active) >= process_limit:
                time.sleep(polling_freq)

        p = AsyncProcess(
            args,
            redirect_stdout_file=f,
            **kwargs)
        p.start()
        if print_commands:
            handler = logging.FileHandler(p.redirect_stdout_file.name)
            handler.setLevel(logging.DEBUG)
            logger.addHandler(handler)
            logger.debug(" ".join(p.args))
            logger.removeHandler(handler)
        active.append(p)

    # Wait for remaining processes
    for p in active:
        p.wait()

    elapsed_time = time.time() - start_time
    logger.info(
        "Ran %d commands in %0.4f seconds",
        len(multiple_args_dict),
        elapsed_time)
