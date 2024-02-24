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

from functools import wraps

def ok_(b, msg=None):
    if not msg:
        assert b
    else:
        assert b, msg

def eq_(x, y, msg=None):
    if msg is None:
        assert x == y
    else:
        assert x == y, msg

def neq_(x, y, msg=None):
    if msg is None:
        assert x != y
    else:
        assert x != y, msg

def gt_(x, y, msg=None):
    if msg is None:
        assert x > y
    else:
        assert x > y, msg

def gte_(x, y, msg=None):
    if msg is None:
        assert x >= y
    else:
        assert x >= y, msg


def lt_(x, y, msg=None):
    if msg is None:
        assert x < y
    else:
        assert x < y, msg

def lte_(x, y, msg=None):
    if msg is None:
        assert x <= y
    else:
        assert x <= y, msg

class assert_raises:
    def __init__(self, exception_type):
        self.exception_type = exception_type

    def __enter__(self):
        pass

    def __exit__(self, type, value, traceback):
        if type is None:
            raise AssertionError("Expected exception %s not raised" % self.exception_type)
        if type != self.exception_type:
            raise AssertionError("Expected exception %s, got %s" % (self.exception_type, type))
        return True
    
def raises(exception_type):    
    def decorator(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            with assert_raises(exception_type):
                f(*args, **kwargs)
        return wrapper
    return decorator
