from nose.tools import eq_
from mhctools.cli.parsing_helpers import parse_int_list

def test_parse_int_list():
    # int by itself
    eq_(parse_int_list("1"), [1])
    # range of integers
    eq_(parse_int_list("1-3"), [1, 2, 3])
    # comma separated
    eq_(parse_int_list("9,10"), [9, 10])
    # comma separated with spaces
    eq_(parse_int_list("9, 10"), [9, 10])
    # comma separated and range together
    eq_(parse_int_list("9,10,12-13"), [9, 10, 12, 13])
