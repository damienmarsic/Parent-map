#!/usr/bin/env python
import pytest
pm=__import__('parent-map')

def test_check_int():
    assert pm.check_int(15,10,'*')

def test_match():
    assert pm.match('gggacgagt','aacgattttttt',4)==[(3,7,1)]

pytest.main()
