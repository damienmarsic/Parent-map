#!/usr/bin/env python
import pytest
pm=__import__('parent-map')

def test_check_file():
    assert not pm.check_file('%1@3$9)*+7','',0)

def test_check_int():
    assert pm.check_int(15,10,'*')

def test_match():
    assert pm.match('gggacgagt','aacgattttttt',4)==[(3,7,1)]

def test_refine():
    x,y=pm.refine(0,11,17,11,17,'VYSEPRPIGTRFLTRNL','VYSEPRPIGTRYLTRNL',[(0,11,0)],[])
    assert (x,y)==([(0, 11, 0), (11, 12, 'F', 11, 'Y'), (12, 17, 12)],[])

#test_refine()
pytest.main()
