#! /usr/bin/python
# -*- coding: utf-8 -*-

# import logging
# logger = logging.getLogger(__name__)
from loguru import logger
import pytest
import os.path

path_to_script = os.path.dirname(os.path.abspath(__file__))

def inc(x):
    return x + 1



# @pytest.mark.interactive
# @pytest.mark.slow
def test_answer():
    assert inc(3) == 5
