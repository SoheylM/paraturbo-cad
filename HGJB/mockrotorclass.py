import cadquery as cq
import numpy as np
import os
import pickle
from math import sin, cos, pi, tan
import timeit
from enum import Enum, auto
import cq_warehouse.extensions

class Rotor():
    def __init__(self):
        self.cwf = os.getcwd().replace("\\", "/")
        