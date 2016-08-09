# -*- coding: utf-8 -*-
from numpy import *
from openravepy import *
from openravepy.misc import SpaceSamplerExtra
import scipy
import copy

env=Environment()
env.SetViewer('qtcoin')
env.Load("../Turbina/env_mh12_0_16.xml")
